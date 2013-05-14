#' Prettifies the table returned by population stats for the pipeline
#'
#' This function updates the rownames to retain only the number of nodes
#' specified. For example, suppose the full path to the node is
#' /singlet/viable/Lymph/CD3/CD8/IFNg. By default, this is updated to IFNg.
#' If instead we specify \code{nodes = 2}, then the corresponding row is
#' CD8/IFNg.
#' 
#' @param population_stats a data frame containing the population summaries,
#' which is returned from the \code{getPopStats} function in the
#' \code{flowWorkspace} package.
#' @param nodes numeric. How many nodes should be preserved in the prettified
#' population statistics data frame returned? Default: 1
#' @return a data frame containing the population statistics but with prettified
#' row names
pretty_popstats <- function(population_stats, nodes = 1) {
  node_names <- sapply(strsplit(rownames(population_stats), split = "/"), tail, n = nodes)
  rownames(population_stats) <- sapply(node_names, paste, collapse = "/")

  population_stats
}

#' Removes commas from a numeric stored as a string.
#'
#' In the summary for the HVTN065 manual gates, the cellular population counts
#' are stored as character strings with commas denoting the thousands place. For
#' instance, 3000 is stored as "3,000". We remove the commas from the strings
#' and convert to the resulting string to a numeric.
#'
#' @param x character string with the nuisance commas
#' @return the updated numeric value
no_commas <- function(x) {
  as.numeric(gsub(",", "", x))
}

#' Constructs a vector of all the combinations of A & B & C
#'
#' The \code{permutations} function is from the \code{gtools} package on CRAN.
#' @param markers character vector of marker names
#' @return vector containing all combinations of the markers
#' @examples
#' polyfunction_nodes(c("IFNg", "IL2", "TNFa", "GzB", "CD57"))
polyfunction_nodes <- function(markers) {
  num_markers <- length(markers)
  plusminus_list <- permutations(n = 2, r = num_markers, c("+", "-"),
                                 repeats = TRUE)
  apply(plusminus_list, 1, function(plusminus_row) {
    paste0(markers, plusminus_row, collapse = "")
  })
}


#' Classification summary for population statistics
#'
#' For a data.frame of population statistics, we conduct a classification summary
#' of the the stimulation groups using the given treatment information.
#'
#' Per Greg, we apply the following classification study:
#' 1. Remove placebos before subsetting treatment group train classifier.
#' 2. Predict placebos (should expect poor classification accuracy because there
#'    should be no separation in the placebos)
#' 3. Predict test data set (should expect good results)
#' 
#' Per Greg:
#' "We also want to do this paired, using the difference in classification
#' probabilities for two samples from the same subject. i.e.
#' d = Pr(sample 1 from subject 1 = post-vaccine) -
#' Pr(sample 2 from subject 1 = post-vaccine).
#' If d > threshold, then classify sample 1 as post-vaccine and sample 2 as
#' pre-vaccine, otherwise if d < threshold classify sample 1 as pre-vaccine and
#' sample 2 as post-vaccine, otherwise mark them as unclassifiable."
#'
#' Hence, we compute classification accuracies by pairing the visits by patient
#' and then computing the proportion of correctly classified patients. Otherwise,
#' we calculate the proportion of visits that are correctly classified.
#'
#' @param popstats a data.frame containing population statistics from
#' \code{getPopStats}
#' @param treatment_info a data.frame containing a lookup of \code{PTID} and
#' placebo/treatment information
#' @param pdata an object returned from \code{pData} from the \code{GatingSet}
#' @param stimulation the stimulation group
#' @param train_pct a numeric value determining the percentage of treated
#' patients used as training data and the remaining patients as test data
#' @param prob_threshold a numeric value above which the difference in
#' classification probabilities for two samples from the same subject indicates
#' that the first sample is classified as post-vaccine. Ignored if \code{paired}
#' is \code{FALSE}. See details.
#' @param ... Additional arguments passed to \code{\link{glmnet}}.
#' @return list containing classification results
classification_summary <- function(popstats, treatment_info, pdata,
                                   stimulation = "GAG-1-PTEG", train_pct = 0.6,
                                   prob_threshold = 0, ...) {
  m_popstats <- melt(popstats)
  colnames(m_popstats) <- c("Marker", "Sample", "Proportion")
  m_popstats$Marker <- as.character(m_popstats$Marker)

  m_popstats <- plyr:::join(m_popstats, pdata, by = "Sample")
  m_popstats$VISITNO <- factor(m_popstats$VISITNO)
  m_popstats$PTID <- factor(m_popstats$PTID)

  # We stored the 'negctrl' with sample numbers appended to the strings so that
  # plotGate could identify unique samples. Here, we strip the sample numbers to
  # summarize the negative controls as a whole.
  m_popstats$Stimulation <- as.character(m_popstats$Stimulation)
  m_popstats$Stimulation <- replace(m_popstats$Stimulation,
                                    grep("^negctrl", m_popstats$Stimulation), "negctrl")

  # Next, we include only the negative controls and the specified stimulation groups.
  m_popstats <- subset(m_popstats, Stimulation %in% c("negctrl", stimulation))
  m_popstats$Stimulation <- factor(m_popstats$Stimulation)

  # Here, we summarize the population statistics within Stimulation group.
  # Effectively, this averages the two negative-control proportions for each
  # marker.
  m_popstats <- ddply(m_popstats, .(PTID, VISITNO, Stimulation, Marker),
                      summarize, Proportion = mean(Proportion))

  # Next, we normalize the population proportions for the stimulated samples to
  # adjust for the background (negative controls) by calculating the difference of
  # the proportions for the stimulated samples and the negative controls.
  m_popstats <- ddply(m_popstats, .(PTID, VISITNO, Marker), summarize,
                    diff_Proportion = diff(Proportion))

  # Converts the melted data.frame to a wider format to continue the classification study.
  m_popstats <- dcast(m_popstats, PTID + VISITNO ~ Marker, value.var = "diff_Proportion")
  m_popstats <- plyr:::join(m_popstats, treatment_info)
  m_popstats$PTID <- as.character(m_popstats$PTID)

  placebo_data <- subset(m_popstats, Treatment == "Placebo", select = -Treatment)
  treatment_data <- subset(m_popstats, Treatment == "Treatment", select = -Treatment)

  # Partitions GAG data for classification study
  treated_patients <- unique(treatment_data$PTID)
  num_treated_patients <- length(treated_patients)
  patients_train <- sample.int(num_treated_patients, train_pct * num_treated_patients)
  
  train_data <- subset(treatment_data, PTID %in% treated_patients[patients_train])
  test_data <- subset(treatment_data, PTID %in% treated_patients[-patients_train])
  
  train_x <- as.matrix(subset(train_data, select = -c(PTID, VISITNO)))
  train_y <- train_data$VISITNO
  
  test_x <- as.matrix(subset(test_data, select = -c(PTID, VISITNO)))
  test_y <- test_data$VISITNO
  
  placebo_x <- as.matrix(subset(placebo_data, select = -c(PTID, VISITNO)))
  placebo_y <- placebo_data$VISITNO
  
  # Trains the 'glmnet' classifier using cross-validation.
  glmnet_cv <- cv.glmnet(x = train_x, y = train_y, family = "binomial", ...)

  # Computes classification accuracies in two different ways. If the visits are
  # paired by patient, then we compute the proportion of correctly classified
  # patients. Otherwise, we calculate the proportion of visits that are correctly
  # classified.
  predictions_treated <- as.vector(predict(glmnet_cv, test_x, s = "lambda.min",
                                           type = "response"))
  predictions_placebo <- as.vector(predict(glmnet_cv, placebo_x,
                                           s = "lambda.min", type = "response"))
   
  # If the difference in classification probabilities exceeds the probability
  # threshold, we assign the first sample as visit 2 and the second as visit 12.
  # Otherwise, we assign the first sample as visit 12 and the second as visit 2.
  correct_patients <- tapply(seq_along(test_data$PTID), test_data$PTID, function(i) {
    if (diff(predictions_treated[i]) > prob_threshold) {
      classification <- c("2", "12")
    } else {
      classification <- c("12", "2")
    }
    all(classification == test_y[i])
  })
  correct_placebo <- tapply(seq_along(placebo_data$PTID), placebo_data$PTID, function(i) {
   if (diff(predictions_placebo[i]) > prob_threshold) {
      classification <- c("2", "12")
    } else {
      classification <- c("12", "2")
    }
    all(classification == placebo_y[i])
  })

  accuracy_treated <- mean(correct_patients)
  accuracy_placebo <- mean(correct_placebo)
  
  # Determines which features should be kept for classification
  # If present, we manually remove the "(Intercept)"
  # In the case that all markers are removed and only the intercept remains, we
  # add the "(Intercept)" back, for clarity.
  coef_glmnet <- coef(glmnet_cv)
  markers_kept <- rownames(coef_glmnet)[as.vector(coef_glmnet) != 0]
  markers_kept <- markers_kept[!grepl("(Intercept)", markers_kept)]
  if (length(markers_kept) == 0) {
    markers_kept <- "(Intercept)"
  }

  # We return the classification accuracies, the classification probabilities and
  # markers kept using 'glmnet', the corresponding test data sets for further
  # analysis (e.g., ROC curves), and the markers kept by 'glmnet'.
  list(accuracy = list(treatment = accuracy_treated, placebo = accuracy_placebo),
       classification_probs = list(treated = predictions_treated,
         placebo = predictions_placebo),
       markers = markers_kept,
       test_data = list(treated = test_data, placebo = placebo_data))
}

#' Scales a vector of data using the Huber robust estimator for mean and
#' standard deviation
#'
#' This function is an analog to \code{\link{scale}} but using Huber robust
#' estimators instead of the usual sample mean and standard deviation.
#'
#' @param x numeric vector
#' @param center logical value. Should \code{x} be centered?
#' @param scale logical value. Should \code{x} be scaled?
#' @return numeric vector containing the scaled data
scale_huber <- function(x, center = TRUE, scale = TRUE) {
  x <- as.vector(x)
  huber_x <- huber(x)

  # If 'center' is set to TRUE, we center 'x' by the Huber robust location
  # estimator.
  center_x <- FALSE
  if (center) {
    center_x <- huber_x$mu
  }

  # If 'scale' is set to TRUE, we scale 'x' by the Huber robust standard
  # deviation estimator.
  scale_x <- FALSE
  if (scale) {
    scale_x <- huber_x$s
  }

  as.vector(base:::scale(x, center = center_x, scale = scale_x))
}

#' Centers a vector of data using the mode of the kernel density estimate
#'
#' @param x numeric vector
#' @param ... additional arguments passed to \code{\link{density}}
#' @return numeric vector containing the centered data
center_mode <- function(x, ...) {
  x <- as.vector(x)
  density_x <- density(x)
  mode <- density_x$x[which.max(density_x$y)]

  as.vector(scale(x, center = mode, scale = FALSE))
}

#' Extracts legend from ggplot2 object
#'
#' Code from Hadley Wickham:
#' https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#' 
#' @param a.gplot a \code{ggplot2} object
#' @return a \code{ggplot2} legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Look up channel name for cytokine
#'
#' @param x cytokine marker name
#' @return the corresponding channel name
cytokine2channel <- function(x) {
  switch(x,
         TNFa = "Alexa 680-A",
         IFNg = "PE Cy7-A",
         IL2 = "PE Green laser-A"
  )
}

#' Standardizes a cytokine data with respect to a reference stimulation group
#' 
#' @param x melted data.frame containing cytokine data for each stimulation
#' group (column labeled 'Stim')
#' @param ref_stimulation reference stimulation group
#' @return data.frame with cytokine data standardized with respect to the
#' reference stimulation group
standardize_cytokines <- function(x) {
  # Standardizes the cytokine samples within stimulation group with respect to
  # the reference stimulation group.
  # First, centers the values by the mode of the kernel density estimate for
  # the stimulation group.
  x <- tapply(x$value, x$Stim, center_mode)

  # Scales the stimulation groups by their Huber estimator of the
  # standard deviation to put all stimulations on the same scale.
  x <- lapply(x, scale_huber, center = FALSE)

  x <- melt(x)
  colnames(x) <- c("value", "Stim")
  x
}

#' Reads the specified cytokine from a gating set
#' @param gs gating set object
#' @param PTID patient ID
#' @param VISITNO patient visit number
#' @param tcells which T-cell type?
#' @param cytokine which cytokine?
#' @return a data.frame with the appropriate cytokine data
read_cytokines <- function(gs, PTID, VISITNO = c("2", "12"), tcells = c("cd4", "cd8"),
                           cytokine = c("TNFa", "IFNg", "IL2")) {
  pData_gs <- pData(gs)
  pData_gs <- pData_gs[pData_gs$PTID == PTID & pData_gs$VISITNO == VISITNO, ]

  x <- lapply(seq_along(pData_gs$name), function(i) {
    flow_set <- getData(gs[pData_gs$name[i]], tcells)
    # In the case that 0 or 1 cells are present, we return NULL. The case of 1
    # cell is ignored because a density cannot be estimated from it. In fact,
    # density(...) throws an error in this case.
    if (nrow(flow_set[[1]]) >= 2) {
      cbind("Stim" = as.character(pData_gs$Stim[i]), "value" = exprs(flow_set[[1]])[, cytokine2channel(cytokine)])
    } else {
      NULL
    }
  })
  x <- do.call(rbind, x)
  x <- data.frame(x, stringsAsFactors = FALSE)
  x$value <- as.numeric(x$value)
  x
}

#' Constructs the kernel density estimates for each stimulation group for a
#' given set of cytokine data
#' 
#' @param x melted data.frame containing cytokine data for each stimulation
#' group (column labeled 'Stim')
#' @param ... additional arguments passed to \code{\link{density}}
#' @return data.frame with the cytokine densities for each stimulation group
density_cytokines <- function(x, ...) {
  density_x <- tapply(x$value, x$Stim, density, ...)
  density_x <- lapply(seq_along(density_x), function(i) {
    cbind(Stim = names(density_x)[i], x = density_x[[i]]$x,
          y = density_x[[i]]$y)
  })
  do.call(rbind, density_x)
}

#' Constructs the derivatives of the kernel density estimates for each
#' stimulation group for a given set of cytokine data
#'
#' For guidance on selecting the bandwidth, see this post:
#' http://stats.stackexchange.com/questions/33918/is-there-an-optimal-bandwidth-for-a-kernel-density-estimator-of-derivatives
#' 
#' @param x melted data.frame containing cytokine data for each stimulation
#' group (column labeled 'Stim')
#' @param adjust a numeric weight on the automatic bandwidth
#' @param collapse If \code{TRUE}, the values are collapsed from all stimulation
#' groups before estimating the derivative of the density function
#' @param ... additional arguments passed to \code{\link{density}}
#' @return data.frame with the cytokine densities for each stimulation group
deriv_cytokines <- function(x, deriv = 1, adjust = 1, collapse = FALSE) {
  require('feature')
  require('ks')
  if (collapse) {
    deriv_x <- drvkde(x = x$value, drv = deriv, bandwidth = adjust * hpi(x$value))
    deriv_x <- cbind(x = deriv_x$x.grid[[1]], y = deriv_x$est)
  } else {
    deriv_x <- tapply(seq_along(x$Stim), x$Stim, function(i) {
      drvkde(x = x$value[i], drv = deriv, bandwidth = adjust * hpi(x$value[i]))
    })
    deriv_x <- lapply(seq_along(deriv_x), function(i) {
      cbind(Stim = names(deriv_x)[i], x = deriv_x[[i]]$x.grid[[1]],
            y = deriv_x[[i]]$est)
    })
    deriv_x <- do.call(rbind, deriv_x)
  }
  deriv_x
}

#' Constructs a cutpoint by collapsing the values of a stimulation group and
#' using the first derivative of the kernel density estimate to establish a
#' cutpoint
#'
#' @param x melted data.frame containing cytokine data for each stimulation
#' group (column labeled 'Stim')
#' groups before estimating the derivative of the density function
#' @param ... additional arguments passed to \code{\link{deriv_cytokine}}
#' @return the cutpoint along the x-axis
cytokine_cutpoint <- function(x, tol = 0.001, ...) {
  deriv_out <- deriv_cytokines(x = x, deriv = 1, collapse = TRUE, ...)
  x <- deriv_TNFa[, 1]
  y <- deriv_TNFa[, 2]
  lowest_valley <- x[which.min(y)]
  cutpoint <- x[which(x > lowest_valley & abs(y) < tol)[1]]
  cutpoint
}
