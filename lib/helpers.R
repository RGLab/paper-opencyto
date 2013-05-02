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
#' Hence, if \code{paired} is \code{TRUE}, we compute classification accuracies
#' by pairing the visits by patient and then computing the proportion of
#' correctly classified patients. Otherwise, we calculate the proportion of
#' visits that are correctly classified.
#'
#' @param popstats a data.frame containing population statistics from
#' \code{getPopStats}
#' @param treatment_info a data.frame containing a lookup of \code{PTID} and
#' placebo/treatment information
#' @param pdata an object returned from \code{pData} from the \code{GatingSet}
#' @param train_pct a numeric value determining the percentage of treated
#' patients used as training data and the remaining patients as test data
#' @param paired logical value that determines if the classification accuracies
#' should be calculated by pairing visits for a patient
#' @param prob_threshold a numeric value above which the difference in
#' classification probabilities for two samples from the same subject indicates
#' that the first sample is classified as post-vaccine. Ignored if \code{paired}
#' is \code{FALSE}. See details.
#' @return TODO
classification_summary <- function(popstats, treatment_info, pdata,
                                   train_pct = 0.6, paired = FALSE,
                                   prob_threshold = 0) {
  m_popstats <- melt(popstats)
  colnames(m_popstats) <- c("Marker", "Sample", "Proportion")
  m_popstats$Marker <- as.character(m_popstats$Marker)

  m_popstats <- plyr:::join(m_popstats, pdata, by = "Sample")
  m_popstats$VISITNO <- factor(m_popstats$VISITNO)
  m_popstats$PTID <- factor(m_popstats$PTID)

  # We stored the 'negctrl' with sample numbers appended to the strings so that
  # plotGate could identify unique samples. Here, we strip the sample numbers to
  # summarize the negative controls as a whole.
  m_popstats$Stimulation <- with(m_popstats, replace(Stimulation,
                                                   grep("^negctrl", Stimulation),
                                                   "negctrl"))

  m_popstats <- ddply(m_popstats, .(PTID, VISITNO, Stimulation, Marker),
                      summarize, Proportion = mean(Proportion))

  ENV_data <- subset(m_popstats, Stimulation != "GAG-1-PTEG")
  GAG_data <- subset(m_popstats, Stimulation != "ENV-1-PTEG")

  # Next, we normalize the population proportions for the stimulated samples to
  # adjust for the background (negative controls) by calculating the difference of
  # the proportions for the stimulated samples and the negative controls.
  ENV_data <- ddply(ENV_data, .(PTID, VISITNO, Marker), summarize,
                    diff_Proportion = diff(Proportion))
  GAG_data <- ddply(GAG_data, .(PTID, VISITNO, Marker), summarize,
                    diff_Proportion = diff(Proportion))

  # Converts the melted data.frame to a wider format to continue the classification study.
  ENV_data <- dcast(ENV_data, PTID + VISITNO ~ Marker, value.var = "diff_Proportion")
  ENV_data <- plyr:::join(ENV_data, treatment_info)
  ENV_data$PTID <- as.character(ENV_data$PTID)

  GAG_data <- dcast(GAG_data, PTID + VISITNO ~ Marker, value.var = "diff_Proportion")
  GAG_data <- plyr:::join(GAG_data, treatment_info)
  GAG_data$PTID <- as.character(GAG_data$PTID)

  ENV_placebo_data <- subset(ENV_data, Treatment == "Placebo", select = -Treatment)
  ENV_treatment_data <- subset(ENV_data, Treatment == "Treatment", select = -Treatment)

  GAG_placebo_data <- subset(GAG_data, Treatment == "Placebo", select = -Treatment)
  GAG_treatment_data <- subset(GAG_data, Treatment == "Treatment", select = -Treatment)

  # Partitions ENV data for classification study
  ENV_treated_patients <- unique(ENV_treatment_data$PTID)
  num_ENV_treated_patients <- length(ENV_treated_patients)
  ENV_patients_train <- sample.int(num_ENV_treated_patients,
                                   train_pct * num_ENV_treated_patients)

  ENV_train_data <- subset(ENV_treatment_data,
                           PTID %in% ENV_treated_patients[ENV_patients_train])
  ENV_test_data <- subset(ENV_treatment_data,
                          PTID %in% ENV_treated_patients[-ENV_patients_train])
  
  ENV_train_x <- as.matrix(subset(ENV_train_data, select = -c(PTID, VISITNO)))
  ENV_train_y <- ENV_train_data$VISITNO
  
  ENV_test_x <- as.matrix(subset(ENV_test_data, select = -c(PTID, VISITNO)))
  ENV_test_y <- ENV_test_data$VISITNO
  
  ENV_placebo_x <- as.matrix(subset(ENV_placebo_data, select = -c(PTID, VISITNO)))
  ENV_placebo_y <- ENV_placebo_data$VISITNO
  
  # Partitions GAG data for classification study
  GAG_treated_patients <- unique(GAG_treatment_data$PTID)
  num_GAG_treated_patients <- length(GAG_treated_patients)
  GAG_patients_train <- sample.int(num_GAG_treated_patients,
                                   train_pct * num_GAG_treated_patients)
  
  GAG_train_data <- subset(GAG_treatment_data, PTID %in% GAG_treated_patients[GAG_patients_train])
  GAG_test_data <- subset(GAG_treatment_data,
                          PTID %in% GAG_treated_patients[-GAG_patients_train])
  
  GAG_train_x <- as.matrix(subset(GAG_train_data, select = -c(PTID, VISITNO)))
  GAG_train_y <- GAG_train_data$VISITNO
  
  GAG_test_x <- as.matrix(subset(GAG_test_data, select = -c(PTID, VISITNO)))
  GAG_test_y <- GAG_test_data$VISITNO
  
  GAG_placebo_x <- as.matrix(subset(GAG_placebo_data, select = -c(PTID, VISITNO)))
  GAG_placebo_y <- GAG_placebo_data$VISITNO
  
  # Trains the 'glmnet' classifier using cross-validation.
  ENV_glmnet_cv <- cv.glmnet(x = ENV_train_x, y = ENV_train_y, family = "binomial")
  GAG_glmnet_cv <- cv.glmnet(x = GAG_train_x, y = GAG_train_y, family = "binomial")

  # Computes classification accuracies in two different ways. If the visits are
  # paired by patient, then we compute the proportion of correctly classified
  # patients. Otherwise, we calculate the proportion of visits that are correctly
  # classified.
  if (paired) {
    ENV_predictions_treated <- as.vector(predict(ENV_glmnet_cv, ENV_test_x,
                                                 s = "lambda.min", type = "response"))
    ENV_predictions_placebo <- as.vector(predict(ENV_glmnet_cv, ENV_placebo_x,
                                                 s = "lambda.min", type = "response"))
    GAG_predictions_treated <- as.vector(predict(GAG_glmnet_cv, GAG_test_x,
                                                 s = "lambda.min", type = "response"))
    GAG_predictions_placebo <- as.vector(predict(GAG_glmnet_cv, GAG_placebo_x,
                                                 s = "lambda.min", type = "response"))
   
    # If the difference in classification probabilities exceeds the probability
    # threshold, we assign the first sample as visit 2 and the second as visit 12.
    # Otherwise, we assign the first sample as visit 12 and the second as visit 2.
    ENV_correct_patients <- tapply(seq_along(ENV_test_data$PTID), ENV_test_data$PTID, function(i) {
      if (diff(ENV_predictions_treated[i]) > prob_threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == ENV_test_y[i])
    })
    ENV_correct_placebo <- tapply(seq_along(ENV_placebo_data$PTID), ENV_placebo_data$PTID, function(i) {
      if (diff(ENV_predictions_placebo[i]) > prob_threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == ENV_placebo_y[i])
    })

    GAG_correct_patients <- tapply(seq_along(GAG_test_data$PTID), GAG_test_data$PTID, function(i) {
      if (diff(GAG_predictions_treated[i]) > prob_threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == GAG_test_y[i])
    })
    GAG_correct_placebo <- tapply(seq_along(GAG_placebo_data$PTID), GAG_placebo_data$PTID, function(i) {
     if (diff(GAG_predictions_placebo[i]) > prob_threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == GAG_placebo_y[i])
    })

    ENV_accuracy_treated <- mean(ENV_correct_patients)
    ENV_accuracy_placebo <- mean(ENV_correct_placebo)
    GAG_accuracy_treated <- mean(GAG_correct_patients)
    GAG_accuracy_placebo <- mean(GAG_correct_placebo)
  } else {
    # Computes ENV accuracy
    ENV_predictions_treated <- as.vector(predict(ENV_glmnet_cv, ENV_test_x,
                                                 s = "lambda.min", type = "class"))
    ENV_predictions_placebo <- as.vector(predict(ENV_glmnet_cv, ENV_placebo_x,
                                                 s = "lambda.min", type = "class"))
    ENV_accuracy_treated <- mean(ENV_predictions_treated == ENV_test_y)
    ENV_accuracy_placebo <- mean(ENV_predictions_placebo == ENV_placebo_y)
  
    # Computes GAG accuracy
    GAG_predictions_treated <- as.vector(predict(GAG_glmnet_cv, GAG_test_x,
                                                 s = "lambda.min", type = "class"))
    GAG_predictions_placebo <- as.vector(predict(GAG_glmnet_cv, GAG_placebo_x,
                                                 s = "lambda.min", type = "class"))
    GAG_accuracy_treated <- mean(GAG_predictions_treated == GAG_test_y)
    GAG_accuracy_placebo <- mean(GAG_predictions_placebo == GAG_placebo_y)
  }  
  
  # Determines which features should be kept for classification
  # If present, we manually remove the "(Intercept)"
  # In the case that all markers are removed and only the intercept remains, we
  # add the "(Intercept)" back, for clarity.
  ENV_coef_glmnet <- coef(ENV_glmnet_cv)
  ENV_markers_kept <- rownames(ENV_coef_glmnet)[as.vector(ENV_coef_glmnet) != 0]
  ENV_markers_kept <- ENV_markers_kept[!grepl("(Intercept)", ENV_markers_kept)]
  if (length(ENV_markers_kept) == 0) {
    ENV_markers_kept <- "(Intercept)"
  }
  
  GAG_coef_glmnet <- coef(GAG_glmnet_cv)
  GAG_markers_kept <- rownames(GAG_coef_glmnet)[as.vector(GAG_coef_glmnet) != 0]
  GAG_markers_kept <- GAG_markers_kept[!grepl("(Intercept)", GAG_markers_kept)]
  if (length(GAG_markers_kept) == 0) {
    GAG_markers_kept <- "(Intercept)"
  }

  out <- list(accuracy = list(GAG_treatment = GAG_accuracy_treated,
                GAG_placebo = GAG_accuracy_placebo,
                ENV_treatment = ENV_accuracy_treated,
                ENV_placebo = ENV_accuracy_placebo),
              markers = list(GAG = GAG_markers_kept, ENV = ENV_markers_kept))
  if (paired) {
    # We return the classification probabilities from 'glmnet' and the
    # corresponding test data sets for further analysis (e.g., ROC curves)
    out$classification_probs <- list(GAG_treated = GAG_predictions_treated,
      GAG_placebo = GAG_predictions_placebo,
      ENV_treated = ENV_predictions_treated,
      ENV_placebo = ENV_predictions_placebo)
    out$test_data <- list(GAG_treated = GAG_test_data,
      GAG_placebo = GAG_placebo_data, ENV_treated = ENV_test_data,
      ENV_placebo = ENV_placebo_data)
  }
  out
}


#' Removes quantile information from marker names
#'
#' @param markers character vector of marker names
#' @param quantiles character vector of quantiles to remove from marker names
#' @return character vector of cleaned marker names
strip_quantiles <- function(markers, quantiles = c("9950", "9990", "9999")) {
  for (quant in quantiles) {
    markers <- gsub(quant, "", markers)
  }
  markers
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

#' Computes discrete derivative of a vector based on its kernel density estimate
#'
#' Calculates the kernel density estimate (KDE) of \code{x} using
#' \code{\link{density}} and then calculates the discrete derivative based on
#' the KDE.
#'
#' Based on the following Stack Overflow post:
#' http://bit.ly/Z4g68V
#' 
#' @param obj an object returned from \code{\link{density}}
#' @return list consisting of the sorted values in \code{x} with some centering
#' as well as the discrete slope in \code{y}
deriv_smooth <- function(obj) {
  list(x = rowMeans(embed(obj$x, 2)),
       y = diff(obj$y) / diff(obj$x))
}

#' Computes discrete derivative of a vector based on its kernel density estimate
#'
#' Calculates the kernel density estimate (KDE) of \code{x} using
#' \code{\link{density}} and then calculates the discrete derivative based on
#' the KDE.
#'
#' Based on the following Stack Overflow post:
#' http://bit.ly/Z4g68V
#' 
#' @param obj an object returned from \code{\link{density}}
#' @return list consisting of the sorted values in \code{x} with some centering
#' as well as the discrete slope in \code{y}
second_deriv_smooth <- function(obj, smooth = FALSE, span = 0.2) {
  deriv_out <- deriv_smooth(obj)
  second_deriv <- list(x = rowMeans(embed(deriv_out$x, 2)),
                       y = diff(deriv_out$y) / diff(deriv_out$x))

  # If selected, applies LOESS with minor smoothing to second derivative curve
  if (smooth) {
    loess_out <- loess(y ~ x, data = second_deriv, span = span)
    second_deriv$y <- predict(loess_out, second_deriv$x)
  }

  second_deriv  
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
standardize_cytokines <- function(x, ref_stimulation = "sebctrl") {
  # Scales the data with respect to the negative component of the SEB control
  # sample.
  x_standardize <- x[x$Stim == ref_stimulation, ]$value
  ref_peaks <- openCyto:::find_peaks(x_standardize, order = TRUE, adjust = 3)
  ref_peaks <- ref_peaks[ref_peaks >= 0]
  neg_peak <- ref_peaks[1]
  pos_peak <- ref_peaks[2]

  # If two peaks are present, select cutpoint as mindensity
  # Otherwise, select based on smoothed second derivative
  if (!is.na(pos_peak)) {
    valleys <- openCyto:::find_valleys(x_standardize, adjust = 2)
    cutpoint <- valleys[findInterval(valleys, c(neg_peak, pos_peak)) == 1][1]
  } else {
    second_deriv_standardize <- second_deriv_smooth(x_standardize, n = 4096, adjust = 2)
    second_deriv_standardize <- do.call(cbind.data.frame, second_deriv_standardize)

    # Applies LOESS with minor smoothing to second derivative curve
    loess_out <- loess(y ~ x, data = second_deriv_standardize, span = 0.2)

    # Select cutpoint as first valley after peak
    # Because the LOESS prediction is likely unequal to the negative peak found
    # using openCyto, we find the first valley of the second derivative greater
    # than the valley nearest to the negative_peak
    x_sorted <- sort(x_standardize)
    predict_loess <- predict(loess_out, x_sorted)
    discrete_second_deriv <- diff(sign(diff(predict_loess)))
    which_minima <- which(discrete_second_deriv == 2) + 1
    valleys <- x_sorted[which_minima]
    cutpoint <- valleys[which.min(abs(valleys - neg_peak)) + 1]
  }

  # Calculates standard deviation estimate for the negative-component
  # observations that are greater than the negative peak.
  x_pos <- x_standardize[findInterval(x_standardize, c(neg_peak, cutpoint)) == 1]
  sd_neg <- sqrt(mean((x_pos - neg_peak)^2))

  # Standardizes the cytokine samples within stimulation group with respect to
  # the reference stimulation group.
  # First, centers the values by the mode of the kernel density estimate for
  # the stimulation group.
  x <- tapply(x$value, x$Stim, center_mode)

  # Scales the nonreference stimulation groups by their Huber estimator of the
  # standard deviation and rescales with respect to the reference stimulation
  # sample to put all stimulations on the same scale.
  for(stim in names(x)[!grepl(ref_stimulation, names(x))]) {
    x[[stim]] <- scale_huber(x[[stim]], center = FALSE)
  }

  # For the reference stimulation group (i.e., SEB controls), we scale by the
  # standard deviation of its negative component instead of the Huber estimator.
  x[[ref_stimulation]] <- x[[ref_stimulation]] / sd_neg

  x <- melt(lapply(x, identity))
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
    cbind("Stim" = as.character(pData_gs$Stim[i]), "value" = exprs(flow_set[[1]])[, cytokine2channel(cytokine)])
  })
  x <- do.call(rbind, x)
  x <- data.frame(x, stringsAsFactors = FALSE)
  x$value <- as.numeric(x$value)
  x
}

