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
#' The \code{permutations} function is from the \code{gregmisc} package on CRAN.
#' @param markers character vector of marker names
#' @return vector containing all combinations of the markers
#' @examples
#' polyfunction_nodes(c("IFNg", "IL2", "TNFa", "GzB", "CD57"))
polyfunction_nodes <- function(markers) {
  require('gregmisc')
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
#' @param train_pct a numeric value determining the percentage of treated
#' patients used as training data and the remaining patients as test data
#' @param paired logical value that determines if the classification accuracies
#' should be calculated by pairing visits for a patient
#' @param prob_threshold a numeric value above which the difference in
#' classification probabilities for two samples from the same subject indicates
#' that the first sample is classified as post-vaccine. Ignored if \code{paired}
#' is \code{FALSE}. See details.
#' @return TODO
classification_summary <- function(popstats, treatment_info, train_pct = 0.6, paired = FALSE,
                                   prob_threshold = 0) {
  m_popstats <- melt(popstats)
  colnames(m_popstats) <- c("Marker", "Sample", "Proportion")
  m_popstats$Marker <- as.character(m_popstats$Marker)
  m_popstats <- plyr:::join(m_popstats, pData_HVTN065, by = "Sample")
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

