#' # Classification Study - HVTN065
#'

#' In this document we conduct a classification study based on the automated
#' gating results from the `OpenCyto` R package applied to the HVTN065 data set.
#' In particular, we look at various combinations of quantiles for the TNFa, IFNg,
#' and IL2 cytokines. For each of these three cytokines, we applied the quantiles
#' 0.995, 0.999, and 0.9999. For one of the quantile combinations, we applied
#' `flowClust` to fit the data for the given marker and applied the quantile
#' specified as the gate. Then, we constructed polyfunctional gates for the three
#' cytokines. For instance, suppose that we employed the following 3 quantiles:
#'
#' 1. TNFa: 0.999
#' 2. IFNg: 0.995
#' 3. IL2: 0.999
#' 
#' After constructing the gates for these quantiles, we construct the following
#' $8 = 2^3$ polyfunctional gates.
#'
#' 1. TNFa+IFNg+IL2+
#' 2. TNFa+IFNg+IL2-
#' 3. ...
#' 8. TNFa-IFNg-IL2-
#'
#' Our goal is to determine the quantiles that yield the best classification
#' accuracy of patients' visits (i.e., pre- and post-vaccine).

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE,
               cache = TRUE, echo = FALSE, fig.path = 'figure/HVTN065-',
               cache.path = 'cache/HVTN065-', fig.width = 18, fig.height = 18,
               autodep = TRUE)

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()

#+ prettify_results
colnames(pData_HVTN065) <- c("Sample", "PTID", "Stimulation", "VISITNO")

# Adds treatment group information to the proportion summary
treatment_info <- data.frame(PTID = gsub("-", "", as.character(treatment.HVTN065$Ptid)),
                             Treatment = "Treatment", stringsAsFactors = FALSE)
treatment_info$Treatment <- replace(treatment_info$Treatment, grep("^Placebo", treatment.HVTN065$rx), "Placebo")

# Remove all population statistics for the following markers
markers_remove <- c("root", "cd8gate_pos", "cd4_neg", "cd8gate_neg", "cd4_pos")
popstats_remove <- sapply(strsplit(rownames(popstats_HVTN065), "/"), tail, n = 1)
popstats_remove <- popstats_remove %in% markers_remove
popstats_HVTN065 <- popstats_HVTN065[!popstats_remove, ]

rownames_popstats <- rownames(popstats_HVTN065)

# Updates all cytokine combinations having the name of the form "cd4/TNFa" to "TNFa"
which_combo <- grep("&", rownames_popstats)
rownames_combo <- rownames_popstats[which_combo]
rownames_combo <- gsub("cd4/TNFa", "TNFa", rownames_combo)
rownames_combo <- gsub("cd8/TNFa", "TNFa", rownames_combo)
rownames_combo <- gsub("cd4/IFNg", "IFNg", rownames_combo)
rownames_combo <- gsub("cd8/IFNg", "IFNg", rownames_combo)
rownames_combo <- gsub("cd4/IL2", "IL2", rownames_combo)
rownames_combo <- gsub("cd8/IL2", "IL2", rownames_combo)
rownames_popstats[which_combo] <- rownames_combo

# Updates all cytokines gates to the form "cd4:TNFa"
which_cytokines <- grep("TNFa|IFNg|IL2", rownames_popstats)
rownames_cytokines <- rownames_popstats[which_cytokines]
rownames_cytokines <- sapply(strsplit(rownames_cytokines, "cd3/"), tail, n = 1)
rownames_cytokines <- gsub("/", ":", rownames_cytokines)
rownames_popstats[which_cytokines] <- rownames_cytokines

# Retains the last marker name for all non-cytokine gates
rownames_noncytokines <- rownames_popstats[-which_cytokines]
rownames_noncytokines <- sapply(strsplit(rownames_noncytokines, "/"), tail, n = 1)
rownames_popstats[-which_cytokines] <- rownames_noncytokines

# Reformats cytokine-marker combinations
# The cytokine markers currently have the form:
# !TNFa9950&IFNg9999&IL29990
# We update them to have the form:
# TNFa9950-IFNg9999+IL29990+
TNFa <- IFNg <- IL2 <- c("9950", "9990", "9999")

for (quant in TNFa) {
  quant <- paste0("TNFa", quant)
  rownames_popstats <- gsub(paste0(quant), paste0(quant, "+"), rownames_popstats)
  rownames_popstats <- gsub(paste0("!", quant, "\\+"), paste0(quant, "-"), rownames_popstats)
}

for (quant in IFNg) {
  quant <- paste0("IFNg", quant)
  rownames_popstats <- gsub(paste0(quant), paste0(quant, "+"), rownames_popstats)
  rownames_popstats <- gsub(paste0("!", quant, "\\+"), paste0(quant, "-"), rownames_popstats)
}

for (quant in IL2) {
  quant <- paste0("IL2", quant)
  rownames_popstats <- gsub(paste0(quant), paste0(quant, "+"), rownames_popstats)
  rownames_popstats <- gsub(paste0("!", quant, "\\+"), paste0(quant, "-"), rownames_popstats)
}
rownames_popstats <- gsub("&", "", rownames_popstats)

# Updates popstats rownames
rownames(popstats_HVTN065) <- rownames_popstats

#' ## Classification Accuracies of Visit Times Paired by Patients
#'
#' ## Summary of Simulation Design
#' For each cytokine-quantile combination, we estimate classification accuracies
#' using the following scheme:
#'
#' First, we partition the HVTN065 patients by their treatment status into a
#' treatment group and placebo group. Of the patients in the treatment group, we
#' randomly partition 60% of the them into a training data set and the remaining
#' 40% of the patients into a test data set. For the given cytokine quantiles,
#' we utilize the `glmnet` package to build a classifier from the population
#' proportions for the markers and the polyfunctional gates obtained using the
#' `OpenCyto` package. Next, because there are two visits (i.e., pre- and
#' post-vaccine) for each patient, we pair the visits in the test data set by
#' patient. For each patient-visit pairing, we classify the two samples
#' and calculate the difference in their classification probabilities. Let
#' d = Pr(sample 1 from subject 1 = post-vaccine) - Pr(sample 2 from subject 1 = post-vaccine).
#' For a given probability threshold, if $d > threshold$, then we classify sample
#' 1 as post-vaccine and sample 2 as pre-vaccine. Otherwise, if $d < threshold$,
#' we classify sample 1 as pre-vaccine and sample 2 as post-vaccine. We calculate
#' the classification accuracy as the number of correctly classified patients. In
#' the same manner we calculate the classification accuracy using the placebo
#' patients as a separate test data set.
#'
#' We present both a graphical and tabular summary of the results.
#'

#+ cytokine_quantile_classification_paired
# Per Greg: "We also want to do this paired, using the difference in
# classification probabilities for two samples from the same subject. i.e.
# d = Pr(sample 1 from subject 1 = post-vaccine) - Pr(sample 2 from subject 1 = post-vaccine).
# If d > threshold, then classify sample 1 as post-vaccine and sample 2 as
# pre-vaccine, otherwise if d < threshold classify sample 1 as pre-vaccine and
# sample 2 as post-vaccine, otherwise mark them as unclassifiable."

# Loop through each cytokine combination and construct matrix of features
set.seed(42)

# Extracts the population statistics for the markers upstream (i.e., the gates
# before the cytokines)
which_cytokines <- grep("cd[48]:", rownames_popstats)
popstats_upstream <- popstats_HVTN065[-which_cytokines, ]

TNFa <- IFNg <- IL2 <- c("9950", "9990", "9999")

cytokine_combinations <- expand.grid(TNFa = TNFa, IFNg = IFNg, IL2 = IL2,
                                     stringsAsFactors = FALSE)
# To speed up the processing, we use mclapply with 'num_cores' cores.
num_cores <- 12

results_paired <- mclapply(seq_len(nrow(cytokine_combinations)), function(i) {
  cyto_quantiles <- cytokine_combinations[i, ]

  TNFa <- paste0("TNFa", cyto_quantiles$TNFa)
  IFNg <- paste0("IFNg", cyto_quantiles$IFNg)
  IL2 <- paste0("IL2", cyto_quantiles$IL2)

  # Constructs a lookup table for the current cytokine quantiles
  cytokines <- c(TNFa, IFNg, IL2)
  cytokines <- c(paste0("cd4:", cytokines), paste0("cd8:", cytokines))

  # Constructs a lookup table of the cytokine doubles for the current quantile
  # combination
  cytokine_doubles <- polyfunction_nodes(c(IFNg, IL2))
  cytokine_doubles <- c(paste0("cd4:", cytokine_doubles),
                        paste0("cd8:", cytokine_doubles))

  # Constructs a lookup table of the cytokine triples for the current quantile
  # combination
  cytokine_triples <- polyfunction_nodes(c(TNFa, IFNg, IL2))
  cytokine_triples <- c(paste0("cd4:", cytokine_triples),
                        paste0("cd8:", cytokine_triples))

  # Extracts the population statistics for the current cytokine quantiles and
  # combinations
  popstats_cytokines <- popstats_HVTN065[which(rownames_popstats %in% cytokines), ]
  popstats_doubles <- popstats_HVTN065[which(rownames_popstats %in% cytokine_doubles), ]
  popstats_triples <- popstats_HVTN065[which(rownames_popstats %in% cytokine_triples), ]

  # The population statistics upstream and the current cytokine quantiles
  popstats_combo <- rbind(popstats_upstream, popstats_cytokines, popstats_doubles,
                          popstats_triples)

  classification_summary(popstats_combo, treatment_info, pdata = pData_HVTN065,
                         paired = TRUE, prob_threshold = 0)
}, mc.cores = num_cores)

cytokine_combos <- lapply(seq_len(nrow(cytokine_combinations)), function(i) {
  cyto_combo <- cytokine_combinations[i, ]
  paste(cyto_combo, collapse = ".")
})
names(results_paired) <- do.call(c, cytokine_combos)

#+ classification_results_paired

# Extracts the classification accuracies for each cytokine combination.
accuracy_results <- lapply(results_paired, function(x) x$accuracy)
accuracy_results <- do.call(rbind.data.frame, accuracy_results)

accuracy_results <- cbind(cytokine_combinations, accuracy_results)
rownames(accuracy_results) <- NULL

m_accuracy <- melt(accuracy_results, variable.name = "Treatment_Stimulation",
                   value.name = "Accuracy")
m_accuracy$Treatment_Stimulation <- as.character(m_accuracy$Treatment_Stimulation)

m_accuracy$Stimulation <- "GAG"
m_accuracy$Stimulation <- with(m_accuracy,
                               replace(Stimulation, grep("^ENV", Treatment_Stimulation), "ENV"))

m_accuracy$Treatment <- "Treatment"
m_accuracy$Treatment <- with(m_accuracy,
                             replace(Treatment, grep("placebo", Treatment_Stimulation), "Placebo"))

cytokine_labels <- lapply(seq_len(nrow(m_accuracy)), function(i) {
  cyto_combo <- m_accuracy[i, ]
  cyto_combo <- as.numeric(c(cyto_combo$TNFa[1], cyto_combo$IFNg[1], cyto_combo$IL2[1])) / 1e4
  paste(cyto_combo, collapse = "\n")
})
m_accuracy$Cytokine_Combination <- do.call(c, cytokine_labels)

#+ classification_results_figure, results='asis'

p <- ggplot(m_accuracy, aes(x = Cytokine_Combination, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(Stimulation ~ .) + ylim(0, 1)
p <- p + xlab("Cytokine Quantiles (TNFa, IFNg, IL2)") + ylab("Classification Accuracy")
p <- p + ggtitle("Cytokine-Quantile Classification Accuracy of Visit Numbers Paired by Patient") + theme_bw()
p <- p + theme(plot.title  = element_text(size = 18))
p <- p + theme(strip.text.y = element_text(size = 14))
p + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 16))

#+ classification_results_table, results='asis'

# Extracts the classification accuracies for each cytokine combination.
accuracy_results <- lapply(results_paired, function(x) x$accuracy)
accuracy_results <- do.call(rbind.data.frame, accuracy_results)

cytokine_combinations_numeric <- apply(cytokine_combinations, 2, function(x) {
  as.numeric(x) / 1e4
})

accuracy_results_numeric <- cbind(cytokine_combinations_numeric, accuracy_results)
rownames(accuracy_results_numeric) <- NULL

print(xtable(accuracy_results_numeric, digits = 4), include.rownames = FALSE, type = "html")
 
#' ## Markers for Top 3 Quantile Combinations for Each Stimulation Group
#'
#' We summarize the results for the top 3 cytokine-quantile combinations for each
#' stimulation group. We choose the top 3 to be the largest differences in the
#' classification accuracies between the treatment and placebo groups.

#+ top_markers_summary
seq_top <- seq_len(3)

accuracy_results <- cbind(cytokine_combinations, accuracy_results)

# To ensure that the ordering of the quantile combinations is preserved,
# notice that we reverse the order of cytokines specified in 'ddply'.
accuracy_results <- ddply(accuracy_results, .(IL2, IFNg, TNFa), transform,
                          diff_GAG = GAG_treatment - GAG_placebo,
                          diff_ENV = ENV_treatment - ENV_placebo)

which_GAG_top <- order(accuracy_results$diff_GAG, decreasing = TRUE)[seq_top]
which_ENV_top <- order(accuracy_results$diff_ENV, decreasing = TRUE)[seq_top]

GAG_top <- accuracy_results[which_GAG_top, ]
ENV_top <- accuracy_results[which_ENV_top, ]

#' For the top 3 cytokine-quantile combinations from each stimulation group, we
#' provide the markers that were selected by 'glmnet'. In the case that
#' `(Intercept)` is given, no markers are selected by `glmnet`, leaving only an
#' intercept term.

GAG_top_markers <- lapply(seq_top, function(i) {
  combo <- GAG_top[i, ]
  markers <- results_paired[[with(combo, paste(TNFa, IFNg, IL2, sep = "."))]]$markers$GAG
  strip_quantiles(markers, quantiles = c("9950", "9990", "9999"))
})
GAG_top_markers <- do.call(rbind, lapply(GAG_top_markers, paste, collapse = ", "))
GAG_top_markers <- cbind.data.frame(cytokine_combinations[which_GAG_top, ], Markers = GAG_top_markers)

ENV_top_markers <- lapply(seq_top, function(i) {
  combo <- ENV_top[i, ]
  markers <- results_paired[[with(combo, paste(TNFa, IFNg, IL2, sep = "."))]]$markers$ENV
  strip_quantiles(markers, quantiles = c("9950", "9990", "9999"))
})
ENV_top_markers <- do.call(rbind, lapply(ENV_top_markers, paste, collapse = ", "))
ENV_top_markers <- cbind.data.frame(cytokine_combinations[which_ENV_top, ], Markers = ENV_top_markers)


#' ### Markers Selected by `glmnet` for Top 3 GAG Cytokine-Quantile Combinations
#+ markers_top3_GAG, results='asis'
print(xtable(GAG_top_markers, digits = 4), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for Top 3 ENV Cytokine-Quantile Combinations
#+ markers_top3_ENV, results='asis'
print(xtable(ENV_top_markers, digits = 4), include.rownames = FALSE, type = "html")

#+ accuracy_by_thresholds

prob_thresholds <- seq(0, 0.5, by = 0.01)

#' ### GAG Classification Accuracy by Probability Thresholds

#+ GAG_threshold_plot

GAG_accuracy_threshold <- lapply(seq_top, function(i) {
  combo <- GAG_top[i, ]
  results_markers <- results_paired[[with(combo, paste(TNFa, IFNg, IL2, sep = "."))]]

  # Extracts PTID and VISITNO from the test data used in the classification study
  # for both the treatment and placebo groups
  test_treated <- subset(results_markers$test_data$GAG_treated, select = c(PTID, VISITNO))
  test_placebo <- subset(results_markers$test_data$GAG_placebo, select = c(PTID, VISITNO))

  # Classification probabilities for treatment and placebo groups
  probs_treated <- results_markers$classification_probs$GAG_treated
  probs_placebo <- results_markers$classification_probs$GAG_placebo

  accuracy_threshold <- lapply(prob_thresholds, function(threshold) {
    correct_treated <- tapply(seq_along(test_treated$PTID), test_treated$PTID, function(i) {
      if (diff(probs_treated[i]) > threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == test_treated$VISITNO[i])
    })
    treated_accuracy <- mean(correct_treated)

    correct_placebo <- tapply(seq_along(test_placebo$PTID), test_placebo$PTID, function(i) {
      if (diff(probs_placebo[i]) > threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == test_placebo$VISITNO[i])
    })
    placebo_accuracy <- mean(correct_placebo)

    list(Treated = treated_accuracy, Placebo = placebo_accuracy)
  })
  accuracy_threshold <- do.call(rbind, accuracy_threshold)
  list(markers = subset(combo, select = c(TNFa, IFNg, IL2)),
       accuracy_threshold = cbind.data.frame(Threshold = as.factor(prob_thresholds),
         accuracy_threshold))
})

GAG_markers <- sapply(GAG_accuracy_threshold, function(x) {
  paste0(c("TNFa", "IFNg", "IL2"), x$markers, collapse = ":")
})

GAG_accuracy_threshold <- lapply(GAG_accuracy_threshold, function(x) x$accuracy_threshold)
names(GAG_accuracy_threshold) <- GAG_markers
m_GAG_thresh <- melt(GAG_accuracy_threshold)
colnames(m_GAG_thresh) <- c("Threshold", "Treatment", "Accuracy", "Quantile_Combo")
m_GAG_thresh$Threshold <- as.numeric(as.character(m_GAG_thresh$Threshold))

p <- ggplot(m_GAG_thresh, aes(x = Threshold, y = Accuracy, color = Treatment, linetype = Treatment))
p <- p + geom_line(size = 2) + facet_wrap(~ Quantile_Combo)
p <- p + ggtitle("Classification Accuracy by Probability Threshold - GAG Stimulation") + theme_bw()
p <- p + theme(legend.title = element_text(size = 14)) + theme(legend.text = element_text(size = 12))
p <- p + theme(plot.title  = element_text(size = 18)) + theme(strip.text.x = element_text(size = 14))
p <- p + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 16))
p + xlab("Probability Threshold")


#' ### ENV Classification Accuracy by Probability Thresholds

#+ ENV_threshold_plot

ENV_accuracy_threshold <- lapply(seq_top, function(i) {
  combo <- ENV_top[i, ]
  results_markers <- results_paired[[with(combo, paste(TNFa, IFNg, IL2, sep = "."))]]

  # Extracts PTID and VISITNO from the test data used in the classification study
  # for both the treatment and placebo groups
  test_treated <- subset(results_markers$test_data$ENV_treated, select = c(PTID, VISITNO))
  test_placebo <- subset(results_markers$test_data$ENV_placebo, select = c(PTID, VISITNO))

  # Classification probabilities for treatment and placebo groups
  probs_treated <- results_markers$classification_probs$ENV_treated
  probs_placebo <- results_markers$classification_probs$ENV_placebo

  accuracy_threshold <- lapply(prob_thresholds, function(threshold) {
    correct_treated <- tapply(seq_along(test_treated$PTID), test_treated$PTID, function(i) {
      if (diff(probs_treated[i]) > threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == test_treated$VISITNO[i])
    })
    treated_accuracy <- mean(correct_treated)

    correct_placebo <- tapply(seq_along(test_placebo$PTID), test_placebo$PTID, function(i) {
      if (diff(probs_placebo[i]) > threshold) {
        classification <- c("2", "12")
      } else {
        classification <- c("12", "2")
      }
      all(classification == test_placebo$VISITNO[i])
    })
    placebo_accuracy <- mean(correct_placebo)

    list(Treated = treated_accuracy, Placebo = placebo_accuracy)
  })
  accuracy_threshold <- do.call(rbind, accuracy_threshold)
  list(markers = subset(combo, select = c(TNFa, IFNg, IL2)),
       accuracy_threshold = cbind.data.frame(Threshold = as.factor(prob_thresholds),
         accuracy_threshold))
})

ENV_markers <- sapply(ENV_accuracy_threshold, function(x) {
  paste0(c("TNFa", "IFNg", "IL2"), x$markers, collapse = ":")
})

ENV_accuracy_threshold <- lapply(ENV_accuracy_threshold, function(x) x$accuracy_threshold)
names(ENV_accuracy_threshold) <- ENV_markers
m_ENV_thresh <- melt(ENV_accuracy_threshold)
colnames(m_ENV_thresh) <- c("Threshold", "Treatment", "Accuracy", "Quantile_Combo")
m_ENV_thresh$Threshold <- as.numeric(as.character(m_ENV_thresh$Threshold))

p <- ggplot(m_ENV_thresh, aes(x = Threshold, y = Accuracy, color = Treatment, linetype = Treatment))
p <- p + geom_line(size = 2) + facet_wrap(~ Quantile_Combo)
p <- p + ggtitle("Classification Accuracy by Probability Threshold - ENV Stimulation") + theme_bw()
p <- p + theme(legend.title = element_text(size = 14)) + theme(legend.text = element_text(size = 12))
p <- p + theme(plot.title  = element_text(size = 18)) + theme(strip.text.x = element_text(size = 14))
p <- p + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 16))
p + xlab("Probability Threshold")


#+ manual_gates, eval=FALSE
HVTN065_manual_gates <- subset(HVTN065_manual_gates, ANTIGEN %in% c("ENV-1-PTEG", "GAG-1-PTEG", "negctrl"))
HVTN065_manual_gates <- subset(HVTN065_manual_gates, PTID %in% levels(m_pop_stats$PTID))
HVTN065_manual_gates <- subset(HVTN065_manual_gates, VISITNO %in% c(2, 12))

manual_counts <- ddply(HVTN065_manual_gates, .(PTID, VISITNO, ANTIGEN), function(x) {
  counts <- data.frame(
    PTID = x$PTID[1],
    VISITNO = x$VISITNO[1],
    ANTIGEN = x$ANTIGEN[1],
    root = no_commas(x$COLLECTCT)[1],
    Singlet = no_commas(x$SUBSET1_NUM)[1],
    Live = no_commas(x$SUBSET2_NUM)[1],
    Lymphocytes = no_commas(x$SUBSET3_NUM)[1],
    CD3 = no_commas(x$SUBSET4_NUM)[1],
    stringsAsFactors = FALSE
  )
  # To grab the appropriate counts from the CD4 and CD8 subtrees, we split the
  # data.
  x_CD4 <- subset(x, SUBSET5 == "CD4+")
  x_CD8 <- subset(x, SUBSET5 == "CD8+")
  
  counts$CD4 <- no_commas(x_CD4$SUBSET5_NUM)[1]
  counts$CD8 <- no_commas(x_CD8$SUBSET5_NUM)[1]
  
  # CD4 Cytokine Counts
  for (i in seq_len(nrow(x_CD4))) {
    cytokine_name <- paste0("CD4:", x_CD4$SUBSET6[i])
    cytokine_count <- no_commas(x_CD4$SUBSET6_NUM[i])
    counts[[cytokine_name]] <- cytokine_count
  }
  
  # CD8 Cytokine Counts
  for (i in seq_len(nrow(x_CD8))) {
    cytokine_name <- paste0("CD8:", x_CD8$SUBSET6[i])
    cytokine_count <- no_commas(x_CD8$SUBSET6_NUM[i])
    counts[[cytokine_name]] <- cytokine_count
  }

  counts
})

# Calculates the proportions for the manual gates
manual_proportions <- ddply(manual_counts, .(PTID, VISITNO, ANTIGEN), function(x) {
  proportions <- with(x, data.frame(
    PTID = PTID,
    VISITNO = VISITNO,
    ANTIGEN = ANTIGEN,
    Singlet = Singlet / root,
    Live = Live / Singlet,
    Lymphocytes = Lymphocytes / Live,
    CD3 = CD3 / Lymphocytes,
    CD4 = CD4 / CD3,
    CD8 = CD8 / CD3))
  
  for (CD4_cytokine in colnames(x)[grep("CD4:", colnames(x))]) {
    proportions[[CD4_cytokine]] <- x[[CD4_cytokine]] / x[["CD4"]]
  }
  
  for (CD8_cytokine in colnames(x)[grep("CD8:", colnames(x))]) {
    proportions[[CD8_cytokine]] <- x[[CD8_cytokine]] / x[["CD8"]]
  }

  proportions
})


# TODO: Wait for full cytokine-gate information
#   On 11 March 2013, Greg said he would look into it.
#   After we have the full information, build a classifier for the manual proportions.

