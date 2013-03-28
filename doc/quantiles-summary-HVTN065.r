#' # Classification Study - HVTN065

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
               cache = FALSE, echo = FALSE, fig.path = 'figure/HVTN065-',
               cache.path = 'cache/HVTN065-', out.width = "7.5in", out.height = "9.25in")

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
# To speed up the processing, we use a combination of plyr and foreach.
registerDoMC(2)
results_paired <- dlply(cytokine_combinations, .(TNFa, IFNg, IL2), function(cyto_quantiles) {
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

  classification_summary(popstats_combo, treatment_info, paired = TRUE,
                                prob_threshold = 0)
}, .parallel = TRUE)


#+ classification_results_paired, results='asis'

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

Cytokine_Combination <- dlply(m_accuracy, .(TNFa, IFNg, IL2), function(cyto_combo) {
  cyto_combo <- as.numeric(c(cyto_combo$TNFa[1], cyto_combo$IFNg[1], cyto_combo$IL2[1])) / 1e4
  paste(cyto_combo, collapse = "\n")
})
m_accuracy$Cytokine_Combination <- do.call(c, Cytokine_Combination)

p <- ggplot(m_accuracy, aes(x = Cytokine_Combination, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(Stimulation ~ .) + ylim(0, 1)
p <- p + xlab("Cytokine Quantiles (TNFa, IFNg, IL2)") + ylab("Classification Accuracy")
p + ggtitle("Cytokine-Quantile Classification Accuracy of Visit Numbers Paired by Patient")

# Extracts the classification accuracies for each cytokine combination.
accuracy_results <- lapply(results_paired, function(x) x$accuracy)
accuracy_results <- do.call(rbind.data.frame, accuracy_results)

cytokine_combinations_numeric <- apply(cytokine_combinations, 2, function(x) {
  as.numeric(x) / 1e4
})

accuracy_results_numeric <- cbind(cytokine_combinations_numeric, accuracy_results)
rownames(accuracy_results_numeric) <- NULL

print(xtable(accuracy_results_numeric, digits = 4,
             caption = "Classification Results"), type = "html")
 

#' We summarize the results for the top 3 cytokine-quantile combinations for each
#' stimulation group. We choose the top 3 to be the largest differences in the
#' classification accuracies between the treatment and placebo groups.

#+ top_results_summary
seq_top <- seq_len(3)

accuracy_results <- ddply(accuracy_results, .(TNFa, IFNg, IL2), transform,
                          diff_GAG = GAG_treatment - GAG_placebo,
                          diff_ENV = ENV_treatment - ENV_placebo)

which_GAG_top <- order(accuracy_results$diff_GAG, decreasing = TRUE)[seq_top]
which_ENV_top <- order(accuracy_results$diff_ENV, decreasing = TRUE)[seq_top]

GAG_top <- accuracy_results[which_GAG_top, ]
ENV_top <- accuracy_results[which_ENV_top, ]

#' For the top 3 cytokine-quantile combinations from each stimulation group, we
#' provide the markers that were selected by 'glmnet'.

#+ top_markers_summary, results="asis"

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

GAG_xtable <- xtable(GAG_top_markers, digits = 4,
                     caption = "Markers Selected by {\tt glmnet} for Top 3 GAG Cytokine-Quantile Combinations")
print(GAG_xtable, include.rownames = FALSE, type = "html")

ENV_xtable <- xtable(ENV_top_markers, digits = 4,
                     caption = "Markers Selected by {\tt glmnet} for Top 3 ENV Cytokine-Quantile Combinations")
print(ENV_xtable, include.rownames = FALSE, type = "html")

#+ ROC

# We use the 'ROCR' package to construct ROC curves for each of the GAG and ENV
# stimulations.
ENV_probs <- 1 - as.vector(predict(ENV_glmnet_cv, ENV_test_x, s = "lambda.min", type = "response"))
ENV_prediction <- prediction(ENV_probs, ENV_test_y)
ENV_performance <- performance(ENV_prediction, "tpr", "fpr")
ENV_cutoffs <- performance(ENV_prediction, "acc")

GAG_probs <- 1 - as.vector(predict(GAG_glmnet_cv, GAG_test_x, s = "lambda.min", type = "response"))
GAG_prediction <- prediction(GAG_probs, GAG_test_y)
GAG_performance <- performance(GAG_prediction, "tpr", "fpr")
GAG_cutoffs <- performance(GAG_prediction, "acc")

# Next, we combine the ROC curves for each stimulation to construct a singlet
# ggplot2 figure to plot the False Positive Rate vs True Positive Rate.
stim_performance <- rbind.data.frame(
  cbind(ENV_performance@x.values[[1]], ENV_performance@y.values[[1]], "ENV"),
  cbind(GAG_performance@x.values[[1]], GAG_performance@y.values[[1]], "GAG")
)
colnames(stim_performance) <- c("FPR", "TPR", "Stimulation")
stim_performance$FPR <- as.numeric(as.character(stim_performance$FPR))
stim_performance$TPR <- as.numeric(as.character(stim_performance$TPR))

# Next, we combine the ROC curves for each stimulation to construct a singlet
# ggplot2 figure to plot the False Positive Rate vs True Positive Rate.
stim_cutoffs <- rbind.data.frame(
  cbind(ENV_cutoffs@x.values[[1]], ENV_cutoffs@y.values[[1]], "ENV"),
  cbind(GAG_cutoffs@x.values[[1]], GAG_cutoffs@y.values[[1]], "GAG")
)
colnames(stim_cutoffs) <- c("Cutoff", "Accuracy", "Stimulation")
stim_cutoffs$Cutoff <- as.numeric(as.character(stim_cutoffs$Cutoff))
stim_cutoffs$Accuracy <- as.numeric(as.character(stim_cutoffs$Accuracy))


p <- ggplot(stim_performance, aes(x = FPR, y = TPR, color = Stimulation))
p + geom_line(aes(linetype = Stimulation)) + ggtitle("ROC Curve by Stimulation")

# The cutoff value on the x-axis is the probability threshold used to determine
# if a sample is classified as post-vaccine. For instance, if the threshold is
# 0.6, then samples with probability estimates of class membership above 0.6 are
# classified as post-vaccine.
# The classification accuracy is given on the y-axis.
p <- ggplot(stim_cutoffs, aes(x = Cutoff, y = Accuracy, color = Stimulation))
p + geom_line(aes(linetype = Stimulation)) + ggtitle("Accuracy by Probability Threshold")

 

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


 

