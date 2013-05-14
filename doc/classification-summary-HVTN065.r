#' # Classification Study - HVTN065
#'
#' In this document we conduct a classification study based on the automated
#' gating results from the `OpenCyto` R package applied to the HVTN065 data set.
#' We have constructed the cytokine gates using the smoothed first derivative
#' gating scheme with `tol = 1e-3`. Then, we construct the following $8 = 2^3$
#' polyfunctional gates:
#'
#' 1. TNFa+IFNg+IL2+
#' 2. TNFa+IFNg+IL2-
#' 3. ...
#' 8. TNFa-IFNg-IL2-
#'
#' Our goal is to examine the classification accuracy of patients' visits (i.e.,
#' pre- and post-vaccine).

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = "default", dev = "png", message = FALSE, warning = FALSE, 
  cache = TRUE, echo = FALSE, fig.path = "figure/HVTN065-", cache.path = "cache/HVTN065-", 
  fig.width = 8, fig.height = 8, autodep = TRUE)

#+ load_data
setwd("..")
library(ProjectTemplate)
load.project()

#+ prettify_results
colnames(pData_HVTN065) <- c("Sample", "PTID", "Stimulation", "VISITNO")

# Adds treatment group information to the proportion summary
treatment_info <- data.frame(PTID = gsub("-", "", as.character(treatment.HVTN065$Ptid)), 
  Treatment = "Treatment", stringsAsFactors = FALSE)
treatment_info$Treatment <- replace(treatment_info$Treatment, grep("^Placebo", treatment.HVTN065$rx), 
  "Placebo")

# Remove all population statistics for the following markers
markers_remove <- c("root", "cd8gate_pos", "cd4_neg", "cd8gate_neg", "cd4_pos")
popstats_remove <- sapply(strsplit(rownames(popstats_HVTN065), "/"), tail, n = 1)
popstats_remove <- popstats_remove %in% markers_remove
popstats_HVTN065 <- popstats_HVTN065[!popstats_remove, ]

rownames_popstats <- rownames(popstats_HVTN065)

# Updates all cytokine combinations having the name of the form 'cd4/TNFa' to
# 'TNFa'
which_combo <- grep("&", rownames_popstats)
rownames_combo <- rownames_popstats[which_combo]
rownames_combo <- gsub("cd[48]:TNFa", "TNFa", rownames_combo)
rownames_combo <- gsub("cd[48]:IFNg", "IFNg", rownames_combo)
rownames_combo <- gsub("cd[48]:IL2", "IL2", rownames_combo)
rownames_popstats[which_combo] <- rownames_combo

# Updates all cytokines gates to the form 'cd4:TNFa'
which_cytokines <- grep("TNFa|IFNg|IL2", rownames_popstats)
rownames_cytokines <- rownames_popstats[which_cytokines]
rownames_cytokines <- sapply(strsplit(rownames_cytokines, "cd3/"), tail, n = 1)
rownames_cytokines <- gsub("/", ":", rownames_cytokines)
rownames_popstats[which_cytokines] <- rownames_cytokines

# Retains the last marker name for all non-cytokine gates
rownames_noncytokines <- rownames_popstats[-which_cytokines]
rownames_noncytokines <- sapply(strsplit(rownames_noncytokines, "/"), tail, n = 1)
rownames_popstats[-which_cytokines] <- rownames_noncytokines

# Reformats cytokine-marker combinations: !TNFa&IFNg&IL2 => TNFa-IFNg+IL2+
rownames_popstats <- gsub("TNFa", "TNFa+", rownames_popstats)
rownames_popstats <- gsub("!TNFa\\+", "TNFa-", rownames_popstats)

rownames_popstats <- gsub("IFNg", "IFNg+", rownames_popstats)
rownames_popstats <- gsub("!IFNg\\+", "IFNg-", rownames_popstats)

rownames_popstats <- gsub("IL2", "IL2+", rownames_popstats)
rownames_popstats <- gsub("!IL2\\+", "IL2-", rownames_popstats)

rownames_popstats <- gsub("&", "", rownames_popstats)

# Updates popstats rownames
rownames(popstats_HVTN065) <- rownames_popstats

#' ## Classification Accuracies of Visit Times Paired by Patients
#'
#' ## Summary of Simulation Design
#'
#' First, we partition the HVTN065 patients by their treatment status into a
#' treatment group and placebo group. Of the patients in the treatment group, we
#' randomly partition 60% of the them into a training data set and the remaining
#' 40% of the patients into a test data set.
#' We utilize the `glmnet` package to build a classifier from the population
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

#+ classification_paired

# Per Greg: 'We also want to do this paired, using the difference in
# classification probabilities for two samples from the same subject. i.e.  d =
# Pr(sample 1 from subject 1 = post-vaccine) - Pr(sample 2 from subject 1 =
# post-vaccine).  If d > threshold, then classify sample 1 as post-vaccine and
# sample 2 as pre-vaccine, otherwise if d < threshold classify sample 1 as
# pre-vaccine and sample 2 as post-vaccine, otherwise mark them as
# unclassifiable.'

# Computes the paired classification results.  The classification accuracies
# computed use a probability threshold of 0 (i.e., d = 0).
# The value of 'alpha = 0.5' is passed to 'glmnet' to indicate the usage of
# elastic net.
set.seed(42)
GAG_results <- classification_summary(popstats_HVTN065, treatment_info, pData_HVTN065, 
  stimulation = "GAG-1-PTEG", alpha = 0.5)
ENV_results <- classification_summary(popstats_HVTN065, treatment_info, pData_HVTN065, 
  stimulation = "ENV-1-PTEG", alpha = 0.5)

# Constructs a ggplot2-friendly results data frame
accuracy_results <- rbind.data.frame(unlist(GAG_results$accuracy), unlist(ENV_results$accuracy))
colnames(accuracy_results) <- c("treatment", "placebo")
accuracy_results$Stimulation <- c("GAG-1-PTEG", "ENV-1-PTEG")
m_accuracy <- melt(accuracy_results)
colnames(m_accuracy) <- c("Stimulation", "Treatment", "Accuracy")

# + classification_paired_figure, results='asis'
p <- ggplot(m_accuracy, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + ylim(0, 1) + xlab("Stimulation Group") + ylab("Classification Accuracy")
p <- p + ggtitle("Classification Accuracy of Visit Numbers Paired by Patient")
p <- p + theme_bw() + theme(plot.title = element_text(size = 18))
p <- p + theme(strip.text.y = element_text(size = 14))
p <- p + theme(axis.text = element_text(size = 12))
p + theme(axis.title = element_text(size = 16))


#' Next, we calculate ROC curves assuming all vaccinees are true positive and
#' the placebos are false positive.
#' For each PTID, we compute the absolute value of the difference in
#' classification probabilties for visits 2 and 12 and then order by the
#' differences.

#+ ROC
# Summarizes the classification probabilties for GAG and ENV.
GAG_treated <- cbind(
                     subset(GAG_results$test_data$treated, select = c(PTID, VISITNO)),
                     Truth = "Treatment",
                     Probability = GAG_results$classification_probs$treated)
GAG_placebo <- cbind(
                     subset(GAG_results$test_data$placebo, select = c(PTID, VISITNO)),
                     Truth = "Placebo",
                     Probability = GAG_results$classification_probs$placebo)
GAG_probs <- rbind(GAG_treated, GAG_placebo)

ENV_treated <- cbind(
                     subset(ENV_results$test_data$treated, select = c(PTID, VISITNO)),
                     Truth = "Treatment",
                     Probability = ENV_results$classification_probs$treated)
ENV_placebo <- cbind(
                     subset(ENV_results$test_data$placebo, select = c(PTID, VISITNO)),
                     Truth = "Placebo",
                     Probability = ENV_results$classification_probs$placebo)
ENV_probs <- rbind(ENV_treated, ENV_placebo)
                   
# For each PTID, we compute the absolute value of the difference in
# classification probabilties for visits 2 and 12 and then order by the
# differences.
GAG_summary <- ddply(GAG_probs, .(PTID), summarize,
                     delta = 1 - abs(diff(Probability[VISITNO %in% c("2", "12")])),
                     Truth = unique(Truth))
GAG_summary <- GAG_summary[with(GAG_summary, order(delta, Truth, decreasing = FALSE)), ]

ENV_summary <- ddply(ENV_probs, .(PTID), summarize,
                     delta = 1 - abs(diff(Probability[VISITNO %in% c("2", "12")])),
                     Truth = unique(Truth))
ENV_summary <- ENV_summary[with(ENV_summary, order(delta, Truth, decreasing = FALSE)), ]

# Calculates true and false positive rates based on Treatment and Placebo
# samples, respectively. Because we are using Treatments and Placebos, we
# calculate TPRs and FPRs differently than usual. The basic idea is that when we
# add to the TPR each time we classify a patient as Treatment and to the FPR
# each time we classify a patient as Placebo. The ordering here is determined by
# the rank of the differences in classification probabilities.
GAG_FPR <- with(GAG_summary, cumsum(Truth == "Placebo") / sum(Truth == "Placebo"))
GAG_TPR <- with(GAG_summary, cumsum(Truth == "Treatment") / sum(Truth == "Treatment"))

ENV_FPR <- with(ENV_summary, cumsum(Truth == "Placebo") / sum(Truth == "Placebo"))
ENV_TPR <- with(ENV_summary, cumsum(Truth == "Treatment") / sum(Truth == "Treatment"))

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
GAG_AUC <- trapz(GAG_FPR, GAG_TPR)
ENV_AUC <- trapz(ENV_FPR, ENV_TPR)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(
                    cbind.data.frame(Stimulation = "GAG-1-PTEG", FPR = GAG_FPR, TPR = GAG_TPR),
                    cbind.data.frame(Stimulation = "ENV-1-PTEG", FPR = ENV_FPR, TPR = ENV_TPR)
              )

#+ ROC_figure, results='asis'
AUC_df <- cbind(Stimulation = c("GAG-1-PTEG", "ENV-1-PTEG"),
                AUC = paste("AUC:", round(c(GAG_AUC, ENV_AUC), 3)))
AUC_df <- data.frame(AUC_df, stringsAsFactors = FALSE)
AUC_df$x <- c(0.5, 0.5)
AUC_df$y <- c(0.55, 0.5)

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation))
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves for ENV and GAG Stimulated Samples")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")

#' Here, we provide the markers that were selected by 'glmnet' for each
#' stimulation group. In the case that `(Intercept)` is given, no markers are
#' selected by `glmnet`, leaving only an intercept term.
#+ markers_selected
GAG_markers <- data.frame(Markers = GAG_results$markers, stringsAsFactors = FALSE)
ENV_markers <- data.frame(Markers = ENV_results$markers, stringsAsFactors = FALSE)

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")



#+ manual_gates, eval=FALSE
HVTN065_manual_gates <- subset(HVTN065_manual_gates, ANTIGEN %in% c("ENV-1-PTEG", 
  "GAG-1-PTEG", "negctrl"))
HVTN065_manual_gates <- subset(HVTN065_manual_gates, PTID %in% levels(m_pop_stats$PTID))
HVTN065_manual_gates <- subset(HVTN065_manual_gates, VISITNO %in% c(2, 12))

manual_counts <- ddply(HVTN065_manual_gates, .(PTID, VISITNO, ANTIGEN), function(x) {
  counts <- data.frame(PTID = x$PTID[1], VISITNO = x$VISITNO[1], ANTIGEN = x$ANTIGEN[1], 
    root = no_commas(x$COLLECTCT)[1], Singlet = no_commas(x$SUBSET1_NUM)[1], 
    Live = no_commas(x$SUBSET2_NUM)[1], Lymphocytes = no_commas(x$SUBSET3_NUM)[1], 
    CD3 = no_commas(x$SUBSET4_NUM)[1], stringsAsFactors = FALSE)
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
  proportions <- with(x, data.frame(PTID = PTID, VISITNO = VISITNO, ANTIGEN = ANTIGEN, 
    Singlet = Singlet/root, Live = Live/Singlet, Lymphocytes = Lymphocytes/Live, 
    CD3 = CD3/Lymphocytes, CD4 = CD4/CD3, CD8 = CD8/CD3))
  
  for (CD4_cytokine in colnames(x)[grep("CD4:", colnames(x))]) {
    proportions[[CD4_cytokine]] <- x[[CD4_cytokine]]/x[["CD4"]]
  }
  
  for (CD8_cytokine in colnames(x)[grep("CD8:", colnames(x))]) {
    proportions[[CD8_cytokine]] <- x[[CD8_cytokine]]/x[["CD8"]]
  }
  
  proportions
})


# TODO: Wait for full cytokine-gate information On 11 March 2013, Greg said he
# would look into it.  After we have the full information, build a classifier
# for the manual proportions.
 
