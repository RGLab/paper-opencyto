#' # Classification Study - HVTN065
#'
#' In this document we conduct a classification study based on the automated
#' gating results from the `OpenCyto` R package applied to the HVTN065 data set.
#' We have constructed the cytokine gates using the smoothed first derivative
#' gating scheme with several candidate values of the cutoff tolerance.
#'
#' For each T-cell subset CD4 and CD8, we construct as features the following 7
#' gates:
#'
#' 1. TNFa+IFNg+IL2+
#' 2. TNFa+IFNg+IL2-
#' 3. TNFa+IFNg-IL2+
#' 4. TNFa+IFNg-IL2-
#' 5. TNFa-IFNg+IL2+
#' 6. TNFa-IFNg+IL2-
#' 7. TNFa-IFNg-IL2+
#'
#' We explicity ignore TNFa-IFNg-IL2- because it is redundant given the other seven
#' features.
#' 
#' Our goal is to examine the classification accuracy of patients' visits (i.e.,
#' pre- and post-vaccine).

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = "default", dev = "png", message = FALSE, warning = FALSE, 
  cache = TRUE, echo = FALSE, fig.path = "figure/HVTN065-", cache.path = "cache/HVTN065-", 
  fig.width = 16, fig.height = 16, autodep = TRUE)
options(stringsAsFactors = FALSE)

#+ load_data
setwd("~/rglab/papers/paper-opencyto")
library(ProjectTemplate)
load.project()

#+ prettify_results
colnames(pData_HVTN065) <- c("Sample", "PTID", "Stimulation", "VISITNO")

# Adds treatment group information to the proportion summary
treatment_info <- data.frame(PTID = gsub("-", "", as.character(treatment.HVTN065$Ptid)), 
  Treatment = "Treatment", stringsAsFactors = FALSE)
treatment_info$Treatment <- replace(treatment_info$Treatment, grep("^Placebo", treatment.HVTN065$rx), 
  "Placebo")

# Prettifies popstats rownames
popstats_cleaned <- pretty_popstats(popstats_HVTN065)
counts_cleaned <- pretty_popstats(counts_HVTN065)

#+ partitions_proportions

# Finally, we remove the proportions of the marginal cytokines because we do not
# wish to include them as features, given that we are using the polyfunctional
# proportions as features.
# For example, we remove "cd4:TNFa+_1e-1"

# Extracts the proportions for the upstream gates (CD3, CD4, and CD8)
which_upstream <- !grepl("TNFa|IFNg|IL2", rownames(popstats_cleaned))
popstats_upstream <- popstats_cleaned[which_upstream, ]
popstats_cleaned <- popstats_cleaned[!which_upstream, ]

# Determines which cytokine tolerances were used.
cytokine_tolerances <- sapply(strsplit(rownames(popstats_cleaned), "_"), tail, n = 1)
cytokine_tolerances <- sort(unique(cytokine_tolerances))

# Extracts the marginal cytokine proportions as well as the proportions for the
# polyfunctional gates.
which_cytokines <- grepl("cd[48]:(TNFa|IFNg|IL2)\\+_1e", rownames(popstats_cleaned))
which_polyfunc <- !which_cytokines
popstats_cytokines <- popstats_cleaned[which_cytokines, ]
popstats_polyfunc <- popstats_cleaned[which_polyfunc, ]

# We update the population statistics for the polyfunctional gates and store
# them in a named list, where each element corresponds to a cytokine tolerance
# value. The tolerance value is then stripped from the marker names.
#
# We also remove TNFa-IFNg-IL2- and IL2+ | IFNg+. The TNFa-IFNg-IL2- is
# redundant because the proportions sum to 1.  Also, we consider IL2+ | IFNg+
# separately in a comparison with manual gates.
popstats_polyfunc <- popstats_polyfunc[!grepl("TNFa-IFNg-IL2-", rownames(popstats_polyfunc)), ]
popstats_polyfunc <- popstats_polyfunc[!grepl("IL2\\+\\|IFNg\\+", rownames(popstats_polyfunc)), ]
popstats_polyfunc <- partition_popstats(popstats_polyfunc,
                                         tolerances = cytokine_tolerances)

#+ partitions_counts, eval=FALSE

# TODO: Update the summary of counts for the MIMOSA summary below.
# Repeats the previous step but with the counts.
rownames_counts <- rownames(counts_cleaned)
which_tcell_subsets <- grep("cd[48]$", rownames_counts)
which_marginal_cytokines <- grep("cd[48]:(TNFa|IFNg|IL2)\\+_1e", rownames_counts)

counts_tcell_subsets <- counts_cleaned[which_tcell_subsets, ]
counts_cytokines <- counts_cleaned[which_marginal_cytokines, ]
counts_cytokines <- partition_popstats(counts_cytokines,
                                       tolerances = cytokine_tolerances)
counts_cytokines <- lapply(counts_cytokines, function(x) {
  rbind(counts_tcell_subsets, x)
})

#' ## Classification Accuracies of Visit Times Paired by Patients
#'
#' ### Summary of Simulation Design
#'
#' First, we partition the HVTN065 patients by their treatment status into a
#' treatment group and placebo group. Of the patients in the treatment group, we
#' randomly partition 60% of the them into a training data set and the remaining
#' 40% of the patients into a test data set.
#' We utilize the `glmnet` package using the elastic net with `alpha = 0.5` to
#' build a classifier from the population proportions for the markers and the
#' polyfunctional gates obtained using the `OpenCyto` package. Next, because
#' there are two visits (i.e., pre- and post-vaccine) for each patient, we pair
#' the visits in the test data set by patient. For each patient-visit pairing, we
#' classify the two samples and calculate the difference in their classification
#' probabilities. Let d = Pr(sample 1 from subject 1 = post-vaccine) -
#' Pr(sample 2 from subject 1 = post-vaccine).
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
GAG_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG", alpha = 0.5)
ENV_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG", alpha = 0.5)
POL_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "POL-1-PTEG", alpha = 0.5)

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
POL_accuracy <- melt(lapply(POL_results, function(x) rbind(unlist(x$accuracy))))[, -1]

accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy),
                          cbind(Stimulation = "POL-1-PTEG", POL_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

#' ### Results

# + classification_paired_figure, results='asis'
p <- ggplot(accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
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
GAG_ROC <- ROC_summary(GAG_results, cytokine_tolerances)
ENV_ROC <- ROC_summary(ENV_results, cytokine_tolerances)
POL_ROC <- ROC_summary(POL_results, cytokine_tolerances)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC),
                    cbind(Stimulation = "POL-1-PTEG", POL_ROC))

#+ ROC_figure, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.75
AUC_df$y <- rep(c(0.55, 0.5, 0.45), each = length(cytokine_tolerances))

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")

#' Here, we provide the markers that were selected by 'glmnet' for each
#' stimulation group. In the case that `(Intercept)` is given, no markers are
#' selected by `glmnet`, leaving only an intercept term.
#+ markers_selected
GAG_markers <- lapply(GAG_results, function(x) x$markers)
num_markers <- max(sapply(GAG_markers, length))
GAG_markers <- do.call(cbind, lapply(GAG_markers, function(x) {
  length(x) <- num_markers
  x
}))
GAG_markers <- data.frame(GAG_markers, check.names = FALSE)

ENV_markers <- lapply(ENV_results, function(x) x$markers)
num_markers <- max(sapply(ENV_markers, length))
ENV_markers <- do.call(cbind, lapply(ENV_markers, function(x) {
  length(x) <- num_markers
  x
}))
ENV_markers <- data.frame(ENV_markers, check.names = FALSE)

POL_markers <- lapply(POL_results, function(x) x$markers)
num_markers <- max(sapply(POL_markers, length))
POL_markers <- do.call(cbind, lapply(POL_markers, function(x) {
  length(x) <- num_markers
  x
}))
POL_markers <- data.frame(POL_markers, check.names = FALSE)

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for POL-1-PTEG
#+ markers_POL, results='asis'
print(xtable(POL_markers), include.rownames = FALSE, type = "html")


### LASSO

#' Next, we repeat the classification study using LASSO rather than elastic net.
#' That is, we set the value of 'alpha = 1' in the 'glmnet' call.

#+ classification_paired_LASSO

set.seed(42)
GAG_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG", alpha = 1)
ENV_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG", alpha = 1)
POL_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "POL-1-PTEG", alpha = 1)

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
POL_accuracy <- melt(lapply(POL_results, function(x) rbind(unlist(x$accuracy))))[, -1]

accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy),
                          cbind(Stimulation = "POL-1-PTEG", POL_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

#' ### LASSO Results

# + classification_paired_figure_LASSO, results='asis'
p <- ggplot(accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
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

#+ ROC_LASSO
# Summarizes the classification probabilties for GAG and ENV.
GAG_ROC <- ROC_summary(GAG_results, cytokine_tolerances)
ENV_ROC <- ROC_summary(ENV_results, cytokine_tolerances)
POL_ROC <- ROC_summary(POL_results, cytokine_tolerances)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC),
                    cbind(Stimulation = "POL-1-PTEG", POL_ROC))

#+ ROC_figure_LASSO, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.75
AUC_df$y <- rep(c(0.55, 0.5, 0.45), each = length(cytokine_tolerances))

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")

#' Here, we provide the markers that were selected by 'glmnet' for each
#' stimulation group. In the case that `(Intercept)` is given, no markers are
#' selected by `glmnet`, leaving only an intercept term.
#+ markers_selected_LASSO
GAG_markers <- lapply(GAG_results, function(x) x$markers)
num_markers <- max(sapply(GAG_markers, length))
GAG_markers <- do.call(cbind, lapply(GAG_markers, function(x) {
  length(x) <- num_markers
  x
}))
GAG_markers <- data.frame(GAG_markers, check.names = FALSE)

ENV_markers <- lapply(ENV_results, function(x) x$markers)
num_markers <- max(sapply(ENV_markers, length))
ENV_markers <- do.call(cbind, lapply(ENV_markers, function(x) {
  length(x) <- num_markers
  x
}))
ENV_markers <- data.frame(ENV_markers, check.names = FALSE)

POL_markers <- lapply(POL_results, function(x) x$markers)
num_markers <- max(sapply(POL_markers, length))
POL_markers <- do.call(cbind, lapply(POL_markers, function(x) {
  length(x) <- num_markers
  x
}))
POL_markers <- data.frame(POL_markers, check.names = FALSE)

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG_LASSO, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV_LASSO, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for POL-1-PTEG
#+ markers_POL_LASSO, results='asis'
print(xtable(POL_markers), include.rownames = FALSE, type = "html")


#' ### Cytokine MFIs as Features

#' Next, we investigate the classification performance by also using median
#' flourescence intensities (MFIs) as features.

#+ classification_paired_MFI
set.seed(42)
GAG_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG", alpha = 0.5, other_features = mfi_cytokines)
ENV_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG", alpha = 0.5, other_features = mfi_cytokines)
POL_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "POL-1-PTEG", alpha = 0.5, other_features = mfi_cytokines)

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
POL_accuracy <- melt(lapply(POL_results, function(x) rbind(unlist(x$accuracy))))[, -1]

accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy),
                          cbind(Stimulation = "POL-1-PTEG", POL_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

#+ classification_paired_figure_MFI, results='asis'
p <- ggplot(accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
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

#+ ROC_MFI
# Summarizes the classification probabilties for GAG and ENV.
GAG_ROC <- ROC_summary(GAG_results, cytokine_tolerances)
ENV_ROC <- ROC_summary(ENV_results, cytokine_tolerances)
POL_ROC <- ROC_summary(POL_results, cytokine_tolerances)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC),
                    cbind(Stimulation = "POL-1-PTEG", POL_ROC))

#+ ROC_figure_MFI, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.75
AUC_df$y <- rep(c(0.55, 0.5, 0.45), each = length(cytokine_tolerances))

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")

#' Here, we provide the markers that were selected by 'glmnet' for each
#' stimulation group. In the case that `(Intercept)` is given, no markers are
#' selected by `glmnet`, leaving only an intercept term.
#+ markers_selected_MFI
GAG_markers <- lapply(GAG_results, function(x) x$markers)
num_markers <- max(sapply(GAG_markers, length))
GAG_markers <- do.call(cbind, lapply(GAG_markers, function(x) {
  length(x) <- num_markers
  x
}))
GAG_markers <- data.frame(GAG_markers, check.names = FALSE)

ENV_markers <- lapply(ENV_results, function(x) x$markers)
num_markers <- max(sapply(ENV_markers, length))
ENV_markers <- do.call(cbind, lapply(ENV_markers, function(x) {
  length(x) <- num_markers
  x
}))
ENV_markers <- data.frame(ENV_markers, check.names = FALSE)

POL_markers <- lapply(POL_results, function(x) x$markers)
num_markers <- max(sapply(POL_markers, length))
POL_markers <- do.call(cbind, lapply(POL_markers, function(x) {
  length(x) <- num_markers
  x
}))
POL_markers <- data.frame(POL_markers, check.names = FALSE)

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG_MFI, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV_MFI, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for POL-1-PTEG
#+ markers_POL_MFI, results='asis'
print(xtable(POL_markers), include.rownames = FALSE, type = "html")


#' ### Cytokine MFIs as Features using LASSO

#' Next, we investigate the classification performance by also using median
#' flourescence intensities (MFIs) as features. Again, we use a LASSO classifier
#' via 'glmnet' by setting 'alpha = 1'.

#+ classification_paired_MFI_LASSO
set.seed(42)
GAG_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG", alpha = 1, other_features = mfi_cytokines)
ENV_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG", alpha = 1, other_features = mfi_cytokines)
POL_results <- lapply(popstats_polyfunc, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "POL-1-PTEG", alpha = 1, other_features = mfi_cytokines)

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
POL_accuracy <- melt(lapply(POL_results, function(x) rbind(unlist(x$accuracy))))[, -1]

accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy),
                          cbind(Stimulation = "POL-1-PTEG", POL_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

#+ classification_paired_figure_MFI_LASSO, results='asis'
p <- ggplot(accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
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

#+ ROC_MFI_LASSO
# Summarizes the classification probabilties for GAG and ENV.
GAG_ROC <- ROC_summary(GAG_results, cytokine_tolerances)
ENV_ROC <- ROC_summary(ENV_results, cytokine_tolerances)
POL_ROC <- ROC_summary(POL_results, cytokine_tolerances)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC),
                    cbind(Stimulation = "POL-1-PTEG", POL_ROC))

#+ ROC_figure_MFI_LASSO, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.75
AUC_df$y <- rep(c(0.55, 0.5, 0.45), each = length(cytokine_tolerances))

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")

#' Here, we provide the markers that were selected by 'glmnet' for each
#' stimulation group. In the case that `(Intercept)` is given, no markers are
#' selected by `glmnet`, leaving only an intercept term.
#+ markers_selected_MFI_LASSO
GAG_markers <- lapply(GAG_results, function(x) x$markers)
num_markers <- max(sapply(GAG_markers, length))
GAG_markers <- do.call(cbind, lapply(GAG_markers, function(x) {
  length(x) <- num_markers
  x
}))
GAG_markers <- data.frame(GAG_markers, check.names = FALSE)

ENV_markers <- lapply(ENV_results, function(x) x$markers)
num_markers <- max(sapply(ENV_markers, length))
ENV_markers <- do.call(cbind, lapply(ENV_markers, function(x) {
  length(x) <- num_markers
  x
}))
ENV_markers <- data.frame(ENV_markers, check.names = FALSE)

POL_markers <- lapply(POL_results, function(x) x$markers)
num_markers <- max(sapply(POL_markers, length))
POL_markers <- do.call(cbind, lapply(POL_markers, function(x) {
  length(x) <- num_markers
  x
}))
POL_markers <- data.frame(POL_markers, check.names = FALSE)

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG_MFI_LASSO, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV_MFI_LASSO, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for POL-1-PTEG
#+ markers_POL_MFI_LASSO, results='asis'
print(xtable(POL_markers), include.rownames = FALSE, type = "html")


#' ## Classification with IFNg+ | IL2+

#' Here, we compare the automated gates constructed using OpenCyto with the
#' manually constructed gates. For both cases, we apply logistic regression to
#' classify patient visits utilizing a single feature, namely IFNg+ | IL2+.

#' ### Automated Gates using OpenCyto

#+ classification_positivity_opencyto
popstats_positivity <- popstats_cleaned[grepl("IL2\\+\\|IFNg\\+", rownames(popstats_cleaned)), ]
popstats_positivity <- partition_popstats(popstats_positivity,
                                         tolerances = cytokine_tolerances)
set.seed(42)
GAG_results <- lapply(popstats_positivity, classification_summary_logistic,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG")
ENV_results <- lapply(popstats_positivity, classification_summary_logistic,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG")
POL_results <- lapply(popstats_positivity, classification_summary_logistic,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "POL-1-PTEG")

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
POL_accuracy <- melt(lapply(POL_results, function(x) rbind(unlist(x$accuracy))))[, -1]
accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy),
                          cbind(Stimulation = "POL-1-PTEG", POL_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

# + classification_positivity_figure, results='asis'
p <- ggplot(accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
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

#+ ROC_positivity
# Summarizes the classification probabilties for GAG and ENV.
GAG_ROC <- ROC_summary(GAG_results, cytokine_tolerances)
ENV_ROC <- ROC_summary(ENV_results, cytokine_tolerances)
POL_ROC <- ROC_summary(POL_results, cytokine_tolerances)

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC),
                    cbind(Stimulation = "POL-1-PTEG", POL_ROC))

#+ ROC_positivity_figure, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.5
AUC_df$y <- rep(c(0.55, 0.5, 0.45), each = length(cytokine_tolerances))

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")


#' ### Manual Gates

#+ classification_positivity_manual
manual_GAG_results <- classification_summary_logistic(HVTN065_manual, treatment_info = treatment_info,
                                                      pdata = pData_HVTN065, stimulation = "GAG-1-PTEG",
                                                      manual_gates = TRUE)
manual_ENV_results <- classification_summary_logistic(HVTN065_manual, treatment_info = treatment_info,
                                                      pdata = pData_HVTN065, stimulation = "ENV-1-PTEG",
                                                      manual_gates = TRUE)

# Constructs a ggplot2-friendly results data frame
GAG_manual_accuracy <- cbind.data.frame(Stimulation = "GAG-1-PTEG", rbind(unlist(manual_GAG_results$accuracy)))
ENV_manual_accuracy <- cbind.data.frame(Stimulation = "ENV-1-PTEG", rbind(unlist(manual_ENV_results$accuracy)))

manual_accuracy_results <- melt(rbind(GAG_manual_accuracy, ENV_manual_accuracy))
colnames(manual_accuracy_results) <- c("Stimulation", "Treatment", "Accuracy")

manual_accuracy_results$Treatment <- factor(as.character(manual_accuracy_results$Treatment), labels = c("placebo", "treatment"))

# + classification_positivity_figure, results='asis'
p <- ggplot(manual_accuracy_results, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + ylim(0, 1) + xlab("Stimulation Group") + ylab("Classification Accuracy")
p <- p + ggtitle("Classification Accuracy of Visit Numbers Paired by Patient")
p <- p + theme_bw() + theme(plot.title = element_text(size = 18))
p <- p + theme(strip.text.y = element_text(size = 14))
p <- p + theme(axis.text = element_text(size = 12))
p + theme(axis.title = element_text(size = 16))


#+ ROC_manual
# Summarizes the classification probabilties for GAG and ENV.
GAG_ROC_manual <- ROC_summary_manual(manual_GAG_results)
ENV_ROC_manual <- ROC_summary_manual(manual_ENV_results)

# Estimates ROC curves for GAG and ENV
ROC_curves_manual <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC_manual),
                           cbind(Stimulation = "ENV-1-PTEG", ENV_ROC_manual))

#+ ROC_manual_figure, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df_manual <- ddply(ROC_curves_manual, .(Stimulation), summarize, AUC = trapz(FPR, TPR))
AUC_df_manual$AUC <- paste("AUC:", round(AUC_df_manual$AUC, 3))
AUC_df_manual$x <- 0.5
AUC_df_manual$y <- replace(rep(0.5, nrow(AUC_df_manual)), AUC_df_manual$Stimulation == "GAG-1-PTEG", 0.55)

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves_manual, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation), size = 1.5)
p <- p + geom_text(data = AUC_df_manual, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves by Stimulation Group")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")


#' ## MIMOSA Results
#' Next, we apply MIMOSA to the gating results.

#+ mimosa_setup, eval=FALSE

# Formats cytokine counts to use with MIMOSA.
m_counts <- melt(counts_cytokines)
colnames(m_counts) <- c("Marker", "Sample", "Count", "Tolerance")
m_counts$Tcell <- "CD4"
m_counts$Tcell <- with(m_counts, replace(Tcell, grepl("cd8", Marker), "CD8"))
m_counts$Marker <- gsub("cd[48]:", "", m_counts$Marker)
m_counts$Marker <- gsub("\\+", "", m_counts$Marker)

counts_cytokines <- ddply(m_counts, .(Tolerance, Sample, Tcell), function(x) {
  which_tcell <- grepl("cd[48]", x$Marker)
  tcells <- x[which_tcell, ]
  cytokines <- x[-which_tcell, ]

  cbind.data.frame(Cytokine = cytokines$Marker,
                   negative = tcells$Count - cytokines$Count,
                   positive = cytokines$Count)
})
counts_cytokines <- plyr:::join(counts_cytokines, pData_HVTN065, by = "Sample")

# We stored the 'negctrl' with sample numbers appended to the strings so that
# plotGate could identify unique samples. Here, we strip the sample numbers to
# summaricounts_cytokinese the negative controls as a whole.
counts_cytokines$Stimulation <- as.character(counts_cytokines$Stimulation)
counts_cytokines$Stimulation <- replace(counts_cytokines$Stimulation,
                                        grep("^negctrl", counts_cytokines$Stimulation), "negctrl")
counts_cytokines <- subset(counts_cytokines, select = -Sample)

# Divvies up the cytokines into GAG and ENV
GAG_counts <- subset(counts_cytokines, Stimulation %in% c("GAG-1-PTEG", "negctrl"))
ENV_counts <- subset(counts_cytokines, Stimulation %in% c("ENV-1-PTEG", "negctrl"))

# Greg provided a fix to handle the case of two negative controls
GAG_counts <- ddply(GAG_counts, .(Tolerance, Tcell, Cytokine, PTID, VISITNO, Stimulation), summarize,
                    negative = sum(negative), positive = sum(positive))
ENV_counts <- ddply(ENV_counts, .(Tolerance, Tcell, Cytokine, PTID, VISITNO, Stimulation), summarize,
                    negative = sum(negative), positive = sum(positive))


# Creates ExpressionSets for MIMOSA
GAG_Eset <- ConstructMIMOSAExpressionSet(GAG_counts, reference = Stimulation %in% "negctrl",
                                      measure.columns = c("negative", "positive"),
                                      other.annotations = c("Tolerance", "Tcell", "Cytokine", "PTID", "Stimulation", "VISITNO"),
                                      default.cast.formula = component ~ Tolerance + Tcell + Cytokine + PTID + Stimulation + VISITNO,
                                      .variables = .(Tolerance, Tcell, Cytokine, PTID, VISITNO))

ENV_Eset <- ConstructMIMOSAExpressionSet(ENV_counts, reference = Stimulation %in% "negctrl",
                                      measure.columns = c("negative", "positive"),
                                      other.annotations = c("Tolerance", "Tcell", "Cytokine", "PTID", "Stimulation", "VISITNO"),
                                      default.cast.formula = component ~ Tolerance + Tcell + Cytokine + PTID + Stimulation + VISITNO,
                                      .variables = .(Tolerance, Tcell, Cytokine, PTID, VISITNO))


# Applies MIMOSA for both stimulation groups
GAG_MIMOSA <- MIMOSA(negative + positive ~ PTID + VISITNO | Tcell + Cytokine + Stimulation + Tolerance,
                     data = GAG_Eset, subset = RefTreat %in% "Treatment",
                     ref = RefTreat %in% "Reference", method = "mcmc", iter = 250000,
                     burn = 50000, pXi = 1, EXPRATE = 1e-5, run.parallel = TRUE)

ENV_MIMOSA <- MIMOSA(negative + positive ~ PTID + VISITNO | Tcell + Cytokine + Stimulation + Tolerance,
                     data = ENV_Eset, subset = RefTreat %in% "Treatment",
                     ref = RefTreat %in% "Reference", method = "mcmc", iter = 350000,
                     burn = 50000, pXi = 1, EXPRATE = 1e-5, run.parallel = TRUE)


#+ MIMOSA_GAG_volcanoplot, eval=FALSE
q <- as.vector(sapply(seq_along(GAG_MIMOSA), function(i) {
  -log10(fdr(GAG_MIMOSA[[i]]@z))
}))
ps <- as.vector(sapply(seq_along(GAG_MIMOSA), function(i) {
  prop.table(as.matrix(GAG_MIMOSA[[i]]@result@n.stim),1)[,2]
}))
pu <- as.vector(sapply(seq_along(GAG_MIMOSA), function(i) {
  prop.table(as.matrix(GAG_MIMOSA[[i]]@result@n.unstim),1)[,2]
}))
p.stim <- as.vector(sapply(seq_along(GAG_MIMOSA), function(i) {
  GAG_MIMOSA[[i]]@z[,2]
}))
mimosa.call<-as.vector(sapply(seq_along(GAG_MIMOSA),function(i){
  MIMOSA:::fdr(GAG_MIMOSA[[i]]@z)
}))
pd<-do.call(rbind,lapply(GAG_MIMOSA,pData))
pd<-data.frame(pd,ps,pu,p.stim,mimosa.call=mimosa.call<0.01)
ggplot(pd)+geom_point(aes(x=ps-pu,y=p.stim,color=mimosa.call))+facet_grid(Tcell~Tolerance+Cytokine,scale="free_x")+theme(axis.text.x=element_text(angle=90))

num_tolerances <- length(cytokine_tolerances)
df <- data.frame(p.stim, q, ps, pu, Tolerance = gl(length(GAG_MIMOSA), length(q) / length(GAG_MIMOSA), labels = seq_along(GAG_MIMOSA)))
p <- ggplot(df) + geom_point(aes(x = ps - pu, y = p.stim)) + theme_bw()
p <- p + scale_y_continuous("Probability of Stimulation") + facet_wrap(~ Tolerance, scale = "free_x")
p + ggtitle("Volcano Plot - GAG-1-PTEG")

#+ MIMOSA_ENV_volcanoplot, eval=FALSE
q <- as.vector(sapply(seq_along(ENV_MIMOSA), function(i) {
  -log10(fdr(ENV_MIMOSA[[i]]@z))
}))
ps <- as.vector(sapply(seq_along(ENV_MIMOSA), function(i) {
  prop.table(as.matrix(ENV_MIMOSA[[i]]@result@data$n.stim),1)[,2]
}))
pu <- as.vector(sapply(seq_along(ENV_MIMOSA), function(i) {
  prop.table(as.matrix(ENV_MIMOSA[[i]]@result@data$n.unstim),1)[,2]
}))
p.stim <- as.vector(sapply(seq_along(ENV_MIMOSA), function(i) {
  ENV_MIMOSA[[i]]@z[,2]
}))
mimosa.call<-as.vector(sapply(seq_along(ENV_MIMOSA),function(i){
  MIMOSA:::fdr(ENV_MIMOSA[[i]]@z)
}))
pd<-do.call(rbind,lapply(ENV_MIMOSA,pData))
pd<-data.frame(pd,ps,pu,p.stim,mimosa.call=mimosa.call<0.01)
ggplot(pd)+geom_point(aes(x=ps-pu,y=p.stim,color=mimosa.call))+facet_grid(Tcell~Tolerance+Cytokine,scale="free_x")+theme(axis.text.x=element_text(angle=90))
ggplot(pd)+geom_boxplot(aes(x=mimosa.call,y=ps-pu,color=Cytokine))+facet_wrap(~Tcell+Tolerance,scale="free_y")

num_tolerances <- length(cytokine_tolerances)
df <- data.frame(p.stim, q, ps, pu, Tolerance = gl(num_tolerances, length(q) / num_tolerances, labels = cytokine_tolerances))
p <- ggplot(df) + geom_point(aes(x = ps - pu, y = p.stim)) + theme_bw()
p <- p + scale_y_continuous("Probability of Stimulation") + facet_wrap(~ Tolerance)
p + ggtitle("Volcano Plot - ENV-1-PTEG")
