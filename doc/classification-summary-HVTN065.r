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
  cache = FALSE, echo = FALSE, fig.path = "figure/HVTN065-", cache.path = "cache/HVTN065-", 
  fig.width = 12, fig.height = 12)
options(stringsAsFactors = FALSE)

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

# Prettifies popstats rownames
popstats_cleaned <- pretty_popstats(popstats_HVTN065)
counts_cleaned <- pretty_popstats(counts_HVTN065)
rownames_popstats <- rownames(popstats_cleaned)
which_cytokines <- grep("TNFa|IFNg|IL2", rownames_popstats)

# Determines which cytokine tolerances were used.
cytokine_tolerances <- sapply(strsplit(rownames_popstats[which_cytokines], "_"), tail, n = 1)
cytokine_tolerances <- sort(unique(cytokine_tolerances))

# Partitions the population statistics into upstream and cytokines for each
# processing below. The cytokine population statistics are stored in a named
# list, where each element corresponds to a cytokine tolerance value.
# The tolerance value is then stripped from the marker names.
#
# We also remove TNFa-IFNg-IL2- and IL2+ | IFNg+
popstats_upstream <- popstats_cleaned[-which_cytokines, ]
popstats_cytokines <- popstats_cleaned[which_cytokines, ]
popstats_cytokines <- popstats_cytokines[!grepl("TNFa-IFNg-IL2-", rownames(popstats_cytokines)), ]
popstats_cytokines <- popstats_cytokines[!grepl("IL2\\+\\|IFNg\\+", rownames(popstats_cytokines)), ]
popstats_cytokines <- partition_popstats(popstats_cytokines,
                                         tolerances = cytokine_tolerances)

# Repeats the previous step but with the counts.
counts_upstream <- counts_cleaned[-which_cytokines, ]
counts_cytokines <- counts_cleaned[which_cytokines, ]
counts_cytokines <- partition_popstats(counts_cytokines,
                                       tolerances = cytokine_tolerances)



#' ## Classification Accuracies of Visit Times Paired by Patients
#'
#' ## Summary of Simulation Design
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
GAG_results <- lapply(popstats_cytokines, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "GAG-1-PTEG", alpha = 0.5)
ENV_results <- lapply(popstats_cytokines, classification_summary,
                        treatment_info = treatment_info, pdata = pData_HVTN065,
                        stimulation = "ENV-1-PTEG", alpha = 0.5)

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy))
colnames(accuracy_results) <- c("Stimulation", "Treatment", "Accuracy", "Tolerance")

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

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC))

#+ ROC_figure, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.5
AUC_df$y <- replace(rep(0.5, nrow(AUC_df)), AUC_df$Stimulation == "GAG-1-PTEG", 0.55)

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation))
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves for ENV and GAG Stimulated Samples")
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

#' ### Markers Selected by `glmnet` for GAG-1-PTEG
#+ markers_GAG, results='asis'
print(xtable(GAG_markers), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for ENV-1-PTEG
#+ markers_ENV, results='asis'
print(xtable(ENV_markers), include.rownames = FALSE, type = "html")

#' ## MIMOSA Results (Coming Soon)
#' Next, we apply MIMOSA to the gating results.

#+ mimosa_setup, eval=FALSE
set.seed(42)
z <- rnorm(5)







#+ manual_gates, eval=FALSE
HVTN065_manual_gates <- subset(HVTN065_manual_gates, ANTIGEN %in% c("ENV-1-PTEG", 
  "GAG-1-PTEG", "negctrl"))
HVTN065_manual_gates <- subset(HVTN065_manual_gates, PTID %in% levels(pData_HVTN065$PTID))
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

# Constructs a ggplot2-friendly results data frame
GAG_accuracy <- melt(lapply(GAG_results, function(x) rbind(unlist(x$accuracy))))[, -1]
ENV_accuracy <- melt(lapply(ENV_results, function(x) rbind(unlist(x$accuracy))))[, -1]
accuracy_results <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_accuracy),
                          cbind(Stimulation = "ENV-1-PTEG", ENV_accuracy))
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

# Estimates ROC curves for GAG and ENV
ROC_curves <- rbind(cbind(Stimulation = "GAG-1-PTEG", GAG_ROC),
                    cbind(Stimulation = "ENV-1-PTEG", ENV_ROC))

#+ ROC_positivity_figure, results='asis'

# The 'pracma:::trapz' function numerically integrates via the trapezoid method
AUC_df <- ddply(ROC_curves, .(Stimulation, Tolerance), summarize,
                AUC = trapz(FPR, TPR))
AUC_df$AUC <- paste("AUC:", round(AUC_df$AUC, 3))
AUC_df$x <- 0.5
AUC_df$y <- replace(rep(0.5, nrow(AUC_df)), AUC_df$Stimulation == "GAG-1-PTEG", 0.55)

# Creates a single plot containing ROC curves for both stimulation groups.
# Also, displays AUC's as text on plot.
p <- ggplot(ROC_curves, aes(x = FPR, y = TPR))
p <- p + geom_line(aes(color = Stimulation))
p <- p + facet_grid(. ~ Tolerance, labeller = label_both)
p <- p + geom_text(data = AUC_df, aes(label = AUC, x = x, y = y, color = Stimulation))
p <- p + theme_bw() + ggtitle("ROC Curves for ENV and GAG Stimulated Samples")
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")


#' ### Manual Gates

#+ classification_positivity_manual, eval=FALSE
