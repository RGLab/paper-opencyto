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
rownames_combo <- gsub("cd[48]:TNFa", "TNFa", rownames_combo)
rownames_combo <- gsub("cd[48]:IFNg", "IFNg", rownames_combo)
rownames_combo <- gsub("cd[48]:IL2", "IL2", rownames_combo)
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

# Reformats cytokine-marker combinations:
# !TNFa&IFNg&IL2 => TNFa-IFNg+IL2+
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

# Computes the paired classification results.
# The classification accuracies computed use a probability threshold of 0 (i.e., d = 0).

# TODO: After gating memory issue is resolved, comptue summaries for each tolerance value.
#   We may need to grab some of the previous code used for quantiles.
set.seed(42)
GAG_results <- classification_summary(popstats_HVTN065, treatment_info,
                                      pData_HVTN065, stimulation = "GAG-1-PTEG")
ENV_results <- classification_summary(popstats_HVTN065, treatment_info,
                                      pData_HVTN065, stimulation = "ENV-1-PTEG")

#+ classification_results_paired

# Constructs a ggplot2-friendly results data frame
accuracy_results <- rbind.data.frame(unlist(GAG_results$accuracy), unlist(ENV_results$accuracy))
colnames(accuracy_results) <- c("treatment", "placebo")
accuracy_results$Stimulation <- c("GAG-1-PTEG", "ENV-1-PTEG")
m_accuracy <- melt(accuracy_results)
colnames(m_accuracy) <- c("Stimulation", "Treatment", "Accuracy")

#+ classification_results_figure, results='asis'

p <- ggplot(m_accuracy, aes(x = Stimulation, fill = Treatment))
p <- p + geom_bar(aes(weight = Accuracy), position = "dodge")
p <- p + ylim(0, 1)
p <- p + xlab("Stimulation Group") + ylab("Classification Accuracy")
p <- p + ggtitle("Classification Accuracy of Visit Numbers Paired by Patient") + theme_bw()
p <- p + theme(plot.title  = element_text(size = 18))
p <- p + theme(strip.text.y = element_text(size = 14))
p + theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 16))






#' ## Markers for Top 3 Quantile Combinations for Each Stimulation Group
#'
#' We summarize the results for the top 3 cytokine-quantile combinations for
#' each stimulation group. First, we calculate ROC curves assuming all vaccinees
#' are true positive and the placebos are false positive. We choose the top 3
#' models that have the largest area under the curve (AUC) scores based on these
#' ROC curves.

#+ top_markers_summary
seq_top <- seq_len(3)

#+ accuracy_by_thresholds

prob_thresholds <- seq(0, 0.5, by = 0.01)

#' ### GAG Classification Accuracy by Probability Thresholds

#+ GAG_threshold_plot
# To ensure that the ordering of the quantile combinations is preserved,
# notice that we reverse the order of cytokines specified in 'ddply'.
GAG_accuracy_threshold <- dlply(cytokine_combinations, .(IL2, IFNg, TNFa), function(combo) {
  marker_combo <- paste0(c("TNFa", "IFNg", "IL2"), combo, collapse = ":")
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
  cbind.data.frame(Marker_Combo = marker_combo,
                   Threshold = as.factor(prob_thresholds), accuracy_threshold)
})

GAG_accuracy_threshold <- do.call(rbind, GAG_accuracy_threshold)
rownames(GAG_accuracy_threshold) <- NULL

GAG_accuracy_threshold$Treated <- unlist(GAG_accuracy_threshold$Treated)
GAG_accuracy_threshold$Placebo <- unlist(GAG_accuracy_threshold$Placebo)


# Calculates AUC for each marker combination
GAG_AUC_combo <- with(GAG_accuracy_threshold, tapply(seq_along(Marker_Combo), Marker_Combo, function(i) {
  pairwise <- expand.grid(Treated = Treated[i], Placebo = Placebo[i])
  mean(pairwise$Treated > pairwise$Placebo)
}))
GAG_AUC_combo <- melt(GAG_AUC_combo)
colnames(GAG_AUC_combo) <- c("Quantile_Combo", "AUC")



#' ### ENV Classification Accuracy by Probability Thresholds

#+ ENV_threshold_plot
# To ensure that the ordering of the quantile combinations is preserved,
# notice that we reverse the order of cytokines specified in 'ddply'.
ENV_accuracy_threshold <- dlply(cytokine_combinations, .(IL2, IFNg, TNFa), function(combo) {
  marker_combo <- paste0(c("TNFa", "IFNg", "IL2"), combo, collapse = ":")
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
  cbind.data.frame(Marker_Combo = marker_combo,
                   Threshold = as.factor(prob_thresholds), accuracy_threshold)
})

ENV_accuracy_threshold <- do.call(rbind, ENV_accuracy_threshold)
rownames(ENV_accuracy_threshold) <- NULL

ENV_accuracy_threshold$Treated <- unlist(ENV_accuracy_threshold$Treated)
ENV_accuracy_threshold$Placebo <- unlist(ENV_accuracy_threshold$Placebo)


# Calculates AUC for each marker combination
ENV_AUC_combo <- with(ENV_accuracy_threshold, tapply(seq_along(Marker_Combo), Marker_Combo, function(i) {
  pairwise <- expand.grid(Treated = Treated[i], Placebo = Placebo[i])
  mean(pairwise$Treated > pairwise$Placebo)
}))
ENV_AUC_combo <- melt(ENV_AUC_combo)
colnames(ENV_AUC_combo) <- c("Quantile_Combo", "AUC")


# Combines AUCs by stimulation to see comparison
AUC_combo <- rbind(cbind(GAG_AUC_combo, Stimulation = "GAG"), cbind(ENV_AUC_combo, Stimulation = "ENV"))
AUC_combo$Combo_Label <- gsub(":", "\n", AUC_combo$Quantile_Combo)

p <- ggplot(AUC_combo, aes(x = Combo_Label, fill = Stimulation)) + geom_bar(aes(weight = AUC), position = "dodge")
p + xlab("Cytokine Quantile Combination") + ylab("AUC") + ggtitle("Area Under the Curve (AUC) by Cytokine Combination")+theme(axis.text.x=element_text(angle=90))

###Begin Greg's code
#'Use results_paired to construct ROC curves.
require(pracma)
GAG_AUC<-sapply(seq_along(results_paired),function(i){
df<-cbind(truth<-rbind(results_paired[[i]]$test_data$GAG_treated[,c("PTID","VISITNO")],results_paired[[i]]$test_data$GAG_placebo[,c("PTID","VISITNO")]),
P=c(results_paired[[i]]$classification_probs$GAG_treated,results_paired[[i]]$classification_probs$GAG_placebo),truth=c(rep(c("Treatment"),nrow(results_paired[[i]]$test_data$GAG_treated)),rep(c("Placebo"),nrow(results_paired[[i]]$test_data$GAG_placebo))))
df<-ddply(df,.(PTID),summarize,delta=1-diff(P[VISITNO%in%c("2","12")]),truth=unique(truth))
df<-df[order(df$delta,df$truth,decreasing=FALSE),]
o<-order(df$delta,decreasing=FALSE)
FPR<-cumsum(df$truth[o]=="Placebo")/sum(df$truth=="Placebo")
TPR<-cumsum(df$truth[o]=="Treatment")/sum(df$truth=="Treatment")
AUC<-trapz(FPR,TPR)
})
GAG_AUC<-data.frame(thresholds=names(results_paired),AUC=GAG_AUC)
GAG_AUC<-GAG_AUC[order(GAG_AUC$AUC,decreasing=TRUE),]
ggplot(GAG_AUC)+geom_bar(aes(x=thresholds,y=AUC),stat='identity')+theme(axis.text.x=element_text(angle=90))+scale_y_continuous(lim=c(-0.05,1),breaks=seq(0,1,l=21))+scale_x_discrete("Quantile Thresholds")+labs(title="GAG AUC Values")
#'Draw ROC for number 7, the top performer
i<-14
df<-cbind(truth<-rbind(results_paired[[i]]$test_data$GAG_treated[,c("PTID","VISITNO")],results_paired[[i]]$test_data$GAG_placebo[,c("PTID","VISITNO")]),
          P=c(results_paired[[i]]$classification_probs$GAG_treated,results_paired[[i]]$classification_probs$GAG_placebo),truth=c(rep(c("Treatment"),nrow(results_paired[[i]]$test_data$GAG_treated)),rep(c("Placebo"),nrow(results_paired[[i]]$test_data$GAG_placebo))))
df<-ddply(df,.(PTID),summarize,delta=1-diff(P[VISITNO%in%c("2","12")]),truth=unique(truth))
df<-df[order(df$delta,df$truth,decreasing=FALSE),]
o<-order(df$delta,decreasing=FALSE)
FPR<-cumsum(df$truth[o]=="Placebo")/sum(df$truth=="Placebo")
TPR<-cumsum(df$truth[o]=="Treatment")/sum(df$truth=="Treatment")
GAG_ROC<-data.frame(FPR,TPR)
ggplot(GAG_ROC)+geom_line(aes(x=FPR,y=TPR))

#'ENV
ENV_AUC<-sapply(seq_along(results_paired),function(i){
df<-cbind(truth<-rbind(results_paired[[i]]$test_data$ENV_treated[,c("PTID","VISITNO")],results_paired[[i]]$test_data$ENV_placebo[,c("PTID","VISITNO")]),
          P=c(results_paired[[i]]$classification_probs$ENV_treated,results_paired[[i]]$classification_probs$ENV_placebo),truth=c(rep(c("Treatment"),nrow(results_paired[[i]]$test_data$ENV_treated)),rep(c("Placebo"),nrow(results_paired[[i]]$test_data$ENV_placebo))))
df<-ddply(df,.(PTID),summarize,delta=1-diff(P[VISITNO%in%c("2","12")]),truth=unique(truth))
df<-df[order(df$delta,df$truth,decreasing=FALSE),]
o<-order(df$delta,decreasing=FALSE)
FPR<-cumsum(df$truth[o]=="Placebo")/sum(df$truth=="Placebo")
TPR<-cumsum(df$truth[o]=="Treatment")/sum(df$truth=="Treatment")
AUC<-trapz(FPR,TPR)
})
ENV_AUC<-data.frame(thresholds=names(results_paired),AUC=ENV_AUC)
ENV_AUC<-ENV_AUC[order(ENV_AUC$AUC,decreasing=TRUE),]
ggplot(ENV_AUC)+geom_bar(aes(x=thresholds,y=AUC),stat='identity')+theme(axis.text.x=element_text(angle=90))+scale_y_continuous(lim=c(-0.05,1),breaks=seq(0,1,l=21))+scale_x_discrete("Quantile Thresholds")+labs(title="ENV AUC Values")

#'ROC for top performer in ENV, 14
i<-14
df<-cbind(truth<-rbind(results_paired[[i]]$test_data$ENV_treated[,c("PTID","VISITNO")],results_paired[[i]]$test_data$ENV_placebo[,c("PTID","VISITNO")]),
          P=c(results_paired[[i]]$classification_probs$ENV_treated,results_paired[[i]]$classification_probs$ENV_placebo),truth=c(rep(c("Treatment"),nrow(results_paired[[i]]$test_data$ENV_treated)),rep(c("Placebo"),nrow(results_paired[[i]]$test_data$ENV_placebo))))
df<-ddply(df,.(PTID),summarize,delta=1-diff(P[VISITNO%in%c("2","12")]),truth=unique(truth))
df<-df[order(df$delta,df$truth,decreasing=FALSE),]
o<-order(df$delta,decreasing=FALSE)
FPR<-cumsum(df$truth[o]=="Placebo")/sum(df$truth=="Placebo")
TPR<-cumsum(df$truth[o]=="Treatment")/sum(df$truth=="Treatment")
ENV_ROC<-data.frame(FPR,TPR)
ggplot(ENV_ROC)+geom_line(aes(x=FPR,y=TPR))

#'Combined ROC for ENV and GAG
COMB_ROC<-data.frame(rbind(ENV_ROC,GAG_ROC),Stim=c(rep("ENV",nrow(ENV_ROC)),rep("GAG",nrow(GAG_ROC))))
ggplot(COMB_ROC)+geom_line(aes(x=FPR,y=TPR,color=Stim),lwd=2)+labs(title="ROC Curves for best gating of ENV and GAG Stimulated Samples")+theme(axis.title.x=element_text(size=21),axis.title.y=element_text(size=21))

ggplot(data.frame(rbind(ENV_AUC,GAG_AUC),Stim=rep(c("ENV","GAG"),each=nrow(ENV_AUC))))+geom_bar(aes(y=AUC,x=thresholds,fill=Stim),position="identity",stat="identity",alpha=0.5)+theme(axis.text.x=element_text(angle=90))+scale_y_continuous(breaks=seq(0,1,l=21))+labs(title="AUC for ENV and GAG Stimulations\nDifferent Gating Thresholds")
#End Greg's code


#+ ROC_by__top_AUC, eval=FALSE
p <- ggplot(z, aes(x = Placebo, y = Treated)) + geom_line()
p + xlab("FPR (Placebo)") + ylab("TPR (Vaccinated)")


#+ markers_by_top_AUC, eval=FALSE








# TODO: Make sure that this ordering matches the ordering of 'GAG_accuracy_threshold' combo
GAG_markers <- apply(cytokine_combinations, 1, function(combo) {
  paste0(c("TNFa", "IFNg", "IL2"), combo, collapse = ":")
})
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






accuracy_results <- cbind(cytokine_combinations, accuracy_results)

# TODO: Update the 'diff' scores to AUC 
# To ensure that the ordering of the quantile combinations is preserved,
# notice that we reverse the order of cytokines specified in 'ddply'.
accuracy_results <- ddply(accuracy_results, .(IL2, IFNg, TNFa), transform,
                          diff_GAG = GAG_treatment - GAG_placebo,
                          diff_ENV = ENV_treatment - ENV_placebo)

which_GAG_top <- order(accuracy_results$diff_GAG, decreasing = TRUE)[seq_top]
which_ENV_top <- order(accuracy_results$diff_ENV, decreasing = TRUE)[seq_top]

GAG_top <- accuracy_results[which_GAG_top, ]
ENV_top <- accuracy_results[which_ENV_top, ]

#' For the top 6 cytokine-quantile combinations from each stimulation group, we
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


#' ### Markers Selected by `glmnet` for Top 6 GAG Cytokine-Quantile Combinations
#+ markers_top6_GAG, results='asis'
print(xtable(GAG_top_markers, digits = 4), include.rownames = FALSE, type = "html")

#' ### Markers Selected by `glmnet` for Top 6 ENV Cytokine-Quantile Combinations
#+ markers_top6_ENV, results='asis'
print(xtable(ENV_top_markers, digits = 4), include.rownames = FALSE, type = "html")


















# TODO: Wait for full cytokine-gate information
#   On 11 March 2013, Greg said he would look into it.
#   After we have the full information, build a classifier for the manual proportions.
