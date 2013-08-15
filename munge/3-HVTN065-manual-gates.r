# Retrieved data from Greg Finak at /home/gfinak/065.rds

# Some comments from Greg via email on 14 May 2013:
# The relevant columns are ptid, tcellsub, visitno, antigen, response,
# filt_flag, nsub, cytnum, nsub_ref, cytnum_ref (or maybe nsub_neg and
# cytnum_neg). Most are self explanatory, nsub and cytnum are the negative and
# positive cell counts for stimulated samples, nsub_ref and cynum_ref are the
# negative and positive cell counts for unstimulated samples. response is the
# scharp response call, filt_flag indicates if the sample was filtered out for
# various reasons (high background or low T-cell counts).

library(ProjectTemplate)
load.project()

manual_data <- readRDS("data/065.rds")

# Antigens (stimulation groups) that will be extracted
antigens <- c("ENV-1-PTEG", "GAG-1-PTEG")

# Adds additional columns that are formatted for easier usage downstream.
manual_data$VISITNO <- as.character(manual_data$visitno)
manual_data$PTID <- gsub("-", "", as.character(manual_data$ptid))
manual_data$Stimulation <- as.character(manual_data$antigen)
manual_data <- subset(manual_data, PTID %in% as.character(unique(pData_HVTN065$PTID)))

manual_data$Marker <- sapply(strsplit(as.character(manual_data$tcellsub), "/"), tail, n = 1)
manual_data$Marker <- gsub("\\+", "", manual_data$Marker)
manual_data$Marker <- paste0(manual_data$Marker, ":IL2+|IFNg+")

# We filter out samples marked with a filtering flag (i.e., 'filt_flag').
# To so, we determine which PTIDs have 'filt_flag' set to 1 and then remove those PTIDs.
filter_PTIDs <- unique(subset(manual_data, filt_flag == 1)$PTID)
manual_data <- manual_data[!(manual_data$PTID %in% filter_PTIDs), ]

# Construct data.frame of proportions derived from manual gates
# Format: PTID, VISITNO, Stimulation, CD4:IL2+|IFNg+, CD8:IL2+|IFNg+
HVTN065_manual <- ddply(manual_data, .(PTID, VISITNO, Stimulation), function(x) {
  PTID <- x$PTID[1]
  Stimulation <- x$Stimulation[1]
  VISITNO <- x$VISITNO[1]

  # First, we obtain the percentages for the stimulation group.
  stimulated_pct <- with(x, cytnum / nsub)
  names(stimulated_pct) <- x$Marker

  # In the same row of 'manual_data' the negative-control percentages are also
  # given. We extract those here.
  negctrl_pct <- with(x, cytnum_neg / nsub_neg)
  names(negctrl_pct) <- x$Marker

  proportions <- data.frame(PTID, VISITNO, Stimulation = c(Stimulation, "negctrl"),
                            rbind(stimulated_pct, negctrl_pct),
                            check.names = FALSE)
  rownames(proportions) <- NULL

  proportions
})

# Because the negative controls are repeated twice given our construct above, we
# average the proportions within PTID, VISITNO, and Stimulation group.
HVTN065_manual <- ddply(HVTN065_manual, .(PTID, VISITNO, Stimulation), summarize,
                        `cd4:IL2+|IFNg+` = mean(`cd4:IL2+|IFNg+`),
                        `cd8:IL2+|IFNg+` = mean(`cd8:IL2+|IFNg+`))

# Ensures that 2 precedes 12 in the factor level ordering for VISITNO
HVTN065_manual$VISITNO <- relevel(HVTN065_manual$VISITNO, "2")

# Relevels Stimulation group so that 'negctrl' appears first because it is the
# baseline level.
HVTN065_manual$Stimulation <- relevel(HVTN065_manual$Stimulation, "negctrl")

save(HVTN065_manual, file = "cache/HVTN065-manual-proportions.RData")

# NOTE: It is not clear to me how to compute proportions for upstream gates if
# we wish to do this.  In the 'manual_data' data.frame above, each row gives
# both the counts and proportions for the cytokines for both the stimulated
# sample and negative control. However, each row has only one of each upstream
# gate, so it is not clear if this should correspond to stimulation or negative
# control.
