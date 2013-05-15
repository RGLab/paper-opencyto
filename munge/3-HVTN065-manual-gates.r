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
manual_data <- subset(manual_data, PTID %in% levels(pData_HVTN065$PTID))

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
  stimulated_pct <- x$pctpos
  names(stimulated_pct) <- x$Marker

  # In the same row of 'manual_data' the negative-control percentages are also
  # given. We extract those here.
  negctrl_pct <- x$pctpos_neg
  names(negctrl_pct) <- x$Marker

  proportions <- data.frame(PTID, VISITNO, Stimulation = c(Stimulation, "negctrl"),
                            rbind(stimulated_pct, negctrl_pct),
                            check.names = FALSE)
  rownames(proportions) <- NULL

  proportions
})

save(HVTN065_manual, file = "cache/HVTN065-manual-proportions.RData")

