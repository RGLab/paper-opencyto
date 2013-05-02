library(ProjectTemplate)
load.project()

num_cores <- 12

# Antigens (stimulation groups) that will be extracted
antigens <- c("ENV-1-PTEG", "GAG-1-PTEG", "sebctrl", "negctrl")

# The path of the root directory containing the FCS files
fcs_path <- '/shared/silo_researcher/Gottardo_R/gfinak_working/HVTN065/FACSData'

# The path of the compensation matrices
comp_path <- '/loc/no-backup/HVTN/065/compensation_matrices'

# Subdirectories containing FCS files
# BUG: The 'list.dirs' function returns only the full.names.
# We want only the folder name -- not the entire path.
fcs_subdirs <- list.dirs(fcs_path, full.names = FALSE, recursive = FALSE)
fcs_subdirs <- sapply(strsplit(fcs_subdirs, "/"), tail, n = 1)

# We traverse each of the subdirectories above and read the headers of each FCS
# file.  In particular, we care about the headers: 'EXPERIMENT NAME' and 'Sample
# Order' We extract these keywords and construct a data.frame summarizing both
# headers for all FCS files.
fcs_headers <- mclapply(fcs_subdirs, function(current_subdir) {
  subdir_fullpath <- file.path(fcs_path, current_subdir)
  fcs_subdir <- dir(subdir_fullpath, "fcs")

  # If no FCS files are in the directory, we return NULL.
  fcs_headers <- NULL
  if (length(fcs_subdir) > 0) {
    fcs_headers <- read.FCSheader(files = fcs_subdir, path = subdir_fullpath,
                                  keyword = c("EXPERIMENT NAME", "Sample Order", "Stim"))
    fcs_headers <- data.frame(do.call(rbind, fcs_headers))
    fcs_headers$FCS_file <- file.path(subdir_fullpath, rownames(fcs_headers))
    fcs_headers
  }
  fcs_headers
}, mc.cores = num_cores)

# Remove the NULL entries.
fcs_headers <- fcs_headers[!sapply(fcs_headers, is.null)]

# Some of the folders have no "Sample Order" keyword. In this case, a column of
# NA values is stored. In this case, we remove them.
fcs_headers <- do.call(rbind, lapply(fcs_headers, na.omit))

# Filter out columns of analysis.plan that we don't need
# We are only interested in mapping fcs files to PTIDs.
analysis_plan <- analysis.plan[, c("ASSAYID", "PTID", "SAMP_ORD", "ANTIGEN", "VISITNO")]
analysis_plan <- na.omit(unique(analysis_plan))
analysis_plan <- subset(analysis_plan, VISITNO %in% c("2", "12"))

# Converts data.frame's to data.table's for faster joins
analysis_plan <- data.table(analysis_plan)
fcs_headers <- data.table(fcs_headers)

# Merge 'fcs_headers' and 'analysis_plan' 
setkey(analysis_plan, SAMP_ORD, ASSAYID, ANTIGEN)
analysis_plan[, `:=`(PTID, factor(PTID))]
analysis_plan[, `:=`(SAMP_ORD, factor(SAMP_ORD))]
setnames(fcs_headers, c("ASSAYID", "SAMP_ORD", "ANTIGEN", "FCS_file"))
setkey(fcs_headers, SAMP_ORD, ASSAYID, ANTIGEN)
merged <- base:::merge(analysis_plan, fcs_headers, by = c("SAMP_ORD", "ASSAYID", "ANTIGEN"))
setkey(merged, PTID) 
merged_sub <- merged[ANTIGEN %in% antigens]

# Below, we wish to ensure that there are two negative controls for each PTID.
# To do this, we update rows with ANTIGEN of 'negctrl' to 'negctrl1' and
# 'negctrl2'.
merged_negctrl <- merged_sub[ANTIGEN == "negctrl"]
merged_negctrl <- ddply(merged_negctrl, .(PTID, VISITNO, SAMP_ORD), transform,
                       ANTIGEN = paste0(ANTIGEN, seq_along(ANTIGEN)))
merged_sub[merged_sub$ANTIGEN == "negctrl", ] <- merged_negctrl
merged_sub <- merged_sub[!(ANTIGEN %in% c("negctrl3", "negctrl4"))]
setkey(merged_sub, PTID, ANTIGEN, VISITNO)
merged_sub <- unique(merged_sub)

# Now, we replace "negctrl1" and "negctrl2" with "negctrl" to ensure they are
# grouped together.
merged_sub$ANTIGEN <- gsub("negctrl[12]", "negctrl", as.character(merged_sub$ANTIGEN))
merged_sub$ANTIGEN <- factor(merged_sub$ANTIGEN)

# Grab filenames for compensation matrices from 'comp_path'
# Simply append '-comp' to the ASSAYID and extract
# We remove any patients for which no file with a compensation matrix exists.
merged_sub$comp_file <- file.path(comp_path, paste0(merged_sub$ASSAYID, "-comp"))
merged_sub <- merged_sub[file.exists(comp_file)]

# Remove records with only one visit
merged_sub <- merged_sub[, .SD[length(unique(VISITNO)) == 2], by = c("PTID,ANTIGEN")]

# Create analysis_plan
analysis_plan <- merged_sub
analysis_plan <- as.data.frame(lapply(analysis_plan, factor)) 
analysis_plan$FCS_file <- as.character(analysis_plan$FCS_file)
analysis_plan$comp_file <- as.character(analysis_plan$comp_file)

# For each unique compensation file, we read in the FCS files are to be
# compensated with this matrix and then construct a flowSet object.
fs_list <- tapply(seq_along(analysis_plan$comp_file), analysis_plan$comp_file, function(i) {
  comp_file <- unique(analysis_plan$comp_file[i])

  # Reads the compensation matrix and the column names separately to avoid the
  # R-friendly dots that are added to some of the channels
  comp_colnames <- as.character(read.csv(comp_file, header = FALSE,
                                         stringsAsFactors = FALSE, nrows = 1))
  comp_matrix <- read.csv(comp_file, header = FALSE, stringsAsFactors = FALSE,
                          skip = 1)
  colnames(comp_matrix) <- comp_colnames

  # Reads the FCS files as a flowSet and then apply compensation
  flow_set <- read.flowSet(analysis_plan$FCS_file[i], min.limit = -100)
  flow_set <- compensate(flow_set, comp_matrix)

  flow_set
})

# Across all samples in the flowSet, we find the common channel names.
common_channels <- colnames(fs_list[[1]])
for (i in seq.int(2, length(fs_list))) {
  common_channels <- intersect(common_channels, colnames(fs_list[[i]]))
}

# Loop through the 'fs_list' and extract the common_channels
fs_list <- lapply(seq_along(fs_list), function(i) {
  fsApply(fs_list[[i]], function(x) {
    x[, common_channels]
  })
})

# Merges the list of flowSet objects into a single flowSet object. This code is
# verbose but it circumvents an issue introduced recently in flowIncubator.
flow_set <- fs_list[[1]]
for (i in seq.int(2, length(fs_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]])
}

# Determines transformation for all channels except for "Time" and the sidescatter channels
logicle_trans <- openCyto:::estimateMedianLogicle(flow_set,
                                                  channels = common_channels[-(1:4)])

# Applies the estimated transformation to the ncdfFlow set object
flow_set <- transform(flow_set, logicle_trans)

# Updates pData and varMetadata
pData_HVTN065 <- subset(analysis_plan, select = c(FCS_file, PTID, ANTIGEN, VISITNO))
colnames(pData_HVTN065) <- c("name", "PTID", "Stim", "VISITNO")
pData_HVTN065$name <- sapply(strsplit(pData_HVTN065$name, "/"), tail, n = 1)
pData_HVTN065 <- base:::merge(pData(flow_set), pData_HVTN065, sort = FALSE)
rownames(pData_HVTN065) <- pData_HVTN065$name
pData(flow_set) <- pData_HVTN065

var_meta <- varMetadata(phenoData(flow_set))
var_meta[-1, ] <- rownames(var_meta)[-1]
varMetadata(phenoData(flow_set)) <- var_meta

# Create ncdfFlowSet object from flowSet
ncdf_flowset <- ncdfFlowSet(flow_set, ncdfFile = "/loc/no-backup/ramey/HVTN/065/HVTN065.nc")

# Create GatingSet from ncdfFlowSet
gs_HVTN065 <- GatingSet(ncdf_flowset)

# Archives GatingSet
# save_gs(gs_HVTN065, path = "/loc/no-backup/ramey/HVTN/065/gating-set")
save_gs(gs_HVTN065, path = "/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating_set")
