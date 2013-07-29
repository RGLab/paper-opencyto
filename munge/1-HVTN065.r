library(ProjectTemplate)
load.project()

options(stringsAsFactors = FALSE)
num_cores <- 10

# Antigens (stimulation groups) that will be extracted
antigens <- c("ENV-1-PTEG", "GAG-1-PTEG", "POL-1-PTEG", "negctrl")

# Minimum number of cells per FCS file
# We require the number of cells to be at least 100,000
min_cells <- 1e5

# Constructs a lookup table of the relevant channels and markers
channels <- c("Time", "FSC-A", "FSC-H", "SSC-A", "FITC-A", "Pacific Blue-A",
              "Alexa 680-A", "PE Green laser-A", "PE Tx RD-A", "PE Cy55-A",
              "PE Cy7-A")
markers <- c(NA, NA, NA, NA, "CD4", "ViViD", "TNFa", "IL2", "CD3", "CD8", "IFNg")
marker_lookup <- data.frame(channels, markers, stringsAsFactors = FALSE)

# The path of the root directory containing the FCS files
fcs_path <- '/shared/silo_researcher/Gottardo_R/gfinak_working/HVTN065/FACSData'

# The path of the compensation matrices
comp_path <- '/loc/no-backup/HVTN/065/compensation_matrices'

# The root path of where GatingSets and temp storage
root_path <- '/loc/no-backup/ramey/HVTN/065/'
temp_path <- file.path(root_path, 'temp')
dir.create(temp_path)

# Subdirectories containing FCS files
# BUG: The 'list.dirs' function returns only the full.names.
# We want only the folder name -- not the entire path.
fcs_subdirs <- list.dirs(fcs_path, full.names = FALSE, recursive = FALSE)
fcs_subdirs <- sapply(strsplit(fcs_subdirs, "/"), tail, n = 1)

# We traverse each of the subdirectories above and read the headers of each FCS
# file. In particular, we care about the headers: 'EXPERIMENT NAME', 'Sample
# Order', 'Stim', and "$TOT". We extract these keywords and construct a
# data.frame summarizing both headers for all FCS files.  Some of the folders
# have no "Sample Order" keyword. In this case, we remove them.
#
# Also, we extract a vector contaning the unique channel names, which are stored
# in keywords of the form $P##N.  For example, $P9N or $P11N.
message("Reading headers...")
fcs_headers <- mclapply(fcs_subdirs, function(current_subdir) {
  subdir_fullpath <- file.path(fcs_path, current_subdir)
  fcs_subdir <- dir(subdir_fullpath, "fcs")

  # If no FCS files are in the directory, we return NULL.
  fcs_headers <- NULL
  if (length(fcs_subdir) > 0) {
    fcs_headers <- read.FCSheader(files = fcs_subdir, path = subdir_fullpath,
                                  keyword = c("EXPERIMENT NAME", "Sample Order", "Stim", "$TOT"))
    fcs_headers <- data.frame(do.call(rbind, fcs_headers), check.names = FALSE)
    fcs_headers[["$TOT"]] <- as.numeric(fcs_headers[["$TOT"]])
    fcs_headers$FCS_file <- file.path(subdir_fullpath, rownames(fcs_headers))
  }
  fcs_headers
}, mc.cores = num_cores)

# Remove the NULL entries.
fcs_headers <- fcs_headers[!sapply(fcs_headers, is.null)]

# Remove entries that have no "Sample Order" and then combines the remaining
# headers into a single data.frame.
has_sample_order <- sapply(fcs_headers, function(x) {
  "Sample Order" %in% colnames(x)
})
fcs_headers <- fcs_headers[has_sample_order]
fcs_headers <- data.table(do.call(rbind, fcs_headers))
fcs_headers <- na.omit(fcs_headers)
setnames(fcs_headers, c("ASSAYID", "SAMP_ORD", "ANTIGEN", "num_cells", "FCS_file"))
setkey(fcs_headers, SAMP_ORD, ASSAYID, ANTIGEN)

# Before we map FCS files to PTIDS, we remove the samples that have comments.
# In particular, we want to remove the samples that have the following comment:
# (File:)Unreliable staining, omit from analysis
message("Constructing the analysis plan...")
analysis_plan <- subset(analysis.plan, Comments == "")

# Filter out columns of analysis.plan that we don't need
# We are only interested in mapping fcs files to PTIDs.
analysis_plan <- analysis_plan[, c("ASSAYID", "PTID", "SAMP_ORD", "ANTIGEN", "VISITNO")]
analysis_plan <- na.omit(unique(analysis_plan))
analysis_plan <- subset(analysis_plan, VISITNO %in% c("2", "12"))

# Converts data.frame's to data.table's for faster joins
analysis_plan <- data.table(analysis_plan)
setkey(analysis_plan, SAMP_ORD, ASSAYID, ANTIGEN)
analysis_plan[, `:=`(PTID, factor(PTID))]
analysis_plan[, `:=`(SAMP_ORD, factor(SAMP_ORD))]

# Merge 'fcs_headers' and 'analysis_plan'
message("Merging FCS headers and analysis plan...")
merged <- base:::merge(analysis_plan, fcs_headers, by = c("SAMP_ORD", "ASSAYID", "ANTIGEN"))
setkey(merged, PTID) 
merged_sub <- merged[ANTIGEN %in% antigens]
merged_sub <- merged_sub[num_cells >= min_cells]

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

# We ensure that there are 2 negative controls and 1 sample for each of the
# stimulation groups for each time point.
analysis_plan <- subset(analysis_plan, PTID %in% names(which(table(analysis_plan$PTID) == 10)))

# Saves analysis_plan to temp directory
# This is useful for caching the first part of this script in case the code
# below fails from memory issues
message("Caching analysis plan in temp folder...")
write.csv(analysis_plan, file = file.path(temp_path, "analysis_plan.csv"),
          row.names = FALSE)

# We remove the following objects to free some memory for the flowSet
# construction below
rm(fcs_headers)
rm(analysis.plan)
rm(merged)
rm(merged_sub)
rm(merged_negctrl)
gc(reset = TRUE)

# For each unique compensation file, we read in the FCS files are to be
# compensated with this matrix and then construct a flowSet object.
message("Constructing flowSet objects from FCS files and then compensating...")
fs_list <- lapply(unique(analysis_plan$comp_file), function(comp_file) {
  # Reads the compensation matrix and the column names separately to avoid the
  # R-friendly dots that are added to some of the channels
  comp_colnames <- as.character(read.csv(comp_file, header = FALSE,
                                         stringsAsFactors = FALSE, nrows = 1))
  comp_matrix <- read.csv(comp_file, header = FALSE, stringsAsFactors = FALSE,
                          skip = 1)
  colnames(comp_matrix) <- comp_colnames

  # Reads the FCS files as a flowSet and then compensates
  fcs_files <- analysis_plan$FCS_file[analysis_plan$comp_file == comp_file]
  compensate(read.flowSet(fcs_files, min.limit = -100), comp_matrix)
})

# Combines flowSet objects in list into a single flowSet. We reduce the data to
# a set of relevant channels beforehand.
message("Extracting relevant channels and then constructing a single flowSet object...")
flow_set <- fs_list[[1]][, marker_lookup$channels]
for (i in seq.int(2, length(fs_list))) {
  flow_set <- rbind2(flow_set, fs_list[[i]][, marker_lookup$channels])
}

# Removes the list of large flowSet objects from memory before transforming.
rm(fs_list)
gc(reset = TRUE)

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

message("Caching pData in temp folder...")
write.csv(pData_HVTN065, file = file.path(temp_path, "pData_HVTN065.csv"),
          row.names = FALSE)

write.csv(var_meta, file = file.path(temp_path, "var_meta.csv"))


# Determines transformation for all channels except for "Time" and the sidescatter channels
message("Estimating transformation...")
logicle_trans <- flowCore:::estimateMedianLogicle(flow_set,
                 channels = setdiff(marker_lookup$channels, c("Time", "FSC-A", "FSC-H", "SSC-A")))
flow_set <- transform(flow_set, logicle_trans)

# Saves the flow_set to 'temp_path'
write.flowSet(flow_set, outdir = temp_path)

# Create ncdfFlowSet object from flowSet
# Reads a ncdfFlowSet from the cached FCS files to avoid memory-allocation errors
message("Creating ncdfFlowSet...")
# ncdf_flowset <- ncdfFlowSet(flow_set, ncdfFile = "/loc/no-backup/ramey/HVTN/065/HVTN065.nc")
fcs_files <- dir(path = temp_path, pattern = ".fcs", full.names = TRUE)
ncdf_flowset <- read.ncdfFlowSet(files = fcs_files,
                                 ncdfFile = "/loc/no-backup/ramey/HVTN/065/HVTN065.nc")

# Updates pData, var_meta, and parameters manually
# I tried to do this with read.ncdfFlowSet, but nope...error, as always.
# I posted the following GitHub issue regarding the problem:
# https://github.com/RGLab/ncdfFlow/issues/18
pData_HVTN065 <- read.csv(file.path(temp_path, "pData_HVTN065.csv"))
pData_HVTN065 <- merge(pData(ncdf_flowset), pData_HVTN065)
rownames(pData_HVTN065) <- pData_HVTN065$name
pData(ncdf_flowset) <- pData_HVTN065

var_meta <- read.csv(file.path(temp_path, "var_meta.csv"), row.names = 1)
varMetadata(phenoData(ncdf_flowset)) <- var_meta

# Manually updates parameters. Note that this can take a couple of hours. Sigh.
for (i in seq_along(ncdf_flowset)) {
  temp_pData <- pData(parameters(ncdf_flowset[[i]]))
  channel_match <- match(marker_lookup$channels, temp_pData$name)
  temp_pData$desc <- marker_lookup$markers[channel_match]
  pData(parameters(ncdf_flowset[[i]])) <- temp_pData
}


# Removes the list of large flowSet objects from memory before transforming.
rm(flow_set)
gc(reset = TRUE)

# Create GatingSet from ncdfFlowSet
message("Constructing GatingSet...")
gs_HVTN065 <- GatingSet(ncdf_flowset)

# Archives GatingSet
message("Archiving GatingSet...")
save_gs(gs_HVTN065, path = file.path(root_path, "gating-set"))

unlink(temp_path, recursive = TRUE)
