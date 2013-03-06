library(ProjectTemplate)
load.project()

# Per Mike: "We have 62 subjects *2 visit *(1 env+1gag +2 negctrls)=496 samples (don't have IL4 in this data set)"
ncdf_file <- '/loc/no-backup/ramey/HVTN/065/hvtn065.nc'

ncdf_flowSet <- ncdfFlowSet_open(ncdf_file)
ncdf_flowSet@file <- ncdf_file

# Customizes the path where the GS is unarchived.
# By default, /tmp is used and often can run out of space, due to the large nc files.
archive_path <- '/loc/no-backup/ramey/HVTN/065/'

# The 'pData' function has some extra info parsed from thawlist and FCS keywords,
# but we only need to look at these 4 columns:
#  name, PTID, Stim, VISITNO
pData_fs <- subset(pData(ncdf_flowSet), select = c(name, PTID, Stim, VISITNO))

# TEMP: For the moment, we randomly select 4 patient IDs and gate their samples.
# TODO: Remove this entire code block to omit the random sampling
set.seed(42)
selected_PTIDs <- sample(unique(pData_fs$PTID), 30)

# To overcome some issues with NetCDF files, Mike suggested that I manually clone
# the CDF file before transforming the samples.
ncdf_flowSet <- ncdf_flowSet[which(pData_fs$PTID %in% selected_PTIDs)]
ncdf_clone <- clone.ncdfFlowSet(ncdf_flowSet, ncdfFile = file.path(archive_path, "hvtn065-subset.nc"),
                                isEmpty = FALSE)
pData_fs <- subset(pData_fs, PTID %in% selected_PTIDs)

# Determines transformation for all channels except for "Time" and the sidescatter channels
trans <- openCyto:::estimateMedianLogicle(ncdf_clone, channels = colnames(ncdf_flowSet)[-c(1:3, 5)])

# Applies the estimated transformation to the ncdfFlow set object
ncdf_flowSet_trans <- transform(ncdf_clone, trans)

# Constructs a GatingSet object from the transformed ncdfFlowSet object
gs_HVTN065 <- GatingSet(ncdf_flowSet_trans)

# Loads the GatingTemplate from the CSV file.
gating_template <- gatingTemplate("HVTN065-GatingTemplate.csv", "HVTN065")

# Now we apply the automated pipeline to each gating set
gating(gating_template, gs_HVTN065, prior_group = 'Stim', num_nodes = 12, parallel_type = "multicore")

# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData_HVTN065 <- pData(gs_HVTN065)
pData_HVTN065$Stim <- as.character(pData_HVTN065$Stim)
pData_negctrl <- subset(pData_HVTN065, Stim == "negctrl")
pData_negctrl <- ddply(pData_negctrl, .(PTID, VISITNO), transform,
                       Stim = paste0(Stim, seq_along(Stim)))
pData_HVTN065[pData_HVTN065$Stim == "negctrl", ] <- pData_negctrl
pData_HVTN065 <- subset(pData_HVTN065, select = c(name, PTID, Stim, VISITNO))
pData(gs_HVTN065) <- pData_HVTN065

popstats_HVTN065 <- getPopStats(gs_HVTN065)

save(popstats_HVTN065, pData_HVTN065, file = "data/HVTN065-results.RData")

# Archives the results
archive(gs_HVTN065, file = file.path(archive_path, "HVTN065.tar"))
