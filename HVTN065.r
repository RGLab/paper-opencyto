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
selected_PTIDs <- sample(unique(pData_fs$PTID), 20)

# To overcome some issues with NetCDF files, Mike suggested that I manually clone
# the CDF file before transforming the samples.
ncdf_flowSet <- ncdf_flowSet[which(pData_fs$PTID %in% selected_PTIDs)]
ncdf_clone <- clone.ncdfFlowSet(ncdf_flowSet, ncdfFile = file.path(archive_path, "hvtn065-subset.nc"),
                                isEmpty = FALSE)

# Determines transformation for all channels except for "Time" and the sidescatter channels
trans <- openCyto:::estimateMedianLogicle(ncdf_flowSet, channels = colnames(ncdf_flowSet)[-c(1:3, 5)])

# Applies the estimated transformation to the ncdfFlow set object
ncdf_flowSet_trans <- transform(ncdf_flowSet, trans)

# Constructs a GatingSet object from the ncdfFlowSet object
gs_manual <- GatingSet(ncdf_flowSet_trans)

# Instead of loading the raw FCS files, we clone the manual gates into 3 new gating sets.
# One gating set per stimulation group: 1) negctrl, 2) ENV-1-PTEG, and 3) GAG-1-PTEG
pData_fs <- subset(pData_fs, PTID %in% selected_PTIDs)
gs_negctrl <- flowWorkspace:::clone(gs_manual[grep("negctrl", pData_fs$Stim)])
gs_ENV <- flowWorkspace:::clone(gs_manual[grep("^ENV", pData_fs$Stim)])
gs_GAG <- flowWorkspace:::clone(gs_manual[grep("^GAG", pData_fs$Stim)])

# Loads the GatingTemplate from the CSV file.
gating_template <- gatingTemplate("HVTN065-GatingTemplate.csv", "HVTN065")

# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData(gs_negctrl) <- ddply(pData(gs_negctrl), .(PTID, VISITNO), transform,
                           Stim = paste0(Stim, seq_along(Stim)))

# Now we apply the automated pipeline to each gating set and archive the results
message("Negative Controls")
gating(gating_template, gs_negctrl, num_nodes = 10, parallel_type = "multicore")
message("ENV Stimulated")
gating(gating_template, gs_ENV, num_nodes = 10, parallel_type = "multicore")
message("GAG Stimulated")
gating(gating_template, gs_GAG, num_nodes = 10, parallel_type = "multicore")

archive(gs_negctrl, file = file.path(archive_path, "HVTN065-negctrl.tar"))
archive(gs_ENV, file = file.path(archive_path, "HVTN065-ENV.tar"))
archive(gs_GAG, file = file.path(archive_path, "HVTN065-GAG.tar"))

HVTN065_popstats <- list()
HVTN065_popstats$negctrl <- getPopStats(gs_negctrl)
HVTN065_popstats$ENV <- getPopStats(gs_ENV)
HVTN065_popstats$GAG <- getPopStats(gs_GAG)

HVTN065_pData <- rbind(pData(gs_negctrl), pData(gs_ENV), pData(gs_GAG))
HVTN065_pData <- subset(HVTN065_pData, select = c(name, PTID, Stim, VISITNO))

save(HVTN065_popstats, HVTN065_pData, file = "data/HVTN065-results.RData")



       
