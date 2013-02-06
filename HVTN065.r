library(ProjectTemplate)
load.project()

# Loads the necessary pipeline code
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

# Per Mike: "We have 62 subjects *2 visit *(1 env+1gag +2 negctrls)=496 samples (don't have IL4 in this data set)"
ncdf_file <- '/loc/no-backup/ramey/hvtn065.nc'

ncdf_flowSet <- ncdfFlowSet_open(ncdf_file)
ncdf_flowSet@file <- ncdf_file

# Customizes the path where the GS is unarchived.
# By default, /tmp is used and often can run out of space, due to the large nc files.
archive_path <- '/loc/no-backup/ramey'

# The 'pData' function has some extra info parsed from thawlist and FCS keywords,
# but we only need to look at these 4 columns:
#  name, PTID, Stim, VISITNO
pData_gs_manual <- subset(pData(ncdf_flowSet), select = c(name, PTID, Stim, VISITNO))

# TEMP: For the moment, we randomly select 3 patient IDs and gate their samples.
# TODO: Remove this entire code block to omit the random sampling
set.seed(42)
selected_PTIDs <- sample(unique(pData_gs_manual$PTID), 8)
ncdf_flowSet <- ncdf_flowSet[which(pData_gs_manual$PTID %in% selected_PTIDs)]
pData_gs_manual <- subset(pData_gs_manual, PTID %in% selected_PTIDs)

# Determines transformation for all channels except for "Time" and the sidescatter channels
trans <- estimateMedianLogicle(ncdf_flowSet, channels = colnames(ncdf_flowSet[[1]])[-c(1:3, 5)])

# Applies the estimated transformation to the ncdfFlow set object
ncdf_flowSet_trans <- transform(ncdf_flowSet, trans)

# Constructs a GatingSet object from the ncdfFlowSet object
gs_manual <- GatingSet(ncdf_flowSet_trans)

# Instead of loading the raw FCS files, we clone the manual gates into 3 new gating sets.
# One gating set per stimulation group: 1) negctrl, 2) ENV-1-PTEG, and 3) GAG-1-PTEG
gs_negctrl <- clone(gs_manual[which(pData_gs_manual$Stim == "negctrl")])
gs_ENV <- clone(gs_manual[which(pData_gs_manual$Stim == "ENV-1-PTEG")])
gs_GAG <- clone(gs_manual[which(pData_gs_manual$Stim == "GAG-1-PTEG")])

# Now we apply the automated pipeline to each gating set and archive the results
# in the 'archive_path'.
gating_template <- new("HVTN065")
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

gating(gating_template, gs_negctrl, batch = TRUE, nslaves = 9)
archive(gs_negctrl, file = file.path(archive_path, "HVTN065-negctrl.tar"))

gating(gating_template, gs_ENV, batch = TRUE) #, nslaves = 9)
archive(gs_ENV, file = file.path(archive_path, "HVTN065-ENV.tar"))

gating(gating_template, gs_GAG, batch = TRUE, nslaves = 9)
archive(gs_GAG, file = file.path(archive_path, "HVTN065-GAG.tar"))

population_stats <- list()
population_stats$negctrl <- pretty_popstats(getPopStats(gs_negctrl))
population_stats$ENV <- pretty_popstats(getPopStats(gs_ENV))
population_stats$GAG <- pretty_popstats(getPopStats(gs_GAG))

save(population_stats, pData_gs_manual, file = "data/HVTN065-results.RData")




debug(flowWorkspace:::.addGate)
