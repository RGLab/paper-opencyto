library(ProjectTemplate)
load.project()

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
selected_PTIDs <- sample(unique(pData_gs_manual$PTID), 2)
ncdf_flowSet <- ncdf_flowSet[which(pData_gs_manual$PTID %in% selected_PTIDs)]
pData_gs_manual <- subset(pData_gs_manual, PTID %in% selected_PTIDs)

# Determines transformation for all channels except for "Time" and the sidescatter channels
trans <- estimateMedianLogicle(ncdf_flowSet, channels = colnames(ncdf_flowSet)[-c(1:3, 5)])

# Applies the estimated transformation to the ncdfFlow set object
ncdf_flowSet_trans <- transform(ncdf_flowSet, trans)

# Constructs a GatingSet object from the ncdfFlowSet object
gs_manual <- GatingSet(ncdf_flowSet_trans)

gating_template <- gatingTemplate("HVTN065-GatingTemplate.csv", "HVTN065")

gating(gating_template, gs_manual, num_nodes = 2, parallel_type = "multicore")


# Instead of loading the raw FCS files, we clone the manual gates into 3 new gating sets.
# One gating set per stimulation group: 1) negctrl, 2) ENV-1-PTEG, and 3) GAG-1-PTEG
gs_negctrl <- clone(gs_manual[which(pData_gs_manual$Stim == "negctrl")])
gs_ENV <- clone(gs_manual[which(pData_gs_manual$Stim == "ENV-1-PTEG")])
gs_GAG <- clone(gs_manual[which(pData_gs_manual$Stim == "GAG-1-PTEG")])

# Now we apply the automated pipeline to each gating set and archive the results
# in the 'archive_path'.
gating_template <- gatingTemplate("HVTN065-GatingTemplate.csv", "HVTN065")

gating(gating_template, gs_manual)

gating(gating_template, gs_negctrl, batch = TRUE, nslaves = 9)
# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData(gs_negctrl) <- ddply(pData(gs_negctrl), .(PTID, VISITNO), transform,
                           Stim = paste0(Stim, seq_along(Stim)))
archive(gs_negctrl, file = file.path(archive_path, "HVTN065-negctrl.tar"))

gating(gating_template, gs_ENV, batch = TRUE)
archive(gs_ENV, file = file.path(archive_path, "HVTN065-ENV.tar"))

gating(gating_template, gs_GAG, batch = TRUE, nslaves = 9)
archive(gs_GAG, file = file.path(archive_path, "HVTN065-GAG.tar"))

HVTN065_population_stats <- list()
HVTN065_population_stats$negctrl <- getPopStats(gs_negctrl)
HVTN065_population_stats$ENV <- getPopStats(gs_ENV)
HVTN065_population_stats$GAG <- getPopStats(gs_GAG)

HVTN065_pData_gs_manual <- rbind(pData(gs_negctrl), pData(gs_ENV), pData(gs_GAG))
HVTN065_pData_gs_manual <- subset(HVTN065_pData_gs_manual, select = c(name, PTID, Stim, VISITNO))

save(HVTN065_population_stats, HVTN065_pData_gs_manual, file = "data/HVTN065-results.RData")



       
