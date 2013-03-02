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
selected_PTIDs <- sample(unique(pData_gs_manual$PTID), 4)

# To overcome some issues with NetCDF files, Mike suggested that I manually clone
# the CDF file before cloning it.
ncdf_flowSet <- ncdf_flowSet[which(pData_gs_manual$PTID %in% selected_PTIDs)]
ncdf_clone <- clone.ncdfFlowSet(ncdf_flowSet,
                                ncdfFile = file.path(archive_path, "hvtn065-subset.nc"),
                                isEmpty = FALSE)
trans <- estimateMedianLogicle(ncdf_clone, channels = colnames(ncdf_flowSet[[1]])[-c(1:3, 5)])

# Determines transformation for all channels except for "Time" and the sidescatter
# channels and then applies the estimated transformation to the cloned ncdfFlow
# set object.
ncdf_flowSet_trans <- transform(ncdf_clone, trans)

# Constructs a GatingSet object from the ncdfFlowSet object
gs_manual <- GatingSet(ncdf_flowSet_trans)

# Instead of loading the raw FCS files, we clone the manual gates into 3 new gating sets.
# One gating set per stimulation group: 1) negctrl, 2) ENV-1-PTEG, and 3) GAG-1-PTEG
pData_gs_manual <- subset(pData_gs_manual, PTID %in% selected_PTIDs)
gs_negctrl <- clone(gs_manual[which(pData_gs_manual$Stim == "negctrl")],
                    ncdfFile = file.path(archive_path, "hvtn065-negctrl.nc"))
gs_ENV <- clone(gs_manual[which(pData_gs_manual$Stim == "ENV-1-PTEG")],
                ncdfFile = file.path(archive_path, "hvtn065-ENV.nc"))
gs_GAG <- clone(gs_manual[which(pData_gs_manual$Stim == "GAG-1-PTEG")],
                ncdfFile = file.path(archive_path, "hvtn065-GAG.nc"))

# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData(gs_negctrl) <- ddply(pData(gs_negctrl), .(PTID, VISITNO), transform,
                           Stim = paste0(Stim, seq_along(Stim)))

# Now we apply the automated pipeline to each gating set and archive the results
# in the 'archive_path'.
gating_template <- new("HVTN065")
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

gating(gating_template, gs_negctrl, batch = TRUE, nslaves = 10)
archive(gs_negctrl, file = file.path(archive_path, "HVTN065-negctrl-42patients.tar"))

gating(gating_template, gs_ENV, batch = TRUE, nslaves = 10)
archive(gs_ENV, file = file.path(archive_path, "HVTN065-ENV-42patients.tar"))

gating(gating_template, gs_GAG, batch = TRUE, nslaves = 10)
archive(gs_GAG, file = file.path(archive_path, "HVTN065-GAG-42patients.tar"))

HVTN065_population_stats <- list()
HVTN065_population_stats$negctrl <- getPopStats(gs_negctrl)
HVTN065_population_stats$ENV <- getPopStats(gs_ENV)
HVTN065_population_stats$GAG <- getPopStats(gs_GAG)

HVTN065_pData_gs_manual <- rbind(pData(gs_negctrl), pData(gs_ENV), pData(gs_GAG))
HVTN065_pData_gs_manual <- subset(HVTN065_pData_gs_manual, select = c(name, PTID, Stim, VISITNO))

save(HVTN065_population_stats, HVTN065_pData_gs_manual, file = "data/HVTN065-results-42patients.RData")




gating(gating_template, gs_negctrl, batch = TRUE)

# gs_negctrl <- unarchive(file = file.path(archive_path, "HVTN065-negctrl-42patients.tar"), path = archive_path)
# gs_ENV <- unarchive(file = file.path(archive_path, "HVTN065-ENV-42patients.tar"), path = archive_path)
# gs_GAG <- unarchive(file = file.path(archive_path, "HVTN065-GAG-42patients.tar"), path = archive_path)



png("HVTN065-cytokine-unstimulated.png", width = 360, height = 360)
plotGate(gs_negctrl[37], 8:9, lattice = TRUE, xbin = 128, overlay = 13, cex = 2)
dev.off()

png("HVTN065-cytokine-stimulated.png", width = 360, height = 360)
plotGate(gs_ENV[17], 8:9, lattice = TRUE, xbin = 128, overlay = 13, cex = 2)
dev.off()
