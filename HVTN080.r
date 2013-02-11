library(ProjectTemplate)
load.project()

# Loads the necessary pipeline code
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

# Gating Sets (GSs) parsed from 368 manually gated samples (xml). 184 per/GS.
# We have two GSs because each consists of slightly different channels.
# We want to run the ICS pipeline on each separately.
gatingset_file <- '/loc/no-backup/HVTN/080/output/gs1.tar'
# gatingset_file <- '/loc/no-backup/HVTN/080/output/gs2.tar'

# Customizes the path where the GS is unarchived.
# By default, /tmp is used and often can run out of space, due to the large nc files.
archive_path <- '/loc/no-backup/ramey'

# Restore the GS from the tar file specified in 'gatingset_file'
gs_manual <- unarchive(gatingset_file, path = archive_path)

# The 'pData' function has some extra info parsed from thawlist and FCS keywords,
# but we only need to look at these 4 columns:
#  name, PTID, Stim, VISITNO
# NOTE: Samples for visit #2 were not included for the patient with ID of 158460256.
# TEMP: Because of this, we remove it from consideration for now.
pData_gs_manual <- subset(pData(gs_manual), select = c(name, PTID, Stim, VISITNO))
pData_gs_manual <- subset(pData_gs_manual, PTID != 158460256)

# TEMP: For the moment, we randomly select 3 patient IDs and gate their samples.
# TODO: Remove this entire code block to omit the random sampling
set.seed(42)
selected_PTIDs <- sample(unique(pData_gs_manual$PTID), 8)
gs_manual <- gs_manual[which(pData_gs_manual$PTID %in% selected_PTIDs)]
pData_gs_manual <- subset(pData_gs_manual, PTID %in% selected_PTIDs)

# Instead of loading the raw FCS files, we clone the manual gates into 3 new gating sets.
# One gating set per stimulation group: 1) negctrl, 2) ENV-1-PTEG, and 3) GAG-1-PTEG
gs_negctrl <- clone(gs_manual[which(pData_gs_manual$Stim == "negctrl")])
gs_ENV <- clone(gs_manual[which(pData_gs_manual$Stim == "ENV-1-PTEG")])
gs_GAG <- clone(gs_manual[which(pData_gs_manual$Stim == "GAG-1-PTEG")])

# Then, we remove the manual gates from each of the cloned gating sets.
Rm("S", gs_negctrl); Rm("S", gs_ENV); Rm("S", gs_GAG)

# Now we apply the automated pipeline to each gating set and archive the results
# in the 'archive_path'.
gating_template <- new("HVTN080")
lapply(list.files("~/rglab/HIMCLyoplate/Gottardo/pipeline/R", full = TRUE), source)

gating(gating_template, gs_negctrl, batch = TRUE, nslaves = 9)
# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData(gs_negctrl) <- ddply(pData(gs_negctrl), .(PTID, VISITNO), transform,
                           Stim = paste0(Stim, seq_along(Stim)))
archive(gs_negctrl, file = file.path(archive_path, "HVTN080-negctrl.tar"))

gating(gating_template, gs_ENV, batch = TRUE)
archive(gs_ENV, file = file.path(archive_path, "HVTN080-ENV.tar"))

gating(gating_template, gs_GAG, batch = TRUE, nslaves = 9)
archive(gs_GAG, file = file.path(archive_path, "HVTN080-GAG.tar"))

HVTN080_population_stats <- list()
HVTN080_population_stats$negctrl <- getPopStats(gs_negctrl)
HVTN080_population_stats$ENV <- getPopStats(gs_ENV)
HVTN080_population_stats$GAG <- getPopStats(gs_GAG)

HVTN080_pData_gs_manual <- rbind(pData(gs_negctrl), pData(gs_ENV), pData(gs_GAG))
HVTN080_pData_gs_manual <- subset(HVTN080_pData_gs_manual, select = c(name, PTID, Stim, VISITNO))

save(HVTN080_population_stats, HVTN080_pData_gs_manual, file = "data/HVTN080-results.RData")

