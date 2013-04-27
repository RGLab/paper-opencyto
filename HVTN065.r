start_time <- Sys.time()

library(ProjectTemplate)
load.project()

# Loads Gating Set
gs_HVTN065 <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating_set")

# Loads the GatingTemplate from the CSV file.
gating_template <- gatingTemplate("HVTN065-GatingTemplate.csv", "HVTN065")

# Now we apply the automated pipeline to each gating set
set.seed(42)
gating(gating_template, gs_HVTN065, prior_group = 'Stim', num_cores = 12,
       parallel_type = "multicore")

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

# Archives the results
save_gs(gs_HVTN065,
        path = "/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating-results")

# Extracts the population statistics
popstats_HVTN065 <- getPopStats(gs_HVTN065)

save(popstats_HVTN065, pData_HVTN065, file = "data/HVTN065-results.RData")

finish_time <- Sys.time()
message("Time Elapsed:")
print(finish_time - start_time)

