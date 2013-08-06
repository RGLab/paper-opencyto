library(ProjectTemplate)
load.project()

# Loads Gating Set
gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-set/")
# gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-results/")

# Loads the GatingTemplate from the CSV file.
gating_template <- gatingTemplate("gt-HVTN065.csv")

# Now we apply the automated pipeline to each gating set
set.seed(42)
start_time <- Sys.time()
try(gating(gating_template, gs_HVTN065, prior_group = 'Stim', mc.cores = 10,
       parallel_type = "multicore"))
finish_time <- Sys.time()
message("Time Elapsed:")
print(finish_time - start_time)

# In our summary scripts we require that the 'Stim' string be unique within a
# patient-visit pairing, e.g., plotGate with conditional panels. We ensure that
# the negative control strings are unique within a patient-visit pairing.
pData_HVTN065 <- pData(gs_HVTN065)
pData_HVTN065$Stim <- as.character(pData_HVTN065$Stim)
pData_negctrl <- subset(pData_HVTN065, Stim == "negctrl")
pData_negctrl <- ddply(pData_negctrl, .(PTID, VISITNO), transform,
                       Stim = paste0(Stim, seq_along(Stim)))
pData_HVTN065[pData_HVTN065$Stim == "negctrl", ] <- pData_negctrl
pData(gs_HVTN065) <- pData_HVTN065

# Archives the results
save_gs(gs_HVTN065, path = "/loc/no-backup/ramey/HVTN/065/gating-set", overwrite = TRUE)

# Extracts the population statistics
popstats_HVTN065 <- getPopStats(gs_HVTN065)
counts_HVTN065 <- getPopStats(gs_HVTN065, statistic = "count")

# Extracts MFI for each sample
cytokine_markers <- data.frame(channel = c("Alexa 680-A", "PE Cy7-A", "PE Green laser-A"),
                               marker = c("TNFa", "IFNg", "IL2"), stringsAsFactors = FALSE)

fs_cd4 <- getData(gs_HVTN065, "cd4")[, cytokine_markers$channel]
fs_cd8 <- getData(gs_HVTN065, "cd8")[, cytokine_markers$channel]


mfi_cd4 <- sapply(sampleNames(fs_cd4), function(fcs_file) {
  each_col(fs_cd4[[fcs_file]], median)
})
mfi_cd4 <- data.table(melt(mfi_cd4))
setnames(mfi_cd4, colnames(mfi_cd4), c("Cytokine", "Sample", "MFI"))

mfi_cd8 <- sapply(sampleNames(fs_cd8), function(fcs_file) {
  each_col(fs_cd8[[fcs_file]], median)
})
mfi_cd8 <- data.table(melt(mfi_cd8))
setnames(mfi_cd8, colnames(mfi_cd8), c("Cytokine", "Sample", "MFI"))


save(counts_HVTN065, popstats_HVTN065, pData_HVTN065, mfi_cd4, mfi_cd8,
     file = "data/HVTN065-results.RData")

