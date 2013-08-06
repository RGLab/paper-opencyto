setwd("~/rglab/papers/paper-opencyto")
library(ProjectTemplate)
load.project()

# gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-results")
gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-set")


# Updates pData(...) to factors
pData_HVTN065 <- pData(gs_HVTN065)
pData_HVTN065$PTID <- factor(pData_HVTN065$PTID)
pData_HVTN065$VISITNO <- factor(pData_HVTN065$VISITNO)
pData_HVTN065$Stim <- factor(pData_HVTN065$Stim)
pData(gs_HVTN065) <- pData_HVTN065

# Adds treatment info to pData.
treatment_info <- data.frame(PTID = gsub("-", "", as.character(treatment.HVTN065$Ptid)), 
  Treatment = "Treatment", stringsAsFactors = FALSE)
treatment_info$Treatment <- replace(treatment_info$Treatment,
                                    grep("^Placebo", treatment.HVTN065$rx),
                                    "Placebo")

pData_HVTN065 <- plyr:::join(pData_HVTN065, treatment_info)


# Sets constants for the Shiny UI
levels_PTID <- levels(pData_HVTN065$PTID)

# Extract cytokine tolerance values from gating nodes.
cytokine_tolerances <- grep("_tol", getNodes(gs_HVTN065[[1]]), value = TRUE)
cytokine_tolerances <- unique(sapply(strsplit(cytokine_tolerances, "_tol"), tail, n = 1))
cytokine_tolerances <- sort(as.numeric(cytokine_tolerances), decreasing = TRUE)
