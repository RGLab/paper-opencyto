setwd("~/rglab/papers/paper-opencyto")
library(ProjectTemplate)
load.project()

gs_HVTN065 <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating-results")

# Updates pData(...) to factors
pData_HVTN065 <- pData(gs_HVTN065)
pData_HVTN065$PTID <- factor(pData_HVTN065$PTID)
pData_HVTN065$VISITNO <- factor(pData_HVTN065$VISITNO)
pData_HVTN065$Stim <- factor(pData_HVTN065$Stim)
pData(gs_HVTN065) <- pData_HVTN065

# Sets constants for the Shiny UI
levels_PTID <- levels(pData_HVTN065$PTID)

