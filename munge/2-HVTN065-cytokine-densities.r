# For the HVTN065, we create density summaries of the cytokine cells conditional
# on both CD4+ and CD8+.
# These density summaries are used in the Cytokine Densities Shiny app.

library(ProjectTemplate)
load.project()

num_cores <- 12

gs_HVTN065 <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating-results")

# Updates pData(...) to factors
pData_HVTN065 <- pData(gs_HVTN065)
pData_HVTN065$PTID <- factor(pData_HVTN065$PTID)
pData_HVTN065$VISITNO <- factor(pData_HVTN065$VISITNO)
pData_HVTN065$Stim <- factor(pData_HVTN065$Stim)
pData(gs_HVTN065) <- pData_HVTN065

# Preload the cytokine summaries for all patients and T-cells
# Traverse through each (PTID, VISITNO, Tcell, cytokine) tuple
combos <- expand.grid(PTID = levels(pData_HVTN065$PTID),
                      VISITNO = levels(pData_HVTN065$VISITNO),
                      tcells = c("cd4", "cd8"))

cytokine_densities <- mclapply(seq_along(combos$PTID), function(combo_i) {
  current_PTID <- as.character(combos$PTID[combo_i])
  current_VISITNO <- as.character(combos$VISITNO[combo_i])
  current_tcells <- as.character(combos$tcells[combo_i])

  # Reads in TNFa data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_TNFa <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "TNFa")
  x_TNFa <- standardize_cytokines(x_TNFa)

  density_TNFa <- tapply(x_TNFa$value, x_TNFa$Stim, density, n = 1024)
  density_TNFa <- lapply(seq_along(density_TNFa), function(i) {
    cbind(Stim = names(density_TNFa)[i], x = density_TNFa[[i]]$x,
          y = density_TNFa[[i]]$y)
  })
  density_TNFa <- do.call(rbind, density_TNFa)
  density_TNFa <- cbind(Cytokine = "TNFa", density_TNFa)

  # Reads in IFNg data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_IFNg <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "IFNg")
  x_IFNg <- standardize_cytokines(x_IFNg)
  density_IFNg <- tapply(x_IFNg$value, x_IFNg$Stim, density, n = 1024)
  density_IFNg <- lapply(seq_along(density_IFNg), function(i) {
    cbind(Stim = names(density_IFNg)[i], x = density_IFNg[[i]]$x,
          y = density_IFNg[[i]]$y)
  })
  density_IFNg <- do.call(rbind, density_IFNg)
  density_IFNg <- cbind(Cytokine = "IFNg", density_IFNg)

  # Reads in IL2 data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_IL2 <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "IL2")
  x_IL2 <- standardize_cytokines(x_IL2)
  density_IL2 <- tapply(x_IL2$value, x_IL2$Stim, density, n = 1024)
  density_IL2 <- lapply(seq_along(density_IL2), function(i) {
    cbind(Stim = names(density_IL2)[i], x = density_IL2[[i]]$x,
          y = density_IL2[[i]]$y)
  })
  density_IL2 <- do.call(rbind, density_IL2)
  density_IL2 <- cbind(Cytokine = "IL2", density_IL2)

  density_cytokines <- rbind(density_TNFa, density_IFNg, density_IL2)
  density_cytokines <- data.frame(density_cytokines, stringsAsFactors = FALSE)
  density_cytokines$x <- as.numeric(density_cytokines$x)
  density_cytokines$y <- as.numeric(density_cytokines$y)

  density_cytokines$PTID <- current_PTID
  density_cytokines$VISITNO <- current_VISITNO
  density_cytokines$tcells <- current_tcells

  density_cytokines
}, mc.cores = num_cores)

cytokine_densities <- do.call(rbind, cytokine_densities)

# Stopped the run here.

cytokine_densities <- data.frame(cytokine_densities, stringsAsFactors = FALSE)
cytokine_densities$x <- as.numeric(cytokine_densities$x)
cytokine_densities$y <- as.numeric(cytokine_densities$y)

# TODO: Look at summary here to see if any NAs again from coercion

# TODO: Currently, there are 27 instances that occur when creating the densities.
# These are coerced to NAs when we coerce 'x' and 'y' to numeric. For now, we omit
# these cases.
cytokine_densities <- subset(cytokine_densities, !is.na(x) & !is.na(y))

cytokine_densities <- data.table(cytokine_densities)

save(cytokine_densities, file = "cache/cytokine_densities.RData")




# TEMP: Remove me

# The following patient had the scale issue that Greg, Raphael, and I
# saw in the Shiny app. This patient is one of several.
cytokines_PTID <- read_cytokines(gs_HVTN065, "125290127", "2", "cd4", "TNFa")
cytokines_PTID <- standardize_cytokines(cytokines_PTID)

  density_TNFa <- tapply(x_TNFa$value, x_TNFa$Stim, density, n = 1024)
  density_TNFa <- lapply(seq_along(density_TNFa), function(i) {
    cbind(Stim = names(density_TNFa)[i], x = density_TNFa[[i]]$x,
          y = density_TNFa[[i]]$y)
  })
  density_TNFa <- do.call(rbind, density_TNFa)
  density_TNFa <- cbind(Cytokine = "TNFa", density_TNFa)
