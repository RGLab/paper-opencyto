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

cytokine_summary <- mclapply(seq_along(combos$PTID), function(combo_i) {
  current_PTID <- as.character(combos$PTID[combo_i])
  current_VISITNO <- as.character(combos$VISITNO[combo_i])
  current_tcells <- as.character(combos$tcells[combo_i])

  # Reads in TNFa data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_TNFa <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "TNFa")
  x_TNFa <- standardize_cytokines(x_TNFa)
  density_TNFa <- cbind(Cytokine = "TNFa", density_cytokines(x_TNFa))
  deriv_TNFa <- cbind(Cytokine = "TNFa", deriv_cytokines(x_TNFa, deriv = 1, adjust = 3))
  deriv_collapse_TNFa <- cbind(Cytokine = "TNFa", deriv_cytokines(x_TNFa, deriv = 1, adjust = 3, collapse = TRUE))
  second_deriv_TNFa <- cbind(Cytokine = "TNFa", deriv_cytokines(x_TNFa, deriv = 2, adjust = 3))

  # Reads in IFNg data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_IFNg <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "IFNg")
  x_IFNg <- standardize_cytokines(x_IFNg)
  density_IFNg <- cbind(Cytokine = "IFNg", density_cytokines(x_IFNg))
  deriv_IFNg <- cbind(Cytokine = "IFNg", deriv_cytokines(x_IFNg, deriv = 1, adjust = 3))
  deriv_collapse_IFNg <- cbind(Cytokine = "IFNg", deriv_cytokines(x_IFNg, deriv = 1, adjust = 3, collapse = TRUE))
  second_deriv_IFNg <- cbind(Cytokine = "IFNg", deriv_cytokines(x_IFNg, deriv = 2, adjust = 3))

  # Reads in IL2 data, melts the data, standardizes them, and then computes the
  # kernel density estimates
  x_IL2 <- read_cytokines(gs_HVTN065, PTID = current_PTID, VISITNO = current_VISITNO,
                      tcells = current_tcells, cytokine = "IL2")
  x_IL2 <- standardize_cytokines(x_IL2)
  density_IL2 <- cbind(Cytokine = "IL2", density_cytokines(x_IL2))
  deriv_IL2 <- cbind(Cytokine = "IL2", deriv_cytokines(x_IL2, deriv = 1, adjust = 3))
  deriv_collapse_IL2 <- cbind(Cytokine = "IL2", deriv_cytokines(x_IL2, deriv = 1, adjust = 3, collapse = TRUE))
  second_deriv_IL2 <- cbind(Cytokine = "IL2", deriv_cytokines(x_IL2, deriv = 2, adjust = 3))

  density_cytokines <- rbind(density_TNFa, density_IFNg, density_IL2)
  density_cytokines <- data.frame(density_cytokines, stringsAsFactors = FALSE)
  density_cytokines$x <- as.numeric(density_cytokines$x)
  density_cytokines$y <- as.numeric(density_cytokines$y)
  density_cytokines$PTID <- current_PTID
  density_cytokines$VISITNO <- current_VISITNO
  density_cytokines$tcells <- current_tcells

  deriv_cytokines <- rbind(deriv_TNFa, deriv_IFNg, deriv_IL2)
  deriv_cytokines <- data.frame(deriv_cytokines, stringsAsFactors = FALSE)
  deriv_cytokines$x <- as.numeric(deriv_cytokines$x)
  deriv_cytokines$y <- as.numeric(deriv_cytokines$y)
  deriv_cytokines$PTID <- current_PTID
  deriv_cytokines$VISITNO <- current_VISITNO
  deriv_cytokines$tcells <- current_tcells

  deriv_collapse_cytokines <- rbind(deriv_collapse_TNFa, deriv_collapse_IFNg, deriv_collapse_IL2)
  deriv_collapse_cytokines <- data.frame(deriv_collapse_cytokines, stringsAsFactors = FALSE)
  deriv_collapse_cytokines$x <- as.numeric(deriv_collapse_cytokines$x)
  deriv_collapse_cytokines$y <- as.numeric(deriv_collapse_cytokines$y)
  deriv_collapse_cytokines$PTID <- current_PTID
  deriv_collapse_cytokines$VISITNO <- current_VISITNO
  deriv_collapse_cytokines$tcells <- current_tcells

  second_deriv_cytokines <- rbind(second_deriv_TNFa, second_deriv_IFNg, second_deriv_IL2)
  second_deriv_cytokines <- data.frame(second_deriv_cytokines, stringsAsFactors = FALSE)
  second_deriv_cytokines$x <- as.numeric(second_deriv_cytokines$x)
  second_deriv_cytokines$y <- as.numeric(second_deriv_cytokines$y)
  second_deriv_cytokines$PTID <- current_PTID
  second_deriv_cytokines$VISITNO <- current_VISITNO
  second_deriv_cytokines$tcells <- current_tcells

  list(density = density_cytokines, deriv = deriv_cytokines,
       deriv_collapse = deriv_collapse_cytokines, second_deriv = second_deriv_cytokines)
}, mc.cores = num_cores)

cytokine_densities <- do.call(rbind, lapply(cytokine_summary, function(x) x$density))
cytokine_derivs <- do.call(rbind, lapply(cytokine_summary, function(x) x$deriv))
cytokine_derivs_collapse <- do.call(rbind, lapply(cytokine_summary, function(x) x$deriv_collapse))
cytokine_second_derivs <- do.call(rbind, lapply(cytokine_summary, function(x) x$second_deriv))

# To speed up the lookup of the cytokine densities and derivatives, we convert
# them to data.table objects
cytokine_densities <- data.table(cytokine_densities)
cytokine_derivs <- data.table(cytokine_derivs)
cytokine_derivs_collapse <- data.table(cytokine_derivs_collapse)
cytokine_second_derivs <- data.table(cytokine_second_derivs)


save(cytokine_densities, cytokine_derivs, cytokine_derivs_collapse, cytokine_second_derivs,
     file = "cache/cytokine_densities.RData")


