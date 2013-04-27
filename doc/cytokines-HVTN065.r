#' # Cytokine Plots - HVTN065
#'
#' In this document for each HVTN 065 patient we plot densities and first
#' derivatives for the following three cytokines:
#'
#' 1. TNFa
#' 2. IFNg
#' 3. IL2
#'
#' For each patient, we have 2 visits: 2 (pre-vaccine) and 12
#' (post-vaccine). For each visit, we have 5 samples:
#'
#' * 1 GAG-stimulated
#' * 1 ENV-stimulated
#' * 2 negative controls
#' * 1 SEB control
#'
#' Also, for each patient sample, we consider the cytokines conditional on both
#' CD4 and CD8. Hence, we have 60 (3 x 2 x 5 x 2) cytokine plots per patient.
#'
#' The kernel density estimates are first centered by their mode and
#' then scaled with respect to the GAG stimulation.
#'
#' The plots below do not include SEB controls. These will be included soon
#' after OpenCyto is applied.
#'
#' Our goal is to determine an improved gate for the cytokines.

#+ setup, include=FALSE, cache=FALSE, echo=FALSE, warning=FALSE
opts_chunk$set(fig.align = 'default', dev = 'png', message = FALSE, warning = FALSE, error = FALSE,
               cache = FALSE, echo = FALSE, fig.path = 'figure/cytokines-HVTN065-',
               cache.path = 'cache/cytokines-HVTN065-', fig.width = 18, fig.height = 18,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()
library(flowIncubator)
library(openCyto)
library(MASS)

gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gs-cytokines")

#+ setup_data
 
# Constructs flowSets for the samples after applying the CD4 and CD8 gates
fs_CD4 <- getData(gs_HVTN065, which(getNodes(gs_HVTN065[[1]]) == "cd4"))
fs_CD8 <- getData(gs_HVTN065, which(getNodes(gs_HVTN065[[1]]) == "cd8"))

# Updates pData(...) to factors
pData_CD4 <- pData(fs_CD4)
pData_CD4$PTID <- factor(pData_CD4$PTID)
pData_CD4$VISITNO <- factor(pData_CD4$VISITNO)
pData_CD4$Stim <- factor(pData_CD4$Stim)
pData(fs_CD4) <- pData_CD4

pData_CD8 <- pData(fs_CD8)
pData_CD8$PTID <- factor(pData_CD8$PTID)
pData_CD8$VISITNO <- factor(pData_CD8$VISITNO)
pData_CD8$Stim <- factor(pData_CD8$Stim)
pData(fs_CD8) <- pData_CD8

#' ## CD4 Cytokines

#+ CD4

# For each patient ID, generates a plot of cytokines conditional on CD4.
ggplot_list <- lapply(levels(pData(fs_CD4)$PTID), function(current_PTID) {
  fcs_files <- as.character(subset(pData(fs_CD4), PTID == current_PTID)$name)

  # Create ggplot2-friendly summary of cytokine data conditional on CD4
  cytokine_data <- lapply(fcs_files, function(fcs_file) {
    TNFa <- as.vector(exprs(fs_CD4[[fcs_file]][, "Alexa 680-A"]))
    IFNg <- as.vector(exprs(fs_CD4[[fcs_file]][, "PE Cy7-A"]))
    IL2 <- as.vector(exprs(fs_CD4[[fcs_file]][, "PE Green laser-A"]))

    pData_fcs <- subset(pData(fs_CD4), name == fcs_file)
    cbind.data.frame(name = fcs_file, PTID = pData_fcs$PTID, Stim = pData_fcs$Stim,
                     VISITNO = pData_fcs$VISITNO, TNFa = TNFa, IFNg, IL2)
  })
  cytokine_data <- do.call(rbind, cytokine_data)

  # BUG: For some reason, 'melt' is not renaming the variable. We fix
  # this manually.
  cytokine_data <- reshape2:::melt(cytokine_data, variable.name = "Cytokine")
  colnames(cytokine_data) <- gsub("variable", "Cytokine", colnames(cytokine_data))

  # Scales the data with respect to the GAG sample
  # TODO: Per RG, the scaling should be respect to the negative peak of the SEB controls.
  #   After gating is finished, update the scaling to negative peak of SEB controls.
  stim_standardize <- "GAG-1-PTEG"
  cytokine_data <- ddply(cytokine_data, .(VISITNO, Cytokine), function(x) {
    # Calculates Huber estimates for the cytokine samples from reference
    # stimulation
    ref_huber <- huber(x[x$Stim == stim_standardize, ]$value)

    # Standardizes the cytokine samples within stimulation group with respect to
    # the reference stimulation group.

    # First, centers the values by the mode of the kernel density estimate for
    # the stimulation group. Then, scales by the Huber estimator of the standard
    # deviation
    x <- ddply(x, .(Stim), transform, value = center_mode(value))
    x <- ddply(x, .(Stim), transform, value = scale_huber(value, center = FALSE))
    x$value <- x$value * ref_huber$s
    x
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(deriv_smooth(x$value, n = 1024))
  })

  # Plot of cytokine densities for each stimulation group
  p1 <- ggplot(cytokine_data, aes(x = value, color = Stim, group = Stim))
  p1 <- p1 + geom_density() + theme_bw()
  p1 <- p1 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p1 <- p1 + ggtitle(paste("Scaled Cytokine Densities"))

  # Plot of derivatives of smoothed densities for each stimulation group
  p2 <- ggplot(first_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p2 <- p2 + geom_line() + theme_bw()
  p2 <- p2 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p2 <- p2 + ggtitle(paste("First Derivatives")) + ylab("dy/dx")

  # Creates a single plot containing cytokine densities and derivatives.
  # The plot shares the legend.
  # For more details, see Stack Overflow post: http://bit.ly/15J5nYT
  p2_legend <- g_legend(p2)
  grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                           p2 + theme(legend.position="none"),
                           main = paste("CD4 Cytokines -- Patient:", current_PTID),
                           nrow = 1),
               p2_legend, nrow = 2, heights = c(10, 2))
})

lapply(ggplot_list, plot)

#' ## CD8 Cytokines

#+ CD8

# For each patient ID, generates a plot of cytokines conditional on CD8.
ggplot_list <- lapply(levels(pData(fs_CD8)$PTID), function(current_PTID) {
  fcs_files <- as.character(subset(pData(fs_CD8), PTID == current_PTID)$name)

  # Create ggplot2-friendly summary of cytokine data conditional on CD8
  cytokine_data <- lapply(fcs_files, function(fcs_file) {
    TNFa <- as.vector(exprs(fs_CD8[[fcs_file]][, "Alexa 680-A"]))
    IFNg <- as.vector(exprs(fs_CD8[[fcs_file]][, "PE Cy7-A"]))
    IL2 <- as.vector(exprs(fs_CD8[[fcs_file]][, "PE Green laser-A"]))

    pData_fcs <- subset(pData(fs_CD8), name == fcs_file)
    cbind.data.frame(name = fcs_file, PTID = pData_fcs$PTID, Stim = pData_fcs$Stim,
                     VISITNO = pData_fcs$VISITNO, TNFa = TNFa, IFNg, IL2)
  })
  cytokine_data <- do.call(rbind, cytokine_data)

  # BUG: For some reason, 'melt' is not renaming the variable. We fix
  # this manually.
  cytokine_data <- reshape2:::melt(cytokine_data, variable.name = "Cytokine")
  colnames(cytokine_data) <- gsub("variable", "Cytokine", colnames(cytokine_data))

  # Scales the data with respect to the GAG sample
  # TODO: Per RG, the scaling should be respect to the negative peak of the SEB controls.
  #   After gating is finished, update the scaling to negative peak of SEB controls.
  stim_standardize <- "GAG-1-PTEG"
  cytokine_data <- ddply(cytokine_data, .(VISITNO, Cytokine), function(x) {
    # Calculates Huber estimates for the cytokine samples from reference
    # stimulation
    ref_huber <- huber(x[x$Stim == stim_standardize, ]$value)

    # Standardizes the cytokine samples within stimulation group with respect to
    # the reference stimulation group.
    x <- ddply(x, .(Stim), transform, value = scale_huber(value))
    x$value <- x$value * ref_huber$s
    x
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(deriv_smooth(x$value, n = 1024))
  })

  # Plot of cytokine densities for each stimulation group
  p1 <- ggplot(cytokine_data, aes(x = value, color = Stim, group = Stim))
  p1 <- p1 + geom_density() + theme_bw()
  p1 <- p1 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p1 <- p1 + ggtitle(paste("Scaled Cytokine Densities"))

  # Plot of derivatives of smoothed densities for each stimulation group
  p2 <- ggplot(first_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p2 <- p2 + geom_line() + theme_bw()
  p2 <- p2 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p2 <- p2 + ggtitle(paste("First Derivatives")) + ylab("dy/dx")

  # Creates a single plot containing cytokine densities and derivatives.
  # The plot shares the legend.
  # For more details, see Stack Overflow post: http://bit.ly/15J5nYT
  p2_legend <- g_legend(p2)
  grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                           p2 + theme(legend.position="none"),
                           main = paste("CD8 Cytokines -- Patient:", current_PTID),
                           nrow = 1),
               p2_legend, nrow = 2, heights = c(10, 2))
})
