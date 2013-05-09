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
               cache.path = 'cache/cytokines-HVTN065-', fig.width = 24, fig.height = 15,
               results = 'hide')

#+ load_data  
setwd("..")
library(ProjectTemplate)
load.project()
library(flowIncubator)
library(openCyto)
library(MASS)

gs_HVTN065 <- load_gs("/shared/silo_researcher/Gottardo_R/ramey_working/HVTN/065/gating-results")

#+ setup_data
 
# Constructs flowSets for the samples after applying the CD4 and CD8 gates
fs_CD4 <- getData(gs_HVTN065, "cd4")
fs_CD8 <- getData(gs_HVTN065, "cd8")

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

  # Scales the data with respect to the negative component of the SEB control
  # sample.
  cytokine_data <- ddply(cytokine_data, .(VISITNO, Cytokine), function(x) {
    x_standardize <- x[x$Stim == "sebctrl", ]$value
    sebctrl_peaks <- openCyto:::find_peaks(x_standardize, order = TRUE, adjust = 3)
    neg_peak <- sebctrl_peaks[1]
    pos_peak <- sebctrl_peaks[2]

    # If two peaks are present, select cutpoint as mindensity
    # Otherwise, select based on smoothed second derivative
    if (!is.na(pos_peak)) {
      valleys <- openCyto:::find_valleys(x_standardize, adjust = 2)
      cutpoint <- valleys[findInterval(valleys, c(neg_peak, pos_peak)) == 1][1]
    } else {
      second_deriv_standardize <- second_deriv_smooth(x_standardize, n = 4096, adjust = 2)
      second_deriv_standardize <- do.call(cbind.data.frame, second_deriv_standardize)

      # Applies LOESS with minor smoothing to second derivative curve
      loess_out <- loess(y ~ x, data = second_deriv_standardize, span = 0.2)

      # Select cutpoint as first valley after peak
      # Because the LOESS prediction is likely unequal to the negative peak found
      # using openCyto, we find the first valley of the second derivative greater
      # than the valley nearest to the negative_peak
      x_sorted <- sort(x_standardize)
      predict_loess <- predict(loess_out, x_sorted)
      discrete_second_deriv <- diff(sign(diff(predict_loess)))
      which_minima <- which(discrete_second_deriv == 2) + 1
      valleys <- x_sorted[which_minima]
      cutpoint <- valleys[which.min(abs(valleys - neg_peak)) + 1]
    }

    # Calculates standard deviation estimate for the negative-component
    # observations that are greater than the negative peak.
    x_pos <- x_standardize[findInterval(x_standardize, c(neg_peak, cutpoint)) == 1]
    sd_neg <- sqrt(mean((x_pos - neg_peak)^2))

    # Standardizes the cytokine samples within stimulation group with respect to
    # the reference stimulation group.
    # First, centers the values by the mode of the kernel density estimate for
    # the stimulation group.
    x <- ddply(x, .(Stim), transform, value = center_mode(value))

    # For the reference stimulation group (i.e., SEB contrtols), we scale by the
    # standard deviation of its negative component.
    x_ref <- subset(x, Stim == "sebctrl")
    x_ref$value <- x_ref$value / sd_neg

    # Scales by the Huber estimator of the standard deviation and rescales with respect to the
    # reference stimulation sample to put all stimulations on the same scale.
    x_nonref <- ddply(subset(x, Stim != "sebctrl"), .(Stim), transform,
                      value = scale_huber(value, center = FALSE))
    x_nonref$value <- x_nonref$value * sd_neg

    rbind(x_ref, x_nonref)
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(deriv_smooth(x$value, n = 1024, adjust = 2))
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  second_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(second_deriv_smooth(x$value, n = 1024, adjust = 2))
  })

  # Plot of cytokine densities for each stimulation group
  p1 <- ggplot(cytokine_data, aes(x = value, color = Stim, group = Stim))
  p1 <- p1 + geom_density() + theme_bw()
  p1 <- p1 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p1 <- p1 + ggtitle("Scaled Cytokine Densities")

  # Plot of derivatives of smoothed densities for each stimulation group
  p2 <- ggplot(first_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p2 <- p2 + geom_line() + theme_bw()
  p2 <- p2 + facet_grid(VISITNO ~ Cytokine, scales = "free_x") + ylim(-1, 1)
  p2 <- p2 + ggtitle("First Derivatives") + ylab("dy/dx")

  # Plot of second derivatives of smoothed densities for each stimulation group
  p3 <- ggplot(second_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p3 <- p3 + geom_line() + theme_bw()
  p3 <- p3 + facet_grid(VISITNO ~ Cytokine, scales = "free_x") + ylim(-1, 1)
  p3 <- p3 + ggtitle("Second Derivatives") + ylab("d^2y/dx^2")

  # Creates a single plot containing cytokine densities and derivatives.
  # The plot shares the legend.
  # For more details, see Stack Overflow post: http://bit.ly/15J5nYT
  p2_legend <- g_legend(p2)
  grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                           p2 + theme(legend.position = "none"),
                           p3 + theme(legend.position = "none"),
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

  # Scales the data with respect to the negative component of the SEB control
  # sample.
  cytokine_data <- ddply(cytokine_data, .(VISITNO, Cytokine), function(x) {
    x_standardize <- x[x$Stim == "sebctrl", ]$value
    sebctrl_peaks <- openCyto:::find_peaks(x_standardize, order = TRUE, adjust = 3)
    neg_peak <- sebctrl_peaks[1]
    pos_peak <- sebctrl_peaks[2]

    # If two peaks are present, select cutpoint as mindensity
    # Otherwise, select based on smoothed second derivative
    if (!is.na(pos_peak)) {
      valleys <- openCyto:::find_valleys(x_standardize, adjust = 2)
      cutpoint <- valleys[findInterval(valleys, c(neg_peak, pos_peak)) == 1][1]
    } else {
      second_deriv_standardize <- second_deriv_smooth(x_standardize, n = 4096, adjust = 2)
      second_deriv_standardize <- do.call(cbind.data.frame, second_deriv_standardize)

      # Applies LOESS with minor smoothing to second derivative curve
      loess_out <- loess(y ~ x, data = second_deriv_standardize, span = 0.2)

      # Select cutpoint as first valley after peak
      # Because the LOESS prediction is likely unequal to the negative peak found
      # using openCyto, we find the first valley of the second derivative greater
      # than the valley nearest to the negative_peak
      x_sorted <- sort(x_standardize)
      predict_loess <- predict(loess_out, x_sorted)
      discrete_second_deriv <- diff(sign(diff(predict_loess)))
      which_minima <- which(discrete_second_deriv == 2) + 1
      valleys <- x_sorted[which_minima]
      cutpoint <- valleys[which.min(abs(valleys - neg_peak)) + 1]
    }

    # Calculates standard deviation estimate for the negative-component
    # observations that are greater than the negative peak.
    x_pos <- x_standardize[findInterval(x_standardize, c(neg_peak, cutpoint)) == 1]
    sd_neg <- sqrt(mean((x_pos - neg_peak)^2))

    # Standardizes the cytokine samples within stimulation group with respect to
    # the reference stimulation group.
    # First, centers the values by the mode of the kernel density estimate for
    # the stimulation group.
    x <- ddply(x, .(Stim), transform, value = center_mode(value))

    # For the reference stimulation group (i.e., SEB contrtols), we scale by the
    # standard deviation of its negative component.
    x_ref <- subset(x, Stim == "sebctrl")
    x_ref$value <- x_ref$value / sd_neg

    # Scales by the Huber estimator of the standard deviation and rescales with respect to the
    # reference stimulation sample to put all stimulations on the same scale.
    x_nonref <- ddply(subset(x, Stim != "sebctrl"), .(Stim), transform,
                      value = scale_huber(value, center = FALSE))
    x_nonref$value <- x_nonref$value * sd_neg

    rbind(x_ref, x_nonref)
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(deriv_smooth(x$value, n = 1024, adjust = 2))
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  second_derivs <- ddply(cytokine_data, .(VISITNO, Cytokine, Stim), function(x) {
    as.data.frame(second_deriv_smooth(x$value, n = 1024, adjust = 2))
  })

  # Plot of cytokine densities for each stimulation group
  p1 <- ggplot(cytokine_data, aes(x = value, color = Stim, group = Stim))
  p1 <- p1 + geom_density() + theme_bw()
  p1 <- p1 + facet_grid(VISITNO ~ Cytokine, scales = "free")
  p1 <- p1 + ggtitle("Scaled Cytokine Densities")

  # Plot of derivatives of smoothed densities for each stimulation group
  p2 <- ggplot(first_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p2 <- p2 + geom_line() + theme_bw()
  p2 <- p2 + facet_grid(VISITNO ~ Cytokine, scales = "free_x") + ylim(-1, 1)
  p2 <- p2 + ggtitle("First Derivatives") + ylab("dy/dx")

  # Plot of second derivatives of smoothed densities for each stimulation group
  p3 <- ggplot(second_derivs, aes(x = x, y = y, color = Stim, group = Stim))
  p3 <- p3 + geom_line() + theme_bw()
  p3 <- p3 + facet_grid(VISITNO ~ Cytokine, scales = "free_x") + ylim(-1, 1)
  p3 <- p3 + ggtitle("Second Derivatives") + ylab("d^2y/dx^2")

  # Creates a single plot containing cytokine densities and derivatives.
  # The plot shares the legend.
  # For more details, see Stack Overflow post: http://bit.ly/15J5nYT
  p2_legend <- g_legend(p2)
  grid.arrange(arrangeGrob(p1 + theme(legend.position = "none"),
                           p2 + theme(legend.position = "none"),
                           p3 + theme(legend.position = "none"),
                           main = paste("CD8 Cytokines -- Patient:", current_PTID),
                           nrow = 1),
               p2_legend, nrow = 2, heights = c(10, 2))
})

lapply(ggplot_list, plot)
