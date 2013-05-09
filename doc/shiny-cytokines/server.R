shinyServer(function(input, output) {
  PTID_fcs <- reactive({
    subset(pData_HVTN065, PTID == input$PTID)$name
  })

  cytokine_summary <- reactive({
    subset(cytokine_densities, PTID == input$PTID & tcells == input$tcells)
  })
  
  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- reactive({
    subset(cytokine_derivs, PTID == input$PTID & tcells == input$tcells)
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokine, and Stimulation
  second_derivs <- reactive({
    subset(cytokine_second_derivs, PTID == input$PTID & tcells == input$tcells)
  })

  # Determines the cutpoint based on the derivatives of the kernel density
  # estimates from the cells collapsed across stimulation groups for pair of
  # VISITNO and Cytokine
  cyto_cutpoints <- reactive({
    # In the UI, we display the tolerance values on the log-10 scale. Here, we
    # have to scale the value manually.
    tol <- 10^input$cutpoint_tol

    derivs_collapse <- subset(cytokine_derivs_collapse,
                              PTID == input$PTID & tcells == input$tcells)
    cutpoints <- ddply(derivs_collapse, .(Cytokine, VISITNO), function(foo) {
      lowest_valley <- with(foo, x[which.min(y)])
      cutpoint <- with(foo, x[which(x > lowest_valley & abs(y) < tol)][1])
      cutpoint
    })
    colnames(cutpoints) <- c("Cytokine", "VISITNO", "cutpoint")
    cutpoints
  })

  # Plot of cytokine densities for each stimulation group
  output$densitiesPlot <- renderPlot({
    p1 <- ggplot(cytokine_summary(), aes(x = x, y = y, color = Stim, group = Stim))
    p1 <- p1 + geom_path() + theme_bw()
    p1 <- p1 + facet_grid(VISITNO ~ Cytokine)
    p1 <- p1 + xlim(input$x_range_density[1], input$x_range_density[2])
    p1 <- p1 + ylim(input$y_range_density[1], input$y_range_density[2])
    p1 <- p1 + ggtitle("Scaled Cytokine Densities")

    if (input$display_cutpoint) {
      p1 <- p1 + geom_vline(linetype = "dashed", aes(xintercept = cutpoint), data = cyto_cutpoints())
    }
    
    plot(p1)
  })

  # Plot of derivatives of smoothed densities for each stimulation group
  output$firstDerivPlot <- renderPlot({
    p2 <- ggplot(first_derivs(), aes(x = x, y = y, color = Stim, group = Stim))
    p2 <- p2 + geom_line() + theme_bw()
    p2 <- p2 + facet_grid(VISITNO ~ Cytokine)
    p2 <- p2 + xlim(input$x_range_first_deriv[1], input$x_range_first_deriv[2])
    p2 <- p2 + ylim(input$y_range_first_deriv[1], input$y_range_first_deriv[2])
    p2 <- p2 + ylab("dy/dx")

    if (input$display_cutpoint) {
      p2 <- p2 + geom_vline(linetype = "dashed", aes(xintercept = cutpoint), data = cyto_cutpoints())
    }
    
    plot(p2)
  })

  # Plot of second derivatives of smoothed densities for each stimulation group
  output$secondDerivPlot <- renderPlot({
    p3 <- ggplot(second_derivs(), aes(x = x, y = y, color = Stim, group = Stim))
    p3 <- p3 + geom_line() + theme_bw()
    p3 <- p3 + facet_grid(VISITNO ~ Cytokine, scales = "free_x")
    p3 <- p3 + xlim(input$x_range_second_deriv[1], input$x_range_second_deriv[2])
    p3 <- p3 + ylim(input$y_range_second_deriv[1], input$y_range_second_deriv[2])
    p3 <- p3 + ylab("d^2y/dx^2")

    if (input$display_cutpoint) {
      p3 <- p3 + geom_vline(linetype = "dashed", aes(xintercept = cutpoint), data = cyto_cutpoints())
    }
    
    plot(p3)
  })

  output$TNFaGatesPlot <- renderPlot({
    if (input$TNFa_gates) {
      fcs_files <- PTID_fcs()
      tcells <- input$tcells
      print(plotGate(gs_HVTN065[fcs_files], paste(tcells, "TNFa", sep = "/"),
                     lattice = TRUE, xbin = 128, cond = "factor(Stim):factor(VISITNO)"))
    } else {
      return(NULL)
    }
  })

  output$IFNgGatesPlot <- renderPlot({
    if (input$IFNg_gates) {
      fcs_files <- PTID_fcs()
      tcells <- input$tcells
      print(plotGate(gs_HVTN065[fcs_files], paste(tcells, "IFNg", sep = "/"),
                     lattice = TRUE, xbin = 128, cond = "factor(Stim):factor(VISITNO)"))
    } else {
      return(NULL)
    }
  })

  output$IL2GatesPlot <- renderPlot({
    if (input$IL2_gates) {
      fcs_files <- PTID_fcs()
      tcells <- input$tcells
      print(plotGate(gs_HVTN065[fcs_files], paste(tcells, "IL2", sep = "/"),
                     lattice = TRUE, xbin = 128, cond = "factor(Stim):factor(VISITNO)"))
    } else {
      return(NULL)
    }
  })

})
