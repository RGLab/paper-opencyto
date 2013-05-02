shinyServer(function(input, output) {
  cytokine_summary <- reactive({
    subset(cytokine_densities, PTID == input$PTID & tcells == input$tcells)
  })
  
  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokin, and Stimulation
  first_derivs <- reactive({
    ddply(cytokine_summary(), .(VISITNO, Cytokine, Stim), function(x) {
      as.data.frame(deriv_smooth(x))
    })
  })

  # Calculates a discrete derivative based on the kernel density estimate for
  # each combination of VISITNO, Cytokine, and Stimulation
  second_derivs <- reactive({
    ddply(cytokine_summary(), .(VISITNO, Cytokine, Stim), function(x) {
      as.data.frame(second_deriv_smooth(x, smooth = TRUE))
    })
  })

  # Plot of cytokine densities for each stimulation group
  output$densitiesPlot <- renderPlot({
    p1 <- ggplot(cytokine_summary(), aes(x = x, y = y, color = Stim, group = Stim))
    p1 <- p1 + geom_path() + theme_bw()
    p1 <- p1 + facet_grid(VISITNO ~ Cytokine, scales = "free_x") + ylim(0, input$max_density)
    p1 <- p1 + ggtitle("Scaled Cytokine Densities")
    plot(p1)
  })

  # Plot of derivatives of smoothed densities for each stimulation group
  output$firstDerivPlot <- renderPlot({
    p2 <- ggplot(first_derivs(), aes(x = x, y = y, color = Stim, group = Stim))
    p2 <- p2 + geom_line() + theme_bw()
    p2 <- p2 + facet_grid(VISITNO ~ Cytokine, scales = "free_x")
    p2 <- p2 + ylim(input$range_first_deriv[1], input$range_first_deriv[2])
    p2 <- p2 + ylab("dy/dx")
    plot(p2)
  })

  # Plot of second derivatives of smoothed densities for each stimulation group
  output$secondDerivPlot <- renderPlot({
    p3 <- ggplot(second_derivs(), aes(x = x, y = y, color = Stim, group = Stim))
    p3 <- p3 + geom_line() + theme_bw()
    p3 <- p3 + facet_grid(VISITNO ~ Cytokine, scales = "free_x")
    p3 <- p3 + ylim(input$range_second_deriv[1], input$range_second_deriv[2])
    p3 <- p3 + ylab("d^2y/dx^2")
    plot(p3)
  })

#  output$gatesPlot <- renderPlot({
#    print(plotGate(gs_HVTN065[PTID_fcs()], lattice = TRUE, xbin = 128, cond = "factor(Stim):factor(VISITNO)"))
#  })

})
