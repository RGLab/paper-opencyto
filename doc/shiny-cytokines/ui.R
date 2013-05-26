shinyUI(pageWithSidebar(
  headerPanel("HVTN065 Cytokine Summary"),
  sidebarPanel(
    uiOutput("PTID"),
    selectInput(inputId = "tcells",
                label = "T-Cells",
                choices = c("CD4+" = "cd4", "CD8+" = "cd8")),
    selectInput(inputId = "filter_PTIDs",
              label = "Filter PTIDs by Treatment Status",
              choices = c("Treatment", "Placebo", "Both"),
              selected = "Both"),
    wellPanel(
      p(strong("Options")),
      sliderInput(inputId = "x_range_density", label = "Density Range - x-axis",
                  min = -10, max = 10, value = c(-10, 10), step = 0.25),
      sliderInput(inputId = "y_range_density", label = "Density Range - y-axis",
                  min = 0, max = 2, value = c(0, 1), step = 0.25),
      checkboxInput(inputId = "first_deriv",
                    label = strong("First Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.first_deriv == true",
         sliderInput(inputId = "x_range_first_deriv", label = "Range - x-axis",
                min = -10, max = 10, value = c(-10, 10), step = 0.25),
         sliderInput(inputId = "y_range_first_deriv", label = "Range - y-axis",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)
      ),
      checkboxInput(inputId = "second_deriv",
                    label = strong("Second Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.second_deriv == true",
         sliderInput(inputId = "x_range_second_deriv", label = "Range - x-axis",
                min = -10, max = 10, value = c(-10, 10), step = 0.25),
         sliderInput(inputId = "y_range_second_deriv", label = "Range - y-axis",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)              
      ),
      checkboxInput(inputId = "display_cutpoint",
                    label = strong("Display Cutpoint"),
                    value = FALSE),
      conditionalPanel(condition = "input.display_cutpoint == true",
         sliderInput(inputId = "cutpoint_tol", label = "Cutpoint Tolerance",
                min = -6, max = -1, value = -3, format = "1e#")
      )
    ),
    wellPanel(
      p(strong("Cytokine Gates")),
      checkboxInput(inputId = "TNFa_gates",
                    label = strong("TNFa Gates"),
                    value = FALSE),
      checkboxInput(inputId = "IFNg_gates",
                    label = strong("IFNg Gates"),
                    value = FALSE),
      checkboxInput(inputId = "IL2_gates",
                    label = strong("IL2 Gates"),
                    value = FALSE),
      p(strong("Upstream Gates")),
      checkboxInput(inputId = "singlet_gates",
                    label = strong("Singlet Gates"),
                    value = FALSE),
      checkboxInput(inputId = "viable_gates",
                    label = strong("Viable Gates"),
                    value = FALSE),
      checkboxInput(inputId = "lymph_gates",
                    label = strong("Lymphocyte Gates"),
                    value = FALSE),
      checkboxInput(inputId = "cd3_gates",
                    label = strong("CD3 Gates"),
                    value = FALSE),
      checkboxInput(inputId = "cd4_cd8_gates",
                    label = strong("CD4/CD8 Gates"),
                    value = FALSE)
    )

  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Cytokine Densities",
               plotOutput("densitiesPlot"),
               conditionalPanel(condition = "input.first_deriv == true",
                                p(strong("First Derivatives")),
                                plotOutput("firstDerivPlot")
                                ),
               conditionalPanel(condition = "input.second_deriv == true",
                                p(strong("Second Derivatives")),
                                plotOutput("secondDerivPlot")
                                )
              
      ),
      tabPanel("Cytokine Gates",
               conditionalPanel(condition = "input.TNFa_gates == true",
                                p(strong("TNFa Gates")),
                                plotOutput("TNFaGatesPlot")
                                ),
               conditionalPanel(condition = "input.IFNg_gates == true",
                                p(strong("IFNg Gates")),
                                plotOutput("IFNgGatesPlot")
                                ),
               conditionalPanel(condition = "input.IL2_gates == true",
                                p(strong("IL2 Gates")),
                                plotOutput("IL2GatesPlot")
                                )
      ),
      tabPanel("Upstream Gates",
               conditionalPanel(condition = "input.singlet_gates == true",
                                p(strong("Singlet Gates")),
                                plotOutput("singletGatesPlot")
                                ),
               conditionalPanel(condition = "input.viable_gates == true",
                                p(strong("Viable Gates")),
                                plotOutput("viableGatesPlot")
                                ),
               conditionalPanel(condition = "input.lymph_gates == true",
                                p(strong("Lymphocyte Gates")),
                                plotOutput("lymphGatesPlot")
                                ),

               conditionalPanel(condition = "input.cd3_gates == true",
                                p(strong("CD3 Gates")),
                                plotOutput("cd3GatesPlot")
                                ),
               conditionalPanel(condition = "input.cd4_cd8_gates == true",
                                p(strong("CD4/CD8 Gates")),
                                plotOutput("cd4cd8GatesPlot")
                                )

    )
  ))
))

    
