shinyUI(pageWithSidebar(
  headerPanel("HVTN065 Cytokine Summary"),
  sidebarPanel(
    selectInput(inputId = "PTID",
                label = "PTID",
                choices = levels_PTID), #, selected = "123290487"
    selectInput(inputId = "tcells",
                label = "T-Cells",
                choices = c("CD4+" = "cd4", "CD8+" = "cd8")),
    wellPanel(
      p(strong("Options")),
      sliderInput(inputId = "max_density", label = "Maximum Density",
                  min = 0.25, max = 4, value = 1.0, step = 0.25),
      checkboxInput(inputId = "first_deriv",
                    label = strong("First Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.first_deriv == true",
         sliderInput(inputId = "range_first_deriv", label = "Limit Range",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)
      ),
      checkboxInput(inputId = "second_deriv",
                    label = strong("Second Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.second_deriv == true",
         sliderInput(inputId = "range_second_deriv", label = "Limit Range",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)
      ),
      checkboxInput(inputId = "gates",
                    label = strong("Upstream Gates"),
                    value = FALSE)
    )
  ),
  mainPanel(
    plotOutput("densitiesPlot"),
    conditionalPanel(condition = "input.first_deriv == true",
      p(strong("First Derivatives")),
      plotOutput("firstDerivPlot")
    ),
    conditionalPanel(condition = "input.second_deriv == true",
      p(strong("Second Derivatives")),
      plotOutput("secondDerivPlot")
    ),
    conditionalPanel(condition = "input.gates == true",
      p(strong("Upstream Gates (Inactive currently)")),
      plotOutput("gatesPlot")
    )
  )
))
