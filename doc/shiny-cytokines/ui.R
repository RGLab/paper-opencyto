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
      sliderInput(inputId = "x_range_density", label = "Density Range - x-axis",
                  min = -10, max = 10, value = c(-5, 5), step = 0.25),
      sliderInput(inputId = "y_range_density", label = "Density Range - y-axis",
                  min = 0, max = 2, value = c(0, 1), step = 0.25),
      checkboxInput(inputId = "first_deriv",
                    label = strong("First Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.first_deriv == true",
         sliderInput(inputId = "x_range_first_deriv", label = "Range - x-axis",
                min = -10, max = 10, value = c(-5, 5), step = 0.25),
         sliderInput(inputId = "y_range_first_deriv", label = "Range - y-axis",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)
      ),
      checkboxInput(inputId = "second_deriv",
                    label = strong("Second Derivative"),
                    value = FALSE),
      conditionalPanel(condition = "input.second_deriv == true",
         sliderInput(inputId = "x_range_second_deriv", label = "Range - x-axis",
                min = -10, max = 10, value = c(-5, 5), step = 0.25),
         sliderInput(inputId = "y_range_second_deriv", label = "Range - y-axis",
                min = -3, max = 3, value = c(-1, 1), step = 0.25)              
      ),
      checkboxInput(inputId = "display_cutpoint",
                    label = strong("Display Cutpoint"),
                    value = FALSE),
      conditionalPanel(condition = "input.display_cutpoint == true",
         sliderInput(inputId = "cutpoint_tol", label = "Cutpoint Tolerance",
                min = -6, max = -1, value = -3, format = "1e#")
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
