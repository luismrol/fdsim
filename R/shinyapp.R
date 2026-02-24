library(shiny)

ui <- fluidPage(
  titlePanel("Multivariate Functional Data Simulator"),
  sidebarLayout(
    sidebarPanel(
      numericInput("ncurves", "Number of curves (N):", 10, min = 1, max = 200),
      selectInput("method", "Method:", choices = c("svd", "chol", "eigen")),
      selectInput("covar", "Covariance type:", choices = c("sq", "abs"), multiple = FALSE),
      numericInput("rho", "Correlation (rho):", 0.5, min = -1, max = 1, step = 0.1),
      textInput("mu1_expr", "mu1 (as function of grid):", "-"),
      textInput("mu2_expr", "mu2 (as function of grid):", "-"),
      actionButton("simulate", "Simulate")
    ),
    mainPanel(
      plotOutput("simplot")
    )
  )
)

server <- function(input, output) {
  library(fdsim)
  browser()
  sim_data <- eventReactive(input$simulate, {
    # Evaluate user-specified means safely
    grid <- seq(0, 1, length.out = 100)
    mu1 <- eval(parse(text = input$mu1_expr), envir = list(grid = grid))
    mu2 <- eval(parse(text = input$mu2_expr), envir = list(grid = grid))
    mu <- rbind(t(mu1), t(mu2))
    mfd_sim(input$ncurves, mu, method = input$method,
            covar = input$covar, rho = input$rho)
  })

  output$simplot <- renderPlot({
    Y <- sim_data()
    plot(x = grid, y = Y[[1]][1,], col = "white",
         main = "Multivariate simulated data", ylab = "value",
         ylim = c(min(Y[[1]]-2, Y[[2]]), 2+max(Y[[1]], Y[[2]])))

    for (i in 1:nrow(Y[[1]])) lines(grid, Y[[1]][i,], col = "gray")
    for (i in 1:nrow(Y[[2]])) lines(grid, Y[[2]][i,], col = "pink")

    # replot mean functions
    lines(grid, mu[1,], col = "black")
    lines(grid, mu[2,], col = "red")
  })
}

shinyApp(ui, server)
