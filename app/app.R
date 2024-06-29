library(tidyverse)
library(ggplot2)
library(rmarkdown)
library(broom)

# Define UI
ui <- fluidPage(
  titlePanel("Two subtle issues with over-representation analysis: Demonstration"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fg", "Choose foreground gene list",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".txt")),
      fileInput("bg", "Choose background gene list",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".txt")),
      uiOutput("comparison"),
      uiOutput("genesetlibrary"),
      downloadButton("download_report", "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data",
                 "Foreground", verbatimTextOutput("contents1"),
                 "Background",verbatimTextOutput("contents2"),
                 "No. genes in foregorund", verbatimTextOutput("summary1"),
                 "No. genes in background", verbatimTextOutput("summary2")),
        tabPanel("Scatterplot", plotOutput("scatterplot")),
        tabPanel("Regression Summary", verbatimTextOutput("regression_summary"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  fg <- reactive({
    req(input$fg)
    fg <- readLines(input$fg$datapath)
    return(fg)
  })
  
  bg <- reactive({
    req(input$bg)
    bg <- readLines(input$bg$datapath)
    return(bg)
  })
  
  output$contents1 <- renderPrint({
    req(fg())
    head(fg())
  })
  
  output$contents2 <- renderPrint({
    req(bg())
    head(bg())
  })
  
  output$summary1 <- renderPrint({
    req(fg())
    length(fg())
  })
  
  output$summary2 <- renderPrint({
    req(bg())
    length(bg())
  })
  
  output$comparison <- renderUI({
    selectInput("comparison", "Error Comparison", choices = c("Background error","FDR error","Both errors"), selected = "both")
  })
  
  output$genesetlibrary <- renderUI({
    selectInput("genesetlibrary", "Gene set library", choices = c("Gene Ontology","Reactome Pathways","KEGG Pathways"), selected = names(data())[2])
  })
  
  plot_data <- reactive({
    req(data(), input$x_axis, input$y_axis)
    ggplot(data(), aes_string(x = input$x_axis, y = input$y_axis)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(x = input$x_axis, y = input$y_axis) +
      ggtitle(paste(input$y_axis, "vs", input$x_axis))
  })
  
  output$scatterplot <- renderPlot({
    plot_data()
  })
  
  regression_model <- reactive({
    req(data(), input$x_axis, input$y_axis)
    lm(formula = as.formula(paste(input$y_axis, "~", input$x_axis)), data = data())
  })
  
  regression_stats <- reactive({
    req(regression_model())
    model <- regression_model()
    summary_model <- summary(model)
    
    # Calculate correlation and its p-value
    cor_test <- cor.test(data()[[input$x_axis]], data()[[input$y_axis]])
    
    list(
      r_value = sqrt(summary_model$r.squared) * sign(summary_model$coefficients[2,1]),
      r_squared = summary_model$r.squared,
      adj_r_squared = summary_model$adj.r.squared,
      f_statistic = summary_model$fstatistic[1],
      p_value = summary_model$coefficients[2,4],
      correlation = cor_test$estimate,
      cor_p_value = cor_test$p.value
    )
  })
  
  output$regression_summary <- renderPrint({
    req(regression_model(), regression_stats())
    
    cat("Regression Summary:\n\n")
    cat("R value:", round(regression_stats()$r_value, 4), "\n")
    cat("R-squared:", round(regression_stats()$r_squared, 4), "\n")
    cat("Adjusted R-squared:", round(regression_stats()$adj_r_squared, 4), "\n")
    cat("F-statistic:", round(regression_stats()$f_statistic, 4), "\n")
    cat("P-value:", format.pval(regression_stats()$p_value, digits = 4), "\n\n")
    
    cat("Model Coefficients:\n")
    print(summary(regression_model())$coefficients)
    
    cat("\nAdditional Statistics:\n")
    cat("Correlation coefficient:", round(regression_stats()$correlation, 4), "\n")
    cat("Correlation p-value:", format.pval(regression_stats()$cor_p_value, digits = 4), "\n")
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".html", sep="")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      params <- list(plot = plot_data(),
                     stats = regression_stats(),
                     x_var = input$x_axis,
                     y_var = input$y_axis,
                     coefficients = summary(regression_model())$coefficients)
      
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)