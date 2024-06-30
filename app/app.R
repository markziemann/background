library("clusterProfiler")
library("markdown")
library("rmarkdown")
library("shiny")
library("broom")
library("RhpcBLASctl")
library("eulerr")

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
        tabPanel("Instructions", includeMarkdown("intro.md")),
        tabPanel("Data",
                 "Foreground", verbatimTextOutput("contents1"),
                 "Background",verbatimTextOutput("contents2"),
                 "No. genes in foregorund", verbatimTextOutput("summary1"),
                 "No. genes in background", verbatimTextOutput("summary2")),
        tabPanel("Comparative Analysis", 
                 textOutput("Comparison of original and BG fixed analysis - top 20 pathways with divergent FDR values"),
                 tableOutput("tbl1")),
        tabPanel("Charts", textOutput("counts"),
                 plotOutput("euler"),
                 plotOutput("scatter_es",height=550),
                 plotOutput("scatter_fdr",height=550)),
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
    selectInput("genesetlibrary", "Gene set library", 
                choices = c("CellMarkers", "Gene Ontology","GTRD TF Targets",
                            "Hallmark","Human Phenotype Ontology",
                            "KEGG Pathways","miR targets",
                            "Reactome Pathways","WikiPathways"))
  })
  
  mygs <- reactive({
    gs <- readRDS("genesets/gs.Rds")
    if ( input$genesetlibrary == "CellMarkers") {
      mygs <- gs[["CellMarkers"]]
    }
    if ( input$genesetlibrary == "Gene Ontology") {
      mygs <- gs[["Gene Ontology"]]
    }
    if ( input$genesetlibrary == "GTRD TF Targets") {
      mygs <- gs[["GTRD TF Targets"]]
    }
    if ( input$genesetlibrary == "Hallmark") {
      mygs <- gs[["Hallmark"]]
    }
    if ( input$genesetlibrary == "Human Phenotype Ontology") {
      mygs <- gs[["Human Phenotype Ontology"]]
    }
    if ( input$genesetlibrary == "KEGG Pathways") {
      mygs <- gs[["KEGG Pathways"]]
    }
    if ( input$genesetlibrary == "miR targets") {
      mygs <- gs[["miR targets"]]
    }
    if ( input$genesetlibrary == "Reactome Pathways") {
      mygs <- gs[["Reactome Pathways"]]
    }
    if ( input$genesetlibrary == "WikiPathways") {
      mygs <- gs[["WikiPathways"]]
    }
    return(mygs)
  })
    
  original <- reactive({
    req(fg())
    req(bg())
    req(mygs())
    options(enrichment_force_universe = FALSE)
    ora <- as.data.frame(enricher(gene = fg() ,
                                     universe = bg(), minGSSize = 5, maxGSSize = 500000, TERM2GENE = mygs(),
                                     pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))
    
    gr <- as.numeric(sapply(strsplit(ora$GeneRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora$GeneRatio,"/"),"[[",2))
    
    br <- as.numeric(sapply(strsplit(ora$BgRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora$BgRatio,"/"),"[[",2))
    
    ora$ES <- gr/br
    ora$ID = ora$geneID = ora$p.adjust = ora$Count = NULL
    colnames(ora) <- gsub("qvalue","FDR",colnames(ora))
    colnames(ora) <- gsub("GeneRatio","FgRatio",colnames(ora))
    ora <- ora[,c("Description","FgRatio","BgRatio","ES","pvalue","FDR")]
    return(ora)
  })
  
  bgfix <- reactive({
    req(fg())
    req(bg())
    req(mygs())
    options(enrichment_force_universe = TRUE)
    ora <- as.data.frame(enricher(gene = fg() ,
                                  universe = bg(), minGSSize = 5, maxGSSize = 500000, TERM2GENE = mygs(),
                                  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))
    
    gr <- as.numeric(sapply(strsplit(ora$GeneRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora$GeneRatio,"/"),"[[",2))
    
    br <- as.numeric(sapply(strsplit(ora$BgRatio,"/"),"[[",1)) /
      as.numeric(sapply(strsplit(ora$BgRatio,"/"),"[[",2))
    
    ora$ES <- gr/br
    ora$ID = ora$geneID = ora$p.adjust = ora$Count = NULL
    colnames(ora) <- gsub("qvalue","FDR",colnames(ora))
    colnames(ora) <- gsub("GeneRatio","FgRatio",colnames(ora))
    ora <- ora[,c("Description","FgRatio","BgRatio","ES","pvalue","FDR")]
    return(ora)
  })

  bgfixtbl <- reactive({
    orig_df <- original()
    bgfix_df <- bgfix()
    m <- merge(orig_df,bgfix_df,by="Description")
    m <- m[,c("Description","FgRatio.x","FgRatio.y","BgRatio.x","BgRatio.y","ES.x","ES.y","FDR.x","FDR.y")]
    diff <- abs(m$FDR.x - m$FDR.y)
    m <- m[order(-diff),]
    return(m)
  })
  
  output$tbl1 <- renderTable({
    if ( input$comparison == "Background error" ) {
      tbl <- head(bgfixtbl(),20)
    }
    return(tbl)
  })  
  
  output$counts <- renderText({
    oricnt <- nrow(subset(original(),FDR<0.05))
    bgfixcnt <- nrow(subset(bgfix(),FDR<0.05))
    out <- paste("Original:",oricnt,", and BG fixed:",bgfixcnt,"@FDR<0.05")
    return(out)
  })
  
  output$euler <- renderPlot({
    orig_df <- original()
    bgfix_df <- bgfix()
    orig_sets <- subset(orig_df,FDR < 0.05)$Description
    bgfix_sets <- subset(bgfix_df,FDR < 0.05)$Description
    v1 <- list("Original"=orig_sets, "BG fix"=bgfix_sets)
    plot(euler(v1),quantities = list(cex = 2), labels = list(cex = 2))
  })
  
  output$scatter_es <- renderPlot({
    orig_df <- original()
    bgfix_df <- bgfix()
    m <- merge(orig_df,bgfix_df,by="Description")
    m <- m[,c("ES.x","ES.y")]
    MAX=max(c(m$ES.x,m$ES.y))
    plot(m$ES.x,m$ES.y,xlim=c(0,MAX),ylim=c(0,MAX),
         xlab="Original",ylab="BG Fix")
    mtext("Fold Enrichment Scores")
    abline(a = 0, b = 1,lwd=2,lty=2,col="red")
  })

  output$scatter_fdr <- renderPlot({
    orig_df <- original()
    bgfix_df <- bgfix()
    m <- merge(orig_df,bgfix_df,by="Description")
    m <- m[,c("FDR.x","FDR.y")]
    MAX=max(c(-log10(m$FDR.x),-log10(m$FDR.y)))
    plot(-log10(m$FDR.x),-log10(m$FDR.y),
         xlim=c(0,MAX),ylim=c(0,MAX),
         xlab="Original",ylab="BG Fix")
    mtext("FDR Values")
    abline(a = 0, b = 1,lwd=2,lty=2,col="red")
  })
  
  plot_data <- reactive({
    req(data(), input$x_axis, input$y_axis)
    ggplot(data(), aes_string(x = input$x_axis, y = input$y_axis)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(x = input$x_axis, y = input$y_axis) +
      ggtitle(paste(input$y_axis, "vs", input$x_axis))
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