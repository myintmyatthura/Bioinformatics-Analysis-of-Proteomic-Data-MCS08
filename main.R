library(shiny)
library(shinythemes)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


ui <- fluidPage(
  theme = shinytheme("cerulean"),
  uiOutput("dynamicUI")
)

server <- function(input, output, session) {
  user_auth <- reactiveVal(FALSE)
  tableData <- reactiveVal(NULL)
  login_attempt <- reactiveVal(FALSE)
  statusMessage <- reactiveVal("Please upload a file and click 'Start Analysis'")
  processingMessage <- reactiveVal("Results Not Ready")
  

  
  output$status_message <- renderText({
    message <- statusMessage()
    paste0(message)
  })
  
  output$processingMessage <- renderText({
    message <- processingMessage()
    color <- if (message == "Analysis complete!") "green" else "red"
    paste0(message)
  })
  
  # Ensure UI renders HTML
  outputOptions(output, "status_message", suspendWhenHidden = FALSE)
  outputOptions(output, "processingMessage", suspendWhenHidden = FALSE)
  
  output$dynamicUI <- renderUI({
    if (!login_attempt()) {
      fluidPage(
        titlePanel(div("Proteomic Data Analysis Software", style = "text-align: center; font-size: 32px; font-weight: bold;")),
        br(),
        div(
          actionButton("goToLogin", "Begin", style = "display: block; margin: 0 auto; font-size: 18px; padding: 10px 20px;"),
          style = "text-align: center;"
        ),
        br(), br(),
        div(
          h3("About Us", style = "text-align: center;"),
          p("Our software is designed for quantitative proteomic data analysis, providing visualization and statistical analysis to assist researchers in understanding protein expression changes.", style = "text-align: center;"),
          br(),
          h3("Supported Tools & Formats", style = "text-align: center;"),
          p("We utilize cutting-edge technologies: R, Shiny, ggplot2, limma, and pheatmap.", style = "text-align: center;"),
          p("Supported file formats: CSV, XLSX.", style = "text-align: center;"),
          br(),
          h3("Contributors", style = "text-align: center;"),
          p("Dr. John Doe, Dr. Jane Smith, Alex Johnson.", style = "text-align: center;")
        )
      )
    } else if (!user_auth()) {
      fluidPage(
        titlePanel("Login"),
        sidebarLayout(
          sidebarPanel(
            textInput("username", "Username"),
            passwordInput("password", "Password"),
            actionButton("login", "Login"),
            textOutput("login_message")
          ),
          mainPanel(
            h3("Please enter your credentials to continue.")
          )
        )
      )
    } else {
      navbarPage("Protein Quantitation Analysis",
                 tabPanel("Analysis", 
                          sidebarLayout(
                            sidebarPanel(
                              textOutput("status_message", inline = TRUE),
                              br(),
                              fileInput("file", "Upload CSV File", accept = ".csv"),
                              verbatimTextOutput("fileInfo"),
                              textInput("ci_columns", "Enter CI Columns (comma-separated):", ""),
                              textInput("he_columns", "Enter HE Columns (comma-separated):", ""),
                              textInput("log2_threshold", "Enter Log2 Fold Change Threshold (eg: 1):", ""),
                              textInput("pval_threshold", "Enter Adjusted P-Value Threshold (eg: 0.05):", ""),
                              actionButton("start_analysis", "Start Analysis"),
                              textOutput("processingMessage", inline = TRUE),
                            )
                            ,
                            mainPanel(
                              uiOutput("buttons"),
                            
                              br(), br(),
                              
                              # Volcano Plot Section
                              conditionalPanel(
                                condition = "input.toggle_volcano % 2 == 1",
                                div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                    h4("Volcano Plot"),
                                    p("The volcano plot visualizes the relationship between the log fold change and adjusted p-value of proteins. 
           Points in red indicate significantly upregulated proteins, while points in blue indicate significantly downregulated proteins. 
           Gray points are not significant."),
                                    plotOutput("volcanoPlot")
                                )
                              ),
                              
                              # Heatmap Section
                              conditionalPanel(
                                condition = "input.toggle_heatmap % 2 == 1",
                                div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                    h4("Heatmap Plot"),
                                    p("The heatmap represents the expression levels of significantly different proteins. 
           Red indicates higher expression, blue indicates lower expression, and clustering shows relationships among proteins."),
                                    plotOutput("heatmapPlot")
                                )
                              ),
                              
                              # Filtered Table Section
                              conditionalPanel(
                                condition = "input.toggle_table % 2 == 1",
                                div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                    h4("Filtered Results Table"),
                                    p("This table lists proteins that are significantly different between conditions based on the specified thresholds. 
           It includes log fold change, adjusted p-values, and classification as upregulated, downregulated, or not significant."),
                                    tableOutput("filteredResults"),
                                    actionButton("show_more", "Show More"),
                                    actionButton("show_less", "Show Less")
                                )
                              )
                            )
                            
                          )
                 ),
                 navbarMenu("Options",
                            tabPanel(actionButton("logout", "Logout"))
                 ),
                 tabPanel("About", 
                          fluidPage(
                            titlePanel(div("Proteomic Data Analysis Software", style = "text-align: center; font-size: 32px; font-weight: bold;")),
                            br(), br(),
                            div(
                              h3("About Us", style = "text-align: center;"),
                              p("Our software is designed for quantitative proteomic data analysis, providing visualization and statistical analysis to assist researchers in understanding protein expression changes.", style = "text-align: center;"),
                              br(),
                              h3("Supported Tools & Formats", style = "text-align: center;"),
                              p("We utilize cutting-edge technologies: R, Shiny, ggplot2, limma, and pheatmap.", style = "text-align: center;"),
                              p("Supported file formats: CSV, XLSX.", style = "text-align: center;"),
                              br(),
                              h3("Contributors", style = "text-align: center;"),
                              p("Dr. John Doe, Dr. Jane Smith, Alex Johnson.", style = "text-align: center;")
                            )
                          )
                 )
      )
    }
  })
  output$buttons <- renderUI({
    message <- statusMessage()
    processing <- processingMessage()
    
    # Choose color based on status message
    btn_color <- if (message == "Analysis complete!") {
      "#009914"
    } else if (message == "Please upload a file and click 'Start Analysis'") {
      "#ff0000"
    } else {
      "#FFFFFF"  # Default (white)
    }
    
    # Disable buttons if processing is not "Done!"
    is_disabled <- processing != "Done!"
    
    # Wrap buttons in a div with inline-block styling
    div(
      style = "display: flex; gap: 10px;",  # Arrange buttons side by side with spacing
      actionButton("toggle_volcano", "Volcano Plot", 
                   style = paste("background-color:", btn_color, " !important;",
                                 "border-color:", btn_color, " !important;",
                                 "color: black; width: 150px; height: 40px;"),
                   disabled = is_disabled),
      
      actionButton("toggle_heatmap", "Heatmap Plot", 
                   style = paste("background-color:", btn_color, " !important;",
                                 "border-color:", btn_color, " !important;",
                                 "color: black; width: 150px; height: 40px;"),
                   disabled = is_disabled),
      
      actionButton("toggle_table", "Filtered Table", 
                   style = paste("background-color:", btn_color, " !important;",
                                 "border-color:", btn_color, " !important;",
                                 "color: black; width: 150px; height: 40px;"),
                   disabled = is_disabled)
    )
  })
  
  
  
  
  
  observeEvent(input$goToLogin, {
    login_attempt(TRUE)
  })
  
  observeEvent(input$login, {
    if (input$username == "admin" && input$password == "admin") {
      user_auth(TRUE)
    } else {
      output$login_message <- renderText("Invalid username or password. Please try again.")
    }
  })
  
  observeEvent(input$logout, {
    user_auth(FALSE)
    login_attempt(FALSE)
  })
  
  
  dataset <- reactive({
    req(input$file)
    read.csv(input$file$datapath, row.names = 1, check.names = FALSE)
  })
  
  output$fileInfo <- renderText({
    req(input$file)
    file <- input$file
    data <- dataset()
    paste("File Name:", file$name, "\n",
          "File Type:", file$type, "\n",
          "File Size:", file$size, "bytes\n",
          "Row Count:", nrow(data), "\n",
          "Column Names:", paste(colnames(data), collapse = ", "))
  })
  
  observeEvent(input$start_analysis, {
    statusMessage("Running analysis... Please wait.")
    processingMessage("Running")
    processedData()  # Run the analysis
    Sys.sleep(2)  # Simulate delay
  })
  
  

  
  
  processedData <- eventReactive(input$start_analysis, {
    
    req(input$file, input$ci_columns, input$he_columns, input$log2_threshold, input$pval_threshold)
    log2_threshold <- as.numeric(input$log2_threshold)
    pval_threshold <- as.numeric(input$pval_threshold)
    data <- dataset()
    data <- data[!grepl("immunoglobulin", data$ProteinDescriptions, ignore.case = TRUE), ]
    
    ci_cols <- strsplit(input$ci_columns, ",")[[1]]
    he_cols <- strsplit(input$he_columns, ",")[[1]]
    ci_cols <- trimws(ci_cols)
    he_cols <- trimws(he_cols)
    
    
    if (!all(ci_cols %in% colnames(data)) || !all(he_cols %in% colnames(data))) {
      showNotification("Invalid column names. Please check and try again.", type = "error")
      return(NULL)
    }
    
    CI_data <- data[, ci_cols, drop = FALSE]
    HE_data <- data[, he_cols, drop = FALSE]
    
    impute_mean <- function(group_data) {
      imputed <- apply(group_data, 1, function(impute) {
        impute[is.na(impute)] <- mean(impute, na.rm = TRUE) 
        return(impute)
      })
      return(t(imputed))
    }
    
    CI_data_imputed <- impute_mean(CI_data)
    HE_data_imputed <- impute_mean(HE_data)
    data_preprocessed <- na.omit(cbind(CI_data_imputed, HE_data_imputed))
    
    data_matrix <- as.matrix(data_preprocessed)
    log_data_matrix <- log2(data_matrix + 1)
    normalized_matrix <- normalizeQuantiles(log_data_matrix)
    
    group <- factor(c(rep("CI", length(ci_cols)), rep("HE", length(he_cols))))
    design <- model.matrix(~0 + group)
    contrast.matrix <- makeContrasts(groupCI - groupHE, levels = design)
    
    fit <- lmFit(normalized_matrix, design)
    fitlc <- contrasts.fit(fit, contrast.matrix)
    fitlc <- eBayes(fitlc)
    
    result <- topTable(fitlc, adjust = "BH", number = Inf)
    
    result$Protein <- rownames(result)
    result <- result[, c("Protein", setdiff(colnames(result), "Protein"))] 
    
    result$Significance <- ifelse(result$adj.P.Val < pval_threshold & abs(result$logFC) > log2_threshold, 
                                  ifelse(result$logFC > log2_threshold, "Upregulated", "Downregulated"), 
                                  "Not Significant")
    
    significant_proteins <- subset(result, Significance != "Not Significant")
    statusMessage("Analysis complete!")
    processingMessage("Done!")
    tableData(head(result, 20))
    list(normalized_matrix = normalized_matrix, result = result, significant_proteins = significant_proteins)

  })

  output$volcanoPlot <- renderPlot({
    req(processedData())
    result <- processedData()$result
    
    ggplot(result, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
      theme_minimal() +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
  })
  
  output$heatmapPlot <- renderPlot({
    req(processedData())
    significant_proteins <- processedData()$significant_proteins
    filtered_expression_matrix <- processedData()$normalized_matrix[significant_proteins$Protein, ]
    
    pheatmap(
      filtered_expression_matrix, 
      scale = "row",
      clustering_distance_rows = "euclidean", 
      clustering_method = "ward.D2",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1.5, 1.5, length.out = 101),  # Ensure color mapping covers range
      legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),  # Custom legend values
      legend_labels = c("-1.5 (Strongly Down)", "-1", "-0.5", "0 (No Change)", "0.5", "1", "1.5 (Strongly Up)"),
      main = "logFC Heatmap with Clustering"
    )
    
  })
  
  output$filteredResults <- renderTable({
    req(tableData())
    tableData()
  })
  
  observeEvent(input$show_more, {
    tableData(processedData()$result)
  })
  
  observeEvent(input$show_less, {
    tableData(head(processedData()$result, 20))
  })
  
}

shinyApp(ui = ui, server = server)
