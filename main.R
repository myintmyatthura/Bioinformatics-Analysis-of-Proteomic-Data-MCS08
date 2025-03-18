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
  
  output$dynamicUI <- renderUI({
    if (user_auth()) {
      navbarPage("Protein Quantitation Analysis",
                 tabPanel("Main", 
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("file", "Upload CSV File", accept = ".csv"),
                              verbatimTextOutput("fileInfo")
                            ),
                            mainPanel(
                              div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                  plotOutput("volcanoPlot")
                              ),
                              div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                  plotOutput("heatmapPlot")
                              ),
                              div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                  tableOutput("filteredResults")
                              ),
                              actionButton("show_more", "Show More"),
                              actionButton("show_less", "Show Less")
                            )
                          )
                 ),
                 navbarMenu("Options",
                            tabPanel(actionButton("logout", "Logout"))
                 )
      )
    } else {
      fluidPage(
        titlePanel("Login"),
        sidebarLayout(
          sidebarPanel(
            textInput("username", "Username"),
            passwordInput("password", "Password"),
            actionButton("login", "Login")
          ),
          mainPanel(
            verbatimTextOutput("login_message")
          )
        )
      )
    }
  })
  
  observeEvent(input$login, {
    if (input$username == "admin" && input$password == "admin") {
      user_auth(TRUE)
    } else {
      output$login_message <- renderText("Invalid username or password")
    }
  })
  
  observeEvent(input$logout, {
    user_auth(FALSE)
  })
  
  dataset <- reactive({
    req(input$file)
    read.csv(input$file$datapath, row.names = 1)
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
  
  processedData <- reactive({
    data <- dataset()
    data <- data[!grepl("immunoglobulin", data$ProteinDescriptions, ignore.case = TRUE), ]
    
    CI_cols <- grep("^CI", colnames(data))
    HE_cols <- grep("^HE", colnames(data))
    CI_data <- data[, CI_cols]
    HE_data <- data[, HE_cols]
    
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
    
    group <- factor(c(rep("CI", 4), rep("HE", 4)))
    design <- model.matrix(~0 + group)
    contrast.matrix <- makeContrasts(groupCI - groupHE, levels = design)
    
    fit <- lmFit(normalized_matrix, design)
    fitlc <- contrasts.fit(fit, contrast.matrix)
    fitlc <- eBayes(fitlc)
    
    result <- topTable(fitlc, adjust = "BH", number = Inf)
    result$Significance <- ifelse(result$adj.P.Val < 0.05 & abs(result$logFC) > 1, 
                                  ifelse(result$logFC > 1, "Upregulated", "Downregulated"), 
                                  "Not Significant")
    
    significant_proteins <- subset(result, Significance != "Not Significant")
    significant_proteins$Protein <- rownames(significant_proteins)
    
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