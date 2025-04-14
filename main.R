custom_lib <- "/srv/shiny-server/Bioinformatics-Analysis-of-Proteomic-Data-MCS08/r_libs"
if (dir.exists(custom_lib)) {
  .libPaths(c(custom_lib, .libPaths()))
}

library(shiny)
library(shinythemes)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(httr)
library(jsonlite)
library(DBI)
library(pool)
library(rmarkdown)
library(knitr)
library(tinytex)
library(kableExtra)

config <- fromJSON(readLines("config.json"))


# Create the pool connection using the config values
pool <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = config$db$dbname,
  host = config$db$host,
  port = config$db$port,
  user = config$db$user,
  password = config$db$password
)

initializeDB <- function(pool) {
  
  
  # Create users table
  dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS users (
      user_id SERIAL PRIMARY KEY,
      username VARCHAR(50) UNIQUE NOT NULL,
      password_hash VARCHAR(255) NOT NULL,
      email VARCHAR(100) UNIQUE NOT NULL,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
  ")
  
  # Create analysis_history table
  dbExecute(pool, "
    CREATE TABLE IF NOT EXISTS analysis_history (
      analysis_id SERIAL PRIMARY KEY,
      user_id INTEGER REFERENCES users(user_id),
      filename VARCHAR(255) NOT NULL,
      s3_file_key VARCHAR(255) NOT NULL,
      ci_columns TEXT,
      he_columns TEXT,
      log2_threshold FLOAT,
      pval_threshold FLOAT,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      result_s3_key VARCHAR(255)
    )
  ")
}

initializeDB(pool)

ui <- fluidPage(
  theme = shinytheme("cerulean"),
  uiOutput("dynamicUI")
)

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = config$db$dbname,
  host = config$db$host,
  port = config$db$port,
  user = config$db$user,
  password = config$db$password
)
dbListTables(con)



server <- function(input, output, session) {
  user_auth <- reactiveVal(FALSE)
  tableData <- reactiveVal(NULL)
  login_attempt <- reactiveVal(FALSE)
  show_register <- reactiveVal(FALSE)
  statusMessage <- reactiveVal("Please upload a file and click 'Start Analysis'")
  processingMessage <- reactiveVal("Results Not Ready")
  pubmedLinks <- reactiveVal(NULL)
  pubMessage <- reactiveVal("Not ready!")
  
  
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
      if (!show_register()){
        fluidPage(
          titlePanel("Login"),
          sidebarLayout(
            sidebarPanel(
              textInput("username", "Username"),
              passwordInput("password", "Password"),
              actionButton("login", "Login"),
              textOutput("login_message"),
              hr(),
              actionButton("show_register_form", "Create New Account")
            ),
            mainPanel(
              h3("Please enter your credentials to continue.")
            )
          )
        )
      } else {
        fluidPage(
          titlePanel("Register New Account"),
          sidebarLayout(
            sidebarPanel(
              textInput("reg_username", "Username"),
              textInput("reg_email", "Email"),
              passwordInput("reg_password","Password"),
              passwordInput("reg_confirm_password","Confirm Password"),
              actionButton("register", "Register"),
              textOutput("register_message"),
              hr(),
              actionButton("back_to_login", "Back To Login")
            ),
            mainPanel(
              h3("Create a new account to get started")
            )
          )
        )
      }
    } else {
      navbarPage("Protein Quantitation Analysis",
                 tabPanel("Analysis", 
                          sidebarLayout(
                            sidebarPanel(
                              textOutput("status_message", inline = TRUE),
                              br(),
                              fileInput("file", "Upload CSV File", accept = ".csv"),
                              verbatimTextOutput("fileInfo"),
                              textInput("protein_columns", "Enter Protein Name columns (eg: 'Protein Name')", ""),
                              textInput("gene_columns", "Enter Gene Name columns (eg: 'Gene Name')", ""),
                              textInput("ci_columns", "Enter Patient Sample Columns (comma-separated):", ""),
                              textInput("he_columns", "Enter Healthy Control Sample Columns (comma-separated):", ""),
                              textInput("log2_threshold", "Enter Log2 Fold Change Threshold (eg: 1):", ""),
                              textInput("pval_threshold", "Enter Adjusted P-Value Threshold (eg: 0.05):", ""),
                              actionButton("start_analysis", "Start Analysis"),
                              textOutput("processingMessage", inline = TRUE),
                              br(), br(),
                              downloadButton("download_report", "Download Full Report (PDF)")
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
                                    plotOutput("volcanoPlot"),
                                    downloadButton("download_volcano", "Download Volcano Plot (PNG)")
                                    
                                    
                                )
                              ),
                              
                              # Heatmap Section
                              conditionalPanel(
                                condition = "input.toggle_heatmap % 2 == 1",
                                div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                    h4("Heatmap Plot"),
                                    p("The heatmap represents the expression levels of significantly different proteins. 
           Red indicates higher expression, blue indicates lower expression, and clustering shows relationships among proteins."),
                                    plotOutput("heatmapPlot"),
                                    downloadButton("download_heatmap", "Download Heatmap Plot (PNG)")
                                    
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
                                    actionButton("show_less", "Show Less"),
                                    br(),
                                    downloadButton("download_table", "Download Table (CSV)")
                                    
                                )
                              ),
                              
                              conditionalPanel(
                                condition = "input.fetch_pubmed % 2 == 1",
                                div(style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;", 
                                    h4("PubMed Article Links"),
                                    p("Below are links to PubMed articles related to studying the correlation between each significant protein and strokes."),
                                    uiOutput("pubmedUI")
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
    pubmessage <- pubMessage()
    
    # Choose color based on status message
    btn_color <- if (message == "Analysis complete!") {
      "#009914"
    } else if (message == "Please upload a file and click 'Start Analysis'") {
      "#ff0000"
    } else {
      "#FFFFFF"  # Default (white)
    }
    
    pub_color <- if (pubmessage == "Ok!") {
      "#009914"
    } else {
      "#ff0000"
    }
    
    
    # Disable buttons if processing is not "Done!"
    is_disabled <- processing != "Done!"
    pub_disb <- pubmessage != "Ok!"
    
    # Wrap buttons in a div with inline-block styling
    div(
      style = "display: flex; gap: 10px; flex-wrap: wrap;",  # Responsive layout with spacing
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
                   disabled = is_disabled),
      
      actionButton("fetch_pubmed", "PubMed Articles", 
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
    user_result <- verifyLogin(input$username, input$password)
    if(!is.null(user_result)) {
      user_id <- user_result
      user_auth(TRUE)
    } else {
      output$login_message <- renderText("Invalid username or password. Please try again.")
    }
  })
  
  observeEvent(input$show_register_form, {
    show_register(TRUE)
  })
  
  observeEvent(input$back_to_login, {
    show_register(FALSE)
  })
  
  observeEvent(input$logout, {
    user_auth(FALSE)
    login_attempt(FALSE)
  })
  
  
  dataset <- reactive({
    req(input$file, input$protein_columns)
    
    # Read just the column names
    column_names <- names(read.csv(input$file$datapath, nrows = 1, check.names = FALSE))
    
    # Find the index of the column that matches the user input
    col_num <- which(column_names == input$protein_columns)
    

    # Use the column number for row.names
    read.csv(input$file$datapath, row.names = col_num, check.names = FALSE)
  })
  
  
  output$fileInfo <- renderText({
    req(input$file)
    file <- input$file
    
    # Read the CSV without setting row names
    raw_data <- read.csv(file$datapath, check.names = FALSE)
    
    paste("File Name:", file$name, "\n",
          "File Type:", file$type, "\n",
          "File Size:", file$size, "bytes\n",
          "Row Count:", nrow(raw_data), "\n",
          "Column Names:", paste(colnames(raw_data), collapse = ", "))
  })
  
  
  observeEvent(input$start_analysis, {
    statusMessage("Running analysis... Please wait.")
    processingMessage("Running")
    processedData()  # Run the analysis
    
  })
  
  observe({
    if (processingMessage() == "Done!") {
      # Wait 10 seconds, then enable the PubMed button
      shinyjs::delay(10000, {
        pubMessage("Ok!")
      })
    }
  })
  
  
  processedData <- eventReactive(input$start_analysis, {
    
    req(input$file,input$protein_columns,input$gene_columns, input$ci_columns, input$he_columns, input$log2_threshold, input$pval_threshold)
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
    gene_col <- input$gene_columns
    
    
    result$Gene <- data[rownames(result), gene_col]
    
    result <- result[, c("Protein", "Gene", setdiff(colnames(result), c("Protein", "Gene")))]
    
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
  
  getPubMedLinks <- function(protein_name) {
    query <- URLencode(paste0(protein_name, "protein stroke"))
    esearch_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                          "db=pubmed&retmode=json&retmax=5&term=", query)
    
    search_response <- tryCatch({ GET(esearch_url) }, error = function(e) return(NULL))
    if (is.null(search_response) || http_error(search_response)) return(NULL)
    
    search_result <- content(search_response, as = "parsed", type = "application/json")
    ids <- search_result$esearchresult$idlist
    if (length(ids) == 0) return(NULL)
    
    # Fetch titles using esummary
    esummary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
                           "db=pubmed&retmode=json&id=", paste(ids, collapse = ","))
    
    summary_response <- tryCatch({ GET(esummary_url) }, error = function(e) return(NULL))
    if (is.null(summary_response) || http_error(summary_response)) return(NULL)
    
    summary_data <- content(summary_response, as = "parsed", type = "application/json")
    
    # Combine title + link
    links <- lapply(ids, function(id) {
      title <- summary_data$result[[id]]$title
      url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", id)
      list(title = title, url = url)
    })
    
    return(links)
  }
  
  
  observeEvent(input$fetch_pubmed, {
    req(processedData())
    significant_proteins <- processedData()$significant_proteins$Protein
    pubmed_links_list <- lapply(significant_proteins, function(prot) {
      links <- getPubMedLinks(prot)
      list(protein = prot, links = links)
    })
    pubmedLinks(pubmed_links_list)
  })
  
  output$pubmedUI <- renderUI({
    req(pubmedLinks())
    tagList(
      h4("PubMed Links for Stroke-Related Proteins"),
      lapply(pubmedLinks(), function(entry) {
        if (is.null(entry$links)) return(NULL)
        tagList(
          tags$h5(entry$protein),
          tags$ul(
            lapply(entry$links, function(link_entry) {
              tags$li(
                tags$a(href = link_entry$url, target = "_blank", link_entry$title)
              )
            })
          )
        )
      })
    )
  })
  
  
  observeEvent(input$register, {
    if(input$reg_password == input$reg_confirm_password) {
      success <- registerUser(input$reg_username, input$reg_password, input$reg_email)
      if(success) {
        output$register_message <- renderText("Registration successful! Please log in.")
        show_register(FALSE) # Switch back to login view
      } else {
        output$register_message <- renderText("Username or email already exists!")
      }
    } else {
      output$register_message <- renderText("Passwords do not match!")
    }
  })
  
  # user registration function
  registerUser <- function(username, password, email){
    
    # if username exist
    user_exists <- dbGetQuery(pool, sprintf("SELECT 1 FROM users WHERE username ='%s' OR email ='%s' LIMIT 1", username, email))
    if (nrow((user_exists))>0){
      return (FALSE)
    }
    
    password_hash <- digest::digest(password, algo="sha256")
    
    user_data <- dbGetQuery(pool, sprintf("SELECT user_id FROM users WHERE username ='%s' AND password_hash = '%s'", username, password_hash))
    dbExecute(pool, sprintf("INSERT INTO users (username, password_hash, email) VALUES ('%s', '%s', '%s')",
                            username, password_hash, email))
    return(TRUE)
  }
  
  # User login function
  verifyLogin <- function(username, password) {
    # Hash the provided password
    password_hash <- digest::digest(password, algo = "sha256")
    
    # Query the database
    user_data <- dbGetQuery(pool, sprintf("SELECT user_id FROM users WHERE username = '%s' AND password_hash = '%s'",
                                          username, password_hash))
    
    if(nrow(user_data) == 1) {
      return(user_data$user_id[1]) # Return user ID on success
    } else {
      return(NULL) # Return NULL on failure
    }
    
  }
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = 1000, height = 800)
      result <- processedData()$result
      print(
        ggplot(result, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
          geom_point(alpha = 0.7, size = 2) +
          scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
          theme_minimal() +
          labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
      )
      dev.off()
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("heatmap_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = 1000, height = 800)
      significant_proteins <- processedData()$significant_proteins
      filtered_expression_matrix <- processedData()$normalized_matrix[significant_proteins$Protein, ]
      pheatmap(
        filtered_expression_matrix, 
        scale = "row",
        clustering_distance_rows = "euclidean", 
        clustering_method = "ward.D2",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1.5, 1.5, length.out = 101),
        legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
        legend_labels = c("-1.5 (Strongly Down)", "-1", "-0.5", "0 (No Change)", "0.5", "1", "1.5 (Strongly Up)"),
        main = "logFC Heatmap with Clustering"
      )
      dev.off()
    }
  )
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste("filtered_results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(tableData(), file, row.names = FALSE)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste("Proteomic_Report_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      # Save volcano plot to PNG
      volcano_path <- tempfile(fileext = ".png")
      png(volcano_path, width = 1000, height = 800)
      result <- processedData()$result
      print(
        ggplot(result, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
          geom_point(alpha = 0.7, size = 2) +
          scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
          theme_minimal() +
          labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")
      )
      dev.off()
      
      # Save heatmap plot to PNG
      heatmap_path <- tempfile(fileext = ".png")
      png(heatmap_path, width = 1000, height = 800)
      significant_proteins <- processedData()$significant_proteins
      filtered_expression_matrix <- processedData()$normalized_matrix[significant_proteins$Protein, ]
      pheatmap(
        filtered_expression_matrix, 
        scale = "row",
        clustering_distance_rows = "euclidean", 
        clustering_method = "ward.D2",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1.5, 1.5, length.out = 101),
        legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
        legend_labels = c("-1.5 (Strongly Down)", "-1", "-0.5", "0 (No Change)", "0.5", "1", "1.5 (Strongly Up)"),
        main = "logFC Heatmap with Clustering"
      )
      dev.off()
      
      # Render PDF using RMarkdown
      rmarkdown::render(
        input = "report_template.Rmd",
        output_file = file,
        params = list(
          volcano_path = volcano_path,
          heatmap_path = heatmap_path,
          result_table = processedData()$result
        ),
        envir = new.env(parent = globalenv())
      )
    }
  )
  
  
}

shinyApp(ui = ui, server = server)


