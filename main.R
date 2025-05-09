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
library(readxl)
library(magick)
library(ggrepel)
library(aws.s3)
library(DT)

addResourcePath("assets", ".")

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
  dbExecute(
    pool,
    "
    CREATE TABLE IF NOT EXISTS users (
      user_id SERIAL PRIMARY KEY,
      username VARCHAR(50) UNIQUE NOT NULL,
      password_hash VARCHAR(255) NOT NULL,
      email VARCHAR(100) UNIQUE NOT NULL,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
  "
  )
  
  # Create analysis_history table
  dbExecute(
    pool,
    "
    CREATE TABLE IF NOT EXISTS analysis_history (
      analysis_id SERIAL PRIMARY KEY,
      user_id INTEGER REFERENCES users(user_id),
      filename VARCHAR(255) NOT NULL,
      s3_key VARCHAR(255) NOT NULL,
      ci_columns TEXT,
      he_columns TEXT,
      log2_threshold FLOAT,
      pval_threshold FLOAT,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
  "
  )
}

initializeDB(pool)

ui <- fluidPage(theme = shinytheme("cerulean"), uiOutput("dynamicUI"))

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
  
  user_session_data <- reactiveVal(NULL)
  
  output$status_message <- renderText({
    message <- statusMessage()
    paste0(message)
  })
  
  output$processingMessage <- renderText({
    message <- processingMessage()
    color <- if (message == "Analysis complete!")
      "green"
    else
      "red"
    paste0(message)
  })
  
  # Ensure UI renders HTML
  outputOptions(output, "status_message", suspendWhenHidden = FALSE)
  outputOptions(output, "processingMessage", suspendWhenHidden = FALSE)
  
  output$dynamicUI <- renderUI({
    if (!login_attempt()) {
      fluidPage(
        titlePanel(
          div("Proteomic Data Analysis Software", style = "text-align: center; font-size: 32px; font-weight: bold;")
        ),
        br(),
        div(
          actionButton("goToLogin", "Begin", style = "display: block; margin: 0 auto; font-size: 18px; padding: 10px 20px;"),
          style = "text-align: center;"
        ),
        br(),
        br(),
        div(
          h3("About Us", style = "text-align: center;"),
          p(
            "Our software is designed for quantitative proteomic data analysis, providing visualization and statistical analysis to assist researchers in understanding protein expression changes.",
            style = "text-align: center;"
          ),
          br(),
          h3("Supported Tools & Formats", style = "text-align: center;"),
          p(
            "We utilize cutting-edge technologies: R, Shiny, ggplot2, limma, and pheatmap.",
            style = "text-align: center;"
          ),
          p("Supported file formats: CSV, XLSX.", style = "text-align: center;"),
          br(),
          h3("Contributors", style = "text-align: center;"),
          p("Dr. John Doe, Dr. Jane Smith, Alex Johnson.", style = "text-align: center;")
        )
      )
    } else if (!user_auth()) {
      if (!show_register()) {
        fluidPage(titlePanel("Login"),
                  sidebarLayout(
                    sidebarPanel(
                      textInput("username", "Username"),
                      passwordInput("password", "Password"),
                      actionButton("login", "Login"),
                      textOutput("login_message"),
                      hr(),
                      actionButton("show_register_form", "Create New Account")
                    ),
                    mainPanel(h3(
                      "Please enter your credentials to continue."
                    ))
                  ))
      } else {
        fluidPage(
          titlePanel("Register New Account"),
          sidebarLayout(
            sidebarPanel(
              textInput("reg_username", "Username"),
              textInput("reg_email", "Email"),
              passwordInput("reg_password", "Password"),
              passwordInput("reg_confirm_password", "Confirm Password"),
              actionButton("register", "Register"),
              textOutput("register_message"),
              hr(),
              actionButton("back_to_login", "Back To Login")
            ),
            mainPanel(h3(
              "Create a new account to get started"
            ))
          )
        )
      }
    } else {
      navbarPage(  id = "main_nav",                      # give the navbarPage an id
                   title = "Protein Quantitation Analysis",
        tabPanel("Analysis", sidebarLayout(
          sidebarPanel(
            actionLink(
              inputId = "help_link",
              label   = "If you are unsure about how to use this application, navigate to the Usage page using this link."
            ),
            br(), br(),
            textOutput("status_message", inline = TRUE),
            br(),
            fileInput("file", "Upload CSV or Excel", accept = c(".csv", ".xlsx")),
            
            verbatimTextOutput("fileInfo"),
            textInput(
              "protein_columns",
              "Enter Protein Name columns (eg: 'Protein Name')",
              ""
            ),
            textInput(
              "gene_columns",
              "Enter Gene Name columns (eg: 'Gene Name')",
              ""
            ),
            textInput(
              "ci_columns",
              "Enter Patient Sample Columns (comma-separated):",
              ""
            ),
            textInput(
              "he_columns",
              "Enter Healthy Control Sample Columns (comma-separated):",
              ""
            ),
            textInput(
              "log2_threshold",
              "Enter Log2 Fold Change Threshold (eg: 1):",
              ""
            ),
            textInput(
              "pval_threshold",
              "Enter Adjusted P-Value Threshold (eg: 0.05):",
              ""
            ),
            actionButton("start_analysis", "Start Analysis"),
            textOutput("processingMessage", inline = TRUE),
            br(),
            br(),
            downloadButton("download_report", "Download Full Report (PDF)")
          )
          ,
          mainPanel(
            uiOutput("buttons"),
            
            br(),
            br(),
            
            # Volcano Plot Section
            conditionalPanel(
              condition = "input.toggle_volcano % 2 == 1",
              div(
                style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;",
                h4("Volcano Plot"),
                p(
                  "The volcano plot visualizes the relationship between the log fold change and adjusted p-value of proteins.
           Points in red indicate significantly upregulated proteins, while points in blue indicate significantly downregulated proteins.
           Gray points are not significant."
                ),
                plotOutput("volcanoPlot"),
                downloadButton("download_volcano", "Download Volcano Plot (PNG)")
                
                
              )
            ),
            
            # Heatmap Section
            conditionalPanel(
              condition = "input.toggle_heatmap % 2 == 1",
              div(
                style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;",
                h4("Heatmap Plot"),
                p(
                  "The heatmap represents the expression levels of significantly different proteins.
           Red indicates higher expression, blue indicates lower expression, and clustering shows relationships among proteins."
                ),
                plotOutput("heatmapPlot"),
                downloadButton("download_heatmap", "Download Heatmap Plot (PNG)")
                
              )
            ),
            
            # Filtered Table Section
            conditionalPanel(
              condition = "input.toggle_table % 2 == 1",
              div(
                style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;",
                h4("Filtered Results Table"),
                p(
                  "This table lists proteins that are significantly different between conditions based on the specified thresholds.
           It includes log fold change, adjusted p-values, and classification as upregulated, downregulated, or not significant."
                ),
                tableOutput("filteredResults"),
                actionButton("show_less", "Show Significant Proteins"),
                actionButton("show_more", "Show All Proteins"),
                
                br(),
                downloadButton("download_table_significant", "Download Table Significant Proteins (CSV)"),
                downloadButton("download_table", "Download Table For All Proteins (CSV)"),
                
                
              )
            ),
            
            conditionalPanel(
              condition = "input.fetch_pubmed % 2 == 1",
              div(
                style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;",
                h4("PubMed Article Links"),
                p(
                  "Below are links to PubMed articles related to studying the correlation between each significant protein and strokes."
                ),
                uiOutput("pubmedUI")
              )
            ),
            conditionalPanel(
              condition = "input.toggle_protein2 % 2 == 1",
              div(
                style = "border: 1px solid black; padding: 10px; margin-bottom: 20px;",
                h4("Protein2 Mapping"),
                p("This generates a protein interaction network image from STRING DB based on the significant proteins obtained from the analysis."),
                uiOutput("protein2Mapping"),
                downloadButton("download_protein2", "Download Protein Interaction Map (PNG)")
              )
            )
            
            
            
            
          )
          
        )),
        navbarMenu("Options", 
                   tabPanel(actionButton("history", "Analysis History")),
                   tabPanel(actionButton("logout", "Logout"))
                   ),
        tabPanel("About", fluidPage(
          titlePanel(
            div("Proteomic Data Analysis Software", style = "text-align: center; font-size: 32px; font-weight: bold;")
          ),
          br(),
          br(),
          div(
            h3("About Us", style = "text-align: center;"),
            p(
              "Our software is designed for quantitative proteomic data analysis, providing visualization and statistical analysis to assist researchers in understanding protein expression changes.",
              style = "text-align: center;"
            ),
            br(),
            h3("Supported Tools & Formats", style = "text-align: center;"),
            p(
              "We utilize cutting-edge technologies: R, Shiny, ggplot2, limma, and pheatmap.",
              style = "text-align: center;"
            ),
            p("Supported file formats: CSV, XLSX.", style = "text-align: center;"),
            br(),
            h3("Contributors", style = "text-align: center;"),
            p("Dr. John Doe, Dr. Jane Smith, Alex Johnson.", style = "text-align: center;")
          )
        )),
        tabPanel(title = "Usage", value = "usage",
                 fluidPage(
                   div(style = "text-align: center;",
                       
                       # Main heading
                       h2(style = "font-size:35px;","Usage"),
                       
                       # Intro instructions
                       p(
                         "To use this software, first upload your proteomics data on the Analysis tab.
   If you’d like to try a sample dataset, download it below and then re‑upload it here.",
                         style = "font-size:18px;"
                       ),
                       
                       
                       # Sample file download
                       h4("Sample CSV File"),
                       downloadButton("download_sample_csv", "Download ProteinquantitationHE&CI.csv"),
                       
                       br(), br(),
                       hr(style = "border-top: 1px solid #ccc; margin: 20px 0;"),
                       # Input fields explanation
                       h2(style = "font-size:35px;","Input Fields"),
                       tags$ul(style = "list-style-position: inside; padding-left: 0; text-align: center; font-size: 18px",
                               tags$li(strong("1. File Upload:"), "Select your proteomics dataset in CSV or XLSX format. This should contain your raw expression values."),
                               tags$li(strong("2. Protein Name Column:"), "Specify the header name for the column that lists protein identifiers (e.g. “Protein Name”)."),
                               tags$li(strong("3. Gene Name Column:"), "Enter the header for the column that contains gene symbols (e.g. “Gene Name”)."),
                               tags$li(strong("4. Patient Sample Columns:"), "Comma‑separated list of column names corresponding to patient (case) samples—these values will be compared to controls."),
                               tags$li(strong("5. Healthy Control Sample Columns:"), "Comma‑separated list of column names for your healthy (control) samples—used as the baseline for differential analysis."),
                               tags$li(strong("6. Log2 Fold Change Threshold:"), "A numeric cutoff (e.g. 1) for |log₂(fold change)|; proteins exceeding this change are considered biologically meaningful."),
                               tags$li(strong("7. Adjusted P‑Value Threshold:"), "A significance cutoff (e.g. 0.05) on the multiple‑testing corrected p‑values; proteins below this value are deemed statistically significant.")
                       ),
                       
                       br(),
                       hr(style = "border-top: 1px solid #ccc; margin: 20px 0;"),
                       # New section: How to use with sample file
                       h2(style = "font-size:35px;","How to use with the sample file"),
                       tags$ul(style = "list-style-position: inside; padding-left: 0; text-align: center;",
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   "1. Upload your sample file (or your own data) using the file upload box below:"
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/UPLOAD_BOX.png",
                                   alt   = "File Upload Box",
                                   style = "width:30%; height:50%; border:1px solid #333;"
                                 ),
                                 tags$br()
                               ),
                               tags$br(),
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   '2. Identify and enter the protein‐name column. In the sample CSV it’s labeled "Accession", so type "Accession":'
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/PROTEIN_COL.png",
                                   alt   = "Protein Column Input",
                                   style = "width:30%; height:50%; border:1px solid #333;"
                                 ),
                                 tags$br()
                               ),
                               tags$br(),
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   '3. Identify and enter the gene‐name column. In the sample CSV it’s labeled "Gene Name", so type "Gene Name":'
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/GENE_COL.png",
                                   alt   = "Gene Column Input",
                                   style = "width:30%; height:50%; border:1px solid #333;"
                                 ),
                                 tags$br()
                               ),
                               tags$br(),
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   '4. Specify your sample and control columns. In the sample file, patients are "CI-1,CI-2,CI-3,CI-4" and controls are "HE-1,HE-2,HE-3,HE-4":'
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/SAMPLE_COL.png",
                                   alt   = "Sample Column Input",
                                   style = "width:30%; height:50%; border:1px solid #333;"
                                 ),
                                 tags$br()
                               ),
                               tags$br(),
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   "5. Set your thresholds. For the sample run we recommend Log2 FC = 1 and Adjusted P‑value = 0.05:"
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/ML_VALUES.png",
                                   alt   = "Threshold Values Input",
                                   style = "width:30%; height:50%; border:1px solid #333;"
                                 ),
                                 tags$br()
                               ),
                               tags$br(),
                               tags$li(
                                 tags$div(
                                   style = "font-size:18px; font-weight:bold;",
                                   '6. Click **Start Analysis** and wait for the status to change from "Results not Ready" to "Ready" before exploring your outputs:'
                                 ),
                                 tags$br(),
                                 tags$img(
                                   src   = "assets/images/NOT_READY.png",
                                   alt   = "Not Ready Status",
                                   style = "max-width:45%; height:auto; margin-right:1em; border:1px solid #333;"
                                 ),
                                 tags$img(
                                   src   = "assets/images/READY.png",
                                   alt   = "Ready Status",
                                   style = "max-width:45%; height:auto; border:1px solid #333;"
                                 )
                               )
                       ),
                       
                       br(), br(),
                       hr(style = "border-top: 1px solid #ccc; margin: 20px 0;"),
                       
                       # Output fields explanation
                       # Output fields explanation with images and line breaks
                       h2(style = "font-size:35px;", "Output and Results"),
                       tags$ul(
                         style = "
    list-style-position: inside;
    padding-left: 0;
    text-align: center;
    font-size: 18px;
  ",
                         
                         # Toggle Button
                         tags$li(
                           strong("1. Output Toggles:"), "Can toggle each of the output panels on and off.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/TOGGLE_BUTTON.PNG",
                             alt   = "Toggle Buttons",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Volcano Plot
                         tags$li(
                           strong("2. Volcano Plot (Downloadable in PNG):"),
                           "Visualizes each protein’s log₂ fold change versus –log₁₀(adjusted p‑value). Up‑ and downregulated proteins are highlighted in red/blue",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/VOLCANO_PLOT.PNG",
                             alt   = "Volcano Plot",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Heatmap Plot
                         tags$li(
                           strong("3. Heatmap Plot (Downloadable in PNG):"),
                           "Shows row‑scaled expression of significant proteins across all samples, clustered to reveal patterns of co‑expression.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/HEATMAP_PLOT.PNG",
                             alt   = "Heatmap Plot",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Filtered Table
                         tags$li(
                           strong("4. Filtered Results Table (Downloadable in CSV):"),
                           "Tabular list of proteins passing both your log₂ fold change and adjusted p‑value thresholds, with their numerical statistics and regulation status. You can either download all the proteins collected or only the significant ones.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/FILTERED_TABLE.PNG",
                             alt   = "Filtered Results Table",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # PubMed Article Links
                         tags$li(
                           strong("5. PubMed Article Links:"),
                           "Clickable links to PubMed abstracts for each significant protein, aiding your literature exploration of protein–stroke associations. Clicking on a protein will open up a protein interaction map and brief information on the protein as shown in the next section.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/PUBMED_ARTICLES.PNG",
                             alt   = "PubMed Links",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Popup Modal / Protein Functions
                         tags$li(
                           strong("6. Protein Functions & Info:"),
                           "Network diagram from STRING‑DB showing known interactors for your clicked significant proteins. Includes additional information regarding protein function and catalytic activity.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/POPUP_MODAL.PNG",
                             alt   = "Popup Modal",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Protein Interaction Map
                         tags$li(
                           strong("7. Protein Interaction Map (PNG):"),
                           "Subnetwork focused solely on the significant proteins themselves, illustrating direct interactions within your filtered set.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/PROTEIN_INTERACTION.PNG",
                             alt   = "Protein Interaction Map",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         ),
                         tags$br(),
                         
                         # Download Full Report
                         tags$li(
                           strong("8. Download Full Report (Downloadable in PDF):"),
                           "The button allows you to create a summarized PDF report with all the plots, graphs and tables that you have collected as a result/output.",
                           tags$br(),
                           tags$br(),
                           tags$img(
                             src   = "assets/images/DOWNLOAD_REPORT.PNG",
                             alt   = "Download Report Button",
                             style = "width:30%; height:50%; border:1px solid #333;"
                           ),
                           tags$br()
                         )
                       )
                       
                       
                   )  # end div(style="text-align:center;")
                 )  # end fluidPage()
        )  # end tabPanel("Usage")
        
        
        
        
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
      style = "display: flex; gap: 10px; flex-wrap: wrap;",
      actionButton(
        "toggle_volcano",
        "Volcano Plot",
        style = paste(
          "background-color:", btn_color, " !important;",
          "border-color:", btn_color, " !important;",
          "color: black; width: 150px; height: 40px;"
        ),
        disabled = is_disabled
      ),
      actionButton(
        "toggle_heatmap",
        "Heatmap Plot",
        style = paste(
          "background-color:", btn_color, " !important;",
          "border-color:", btn_color, " !important;",
          "color: black; width: 150px; height: 40px;"
        ),
        disabled = is_disabled
      ),
      actionButton(
        "toggle_table",
        "Filtered Table",
        style = paste(
          "background-color:", btn_color, " !important;",
          "border-color:", btn_color, " !important;",
          "color: black; width: 150px; height: 40px;"
        ),
        disabled = is_disabled
      ),
      actionButton(
        "fetch_pubmed",
        "PubMed Articles",
        style = paste(
          "background-color:", btn_color, " !important;",
          "border-color:", btn_color, " !important;",
          "color: black; width: 150px; height: 40px;"
        ),
        disabled = is_disabled
      ),
      actionButton(
        "toggle_protein2",
        "Protein2 Mapping",
        style = paste(
          "background-color:", btn_color, " !important;",
          "border-color:", btn_color, " !important;",
          "color: black; width: 150px; height: 40px;"
        ),
        disabled = is_disabled
      )
    )
  })
  
  
  
  observeEvent(input$goToLogin, {
    login_attempt(TRUE)
  })
  
  # server.R (or inside server <- function(input, output, session) { … })
  
  output$download_sample_csv <- downloadHandler(
    filename = function() {
      "ProteinquantitationHE&CI.csv"
    },
    content = function(file) {
      # assumes ProteinquantitationHE&CI.csv sits in your app directory
      file.copy("ProteinquantitationHE&CI.csv", file)
    },
    contentType = "text/csv"
  )
  
  
  observeEvent(input$login, {
    user_result <- verifyLogin(input$username, input$password)
    if (!is.null(user_result)) {
      user_id <- user_result
      user_auth(TRUE)
      user_session_data(input$username)
      session$userData$user_id <- user_result
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
    showModal(modalDialog("Logged out..."))
    user_session_data(NULL)
    user_auth(FALSE)
    login_attempt(FALSE)
  })
  
  observeEvent(input$history, {
    
    username <- user_session_data()
    
    if (!is.null(username)) {
      query <- "SELECT a.analysis_id, u.username, a.filename, a.created_at, a.s3_key
                FROM analysis_history a
                JOIN users u ON a.user_id = u.user_id
                WHERE u.username = $1
                ORDER BY a.created_at DESC"
      
      # Execute query with the username
      result <- dbGetQuery(pool, query, params = list(username))
      
      if (nrow(result) > 0) {
        bucket_name <- "fyp-analysis-bucket"
        result$Download <- paste0(
          '<a href="https://', bucket_name, '.s3.amazonaws.com/', result$s3_key, '" target="_blank">Download Report</a>'
        )
        
        result_subset <- result[, c("filename", "created_at", "Download")]
        
        output$analysis_history_table <- renderDT({
          datatable(
            result_subset,
            escape = FALSE,
            options = list(pageLength = 10, lengthMenu = c(5, 10, 25, 50))
          )
        })
      } else {
        output$analysis_history_table <- renderDT({
          datatable(
            data.frame(Message = "No analysis history found."),
            options = list(dom = 't')  # just show table body
          )
        })
      }
      
      showModal(modalDialog(
        title = "Analysis History",
        DTOutput("analysis_history_table"),
        footer = modalButton("Close"),
        easyClose = TRUE,
        size = "l",
        class = "modal-lg-custom"
      ))
      
      
  }})
  
  
  
  dataset <- reactive({
    req(input$file, input$protein_columns)
    
    file_ext <- tools::file_ext(input$file$name)
    
    if (file_ext == "xlsx") {
      # Read the XLSX file
      temp_df <- read_excel(input$file$datapath)
      
      # Convert to data.frame (read_excel returns tibble)
      temp_df <- as.data.frame(temp_df, check.names = FALSE)
      
      # Find the index of the column that matches the user input
      col_num <- which(names(temp_df) == input$protein_columns)
      
      # Set rownames and return
      rownames(temp_df) <- temp_df[[col_num]]
      temp_df[[col_num]] <- NULL  # Optional: remove the rowname column
      return(temp_df)
    } else {
      # For CSV: find column index based on header row
      column_names <- names(read.csv(
        input$file$datapath,
        nrows = 1,
        check.names = FALSE
      ))
      col_num <- which(column_names == input$protein_columns)
      
      read.csv(input$file$datapath,
               row.names = col_num,
               check.names = FALSE)
    }
  })
  
  
  
  output$fileInfo <- renderText({
    req(input$file)
    file <- input$file
    file_ext <- tools::file_ext(file$name)
    
    if (file_ext == "xlsx") {
      raw_data <- read_excel(file$datapath)
      raw_data <- as.data.frame(raw_data, check.names = FALSE)
    } else {
      raw_data <- read.csv(file$datapath, check.names = FALSE)
    }
    
    paste(
      "File Name:",
      file$name,
      "\n",
      "File Type:",
      file$type,
      "\n",
      "File Size:",
      file$size,
      "bytes\n",
      "Row Count:",
      nrow(raw_data),
      "\n",
      "Column Names:",
      paste(colnames(raw_data), collapse = ", ")
    )
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
    req(
      input$file,
      input$protein_columns,
      # input$gene_columns,
      input$ci_columns,
      input$he_columns,
      input$log2_threshold,
      input$pval_threshold
    )
    log2_threshold <- as.numeric(input$log2_threshold)
    pval_threshold <- as.numeric(input$pval_threshold)
    data <- dataset()
    if (input$preprocess_ig && input$protein_desc_column != "") {
      desc_col <- input$protein_desc_column
      if (desc_col %in% colnames(data)) {
        data <- data[!grepl("immunoglobulin", data[[desc_col]], ignore.case = TRUE), ]
      } else {
        showNotification("Protein description column not found.", type = "error")
      }
    }
    ci_cols <- strsplit(input$ci_columns, ",")[[1]]
    he_cols <- strsplit(input$he_columns, ",")[[1]]
    ci_cols <- trimws(ci_cols)
    he_cols <- trimws(he_cols)
    
    if (!all(ci_cols %in% colnames(data)) ||
        !all(he_cols %in% colnames(data))) {
      showNotification("Invalid column names. Please check and try again.", type = "error")
      return(NULL)
    }
    
    CI_data <- data[, ci_cols, drop = FALSE]
    HE_data <- data[, he_cols, drop = FALSE]
    # Convert all columns to numeric
    CI_data[] <- lapply(CI_data, function(x) as.numeric(as.character(x)))
    HE_data[] <- lapply(HE_data, function(x) as.numeric(as.character(x)))
    
    
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
    design <- model.matrix( ~ 0 + group)
    contrast.matrix <- makeContrasts(groupCI - groupHE, levels = design)
    
    fit <- lmFit(normalized_matrix, design)
    fitlc <- contrasts.fit(fit, contrast.matrix)
    fitlc <- eBayes(fitlc)
    
    result <- topTable(fitlc, adjust = "BH", number = Inf)
    
    result$Protein <- rownames(result)
    gene_col <- input$gene_columns
    
    if (!is.null(gene_col) && gene_col != "" && gene_col %in% colnames(data)) {
      # Add Gene column to result
      result$Gene <- data[rownames(result), gene_col]
      
      # Reorder columns: Protein, Gene, rest
      result <- result[, c("Protein", "Gene", setdiff(colnames(result), c("Protein", "Gene")))]
    } else {
      # Just include Protein and rest
      result <- result[, c("Protein", setdiff(colnames(result), c("Protein")))]
    }
    
    result$Significance <- ifelse(
      result$adj.P.Val < pval_threshold &
        abs(result$logFC) > log2_threshold,
      ifelse(
        result$logFC > log2_threshold,
        "Upregulated",
        "Downregulated"
      ),
      "Not Significant"
    )
    
    significant_proteins <- subset(result, Significance != "Not Significant")
    statusMessage("Analysis complete!")
    processingMessage("Done!")
    tableData(head(result, 20))
    list(
      normalized_matrix = normalized_matrix,
      result = result,
      significant_proteins = significant_proteins,
      ci_columns = ci_cols,
      he_columns = he_cols,
      log2_threshold = log2_threshold,
      pval_threshold = pval_threshold
    )

  })
  
  
  
  output$volcanoPlot <- renderPlot({
    req(processedData())
    result <- processedData()$result
    
    ggplot(result, aes(
      x     = logFC,
      y     = -log10(adj.P.Val),
      color = Significance
    )) +
      geom_point(alpha = 0.7, size = 2) +
      # only label the significant proteins to avoid clutter:
      geom_text_repel(
        data    = subset(result, Significance != "Not Significant"),
        aes(label = Protein),
        size    = 3,            # text size
        max.overlaps = 10       # allow up to 10 overlaps before dropping labels
      ) +
      scale_color_manual(
        values = c(
          "Upregulated"   = "red",
          "Downregulated" = "blue",
          "Not Significant" = "gray"
        )
      ) +
      theme_minimal() +
      labs(
        title = "Volcano Plot",
        x     = "Log2 Fold Change",
        y     = "-Log10 Adjusted P‑Value"
      )
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
      breaks = seq(-1.5, 1.5, length.out = 101),
      # Ensure color mapping covers range
      legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
      # Custom legend values
      legend_labels = c(
        "-1.5 (Strongly Down)",
        "-1",
        "-0.5",
        "0 (No Change)",
        "0.5",
        "1",
        "1.5 (Strongly Up)"
      ),
      main = "logFC Heatmap with Clustering",
      fontsize_col = 8
    )
    
  })
  
  output$filteredResults <- renderTable({
    req(processedData())
    significant_results <- processedData()$result[processedData()$result$Significance != "Not Significant", ]
    significant_results
  })
  
  
  observeEvent(input$show_more, {
    tableData(processedData()$result)
  })
  
  observeEvent(input$show_less, {
    filtered_results <- processedData()$result[processedData()$result$Significance != "Not Significant", ]
    tableData(filtered_results)
  })
  
  
  getPubMedLinks <- function(protein_name) {
    query <- URLencode(paste0(protein_name, "protein stroke"))
    esearch_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
      "db=pubmed&retmode=json&retmax=5&term=",
      query
    )
    
    search_response <- tryCatch({
      GET(esearch_url)
    }, error = function(e)
      return(NULL))
    if (is.null(search_response) ||
        http_error(search_response))
      return(NULL)
    
    search_result <- content(search_response, as = "parsed", type = "application/json")
    ids <- search_result$esearchresult$idlist
    if (length(ids) == 0)
      return(NULL)
    
    # Fetch titles using esummary
    esummary_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
      "db=pubmed&retmode=json&id=",
      paste(ids, collapse = ",")
    )
    
    summary_response <- tryCatch({
      GET(esummary_url)
    }, error = function(e)
      return(NULL))
    if (is.null(summary_response) ||
        http_error(summary_response))
      return(NULL)
    
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
  
  observe({
    req(pubmedLinks())
    
    lapply(pubmedLinks(), function(entry) {
      acc <- entry$protein
      
      observeEvent(input[[paste0("show_string_", acc)]], {
        
        meta <- get_uniprot_metadata(acc)
        gene <- if (!is.null(meta$gene) && !is.na(meta$gene)) meta$gene else acc
        
        showModal(modalDialog(
          title = paste0(gene, "'s Protein2 Interaction Map"),
          
          # STRING DB image
          tags$img(
            src = paste0("https://string-db.org/api/highres_image/network?identifiers=", acc, "&species=9606"),
            style = "width: 100%; height: auto; max-width: 1200px;"
          ),
          
          br(),
          
          # FUNCTION section
          if (!is.null(meta$fnc) && length(meta$fnc) > 0)
            tagList(
              tags$h4("Function"),
              tags$ul(
                lapply(meta$fnc, function(f) tags$li(f))
              )
            )
          else NULL,
          
          # CATALYTIC section
          if (!is.null(meta$catalytic) && length(meta$catalytic) > 0)
            tagList(
              tags$h4("Catalytic Activity"),
              tags$ul(
                lapply(meta$catalytic, function(line) {
                  tags$li(line)
                })
              )
            )
          else NULL,
          
          # Modal style
          tags$style(HTML("
          .modal-lg {
            width: 90% !important;
            max-width: 1400px;
          }
        ")),
          
          easyClose = TRUE,
          size = "l"
        ))
      })
    })
  })
  
  
  
  
  output$pubmedUI <- renderUI({
    req(pubmedLinks())
    
    tagList(
      h4("PubMed Links for Stroke-Related Proteins"),
      
      lapply(pubmedLinks(), function(entry) {
        if (is.null(entry$links)) return(NULL)
        
        # Pull UniProt accession
        accession <- entry$protein
        
        # Get metadata
        meta <- get_uniprot_metadata(accession)
        
        # Construct UniProt link
        uniprot_url <- paste0("https://www.uniprot.org/uniprotkb/", accession, "/entry")
        
        # Clean the gene name to remove curly brackets and their content
        if (!is.null(meta$gene) && !is.na(meta$gene)) {
          clean_gene <- trimws(gsub("\\{[^}]*\\}", "", meta$gene))
        } else {
          clean_gene <- ""
        }
        
        tagList(
          tags$h5(
            tags$a(href = uniprot_url, target = "_blank", accession),
            if (!is.null(meta$gene) && !is.na(meta$gene)) paste0(" (Gene: ", meta$gene, ")"),
            if (!is.null(meta$protein) && !is.na(meta$protein)) paste0(" — ", meta$protein)
          ),
          
          # ⇩ Button to trigger the STRING modal using the cleaned gene name in the label
          actionButton(
            inputId = paste0("show_string_", accession),
            label = paste0(clean_gene, "'s Protein2 Interaction Map")
          ),
          br(), br(),
          
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
  
  
  
  
  
  
  get_uniprot_accession <- function(protein_name, organism_id = 9606) {
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0(protein_name, " AND organism_id:", organism_id)
    response <- tryCatch({
      GET(url = base_url, query = list(query = query, format = "json"))
    }, error = function(e) return(NULL))
    
    if (is.null(response) || http_error(response)) return(NA)
    
    content_text <- content(response, as = "text", encoding = "UTF-8")
    
    parsed <- tryCatch({
      fromJSON(content_text)
    }, error = function(e) return(NULL))
    
    if (!is.null(parsed) &&
        is.list(parsed) &&
        "results" %in% names(parsed) &&
        length(parsed$results) > 0 &&
        "primaryAccession" %in% names(parsed$results[[1]])) {
      return(parsed$results[[1]]$primaryAccession)
    }
    
    return(NA)
  }
  
  get_uniprot_metadata <- function(accession) {
    url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".txt")
    
    txt <- tryCatch({
      readLines(url, warn = FALSE)
    }, error = function(e) return(NULL))
    
    if (is.null(txt)) return(NULL)
    
    # Helper: remove evidence tags
    strip_evidence <- function(text) {
      gsub("\\s*\\{[^}]+\\}", "", text)
    }
    
    # Gene Name
    gene_line <- grep("^GN\\s+Name=", txt, value = TRUE)
    gene <- if (length(gene_line) > 0) {
      raw_gene <- sub(".*Name=([^;]+);.*", "\\1", gene_line[1])
      strip_evidence(raw_gene)  # remove any curly bracket evidence tags
    } else NA
    
    
    # Protein Name
    recname_line <- grep("^DE\\s+RecName: Full=", txt, value = TRUE)
    protein_name <- if (length(recname_line) > 0) {
      raw_protein <- sub(".*Full=([^;]+);.*", "\\1", recname_line[1])
      strip_evidence(raw_protein)  # remove any curly brackets and their contents
    } else NA
    
    
    # --- MULTIPLE FUNCTION BLOCKS ---
    fnc <- list()
    fnc_starts <- grep("^CC   -!- FUNCTION:", txt)
    
    for (idx in fnc_starts) {
      lines <- c()
      first_line <- sub("^CC   -!- FUNCTION:\\s*", "", txt[idx])
      lines <- c(lines, first_line)
      
      i <- idx + 1
      while (i <= length(txt) && grepl("^CC       ", txt[i])) {
        lines <- c(lines, sub("^CC       ", "", txt[i]))
        i <- i + 1
      }
      
      block <- paste(lines, collapse = " ")
      block <- strip_evidence(block)
      fnc <- append(fnc, list(block))
    }
    
    # --- MULTIPLE CATALYTIC BLOCKS ---
    catalytic <- list()
    cat_starts <- grep("^CC   -!- CATALYTIC ACTIVITY:", txt)
    
    for (idx in cat_starts) {
      lines <- c()
      first_line <- sub("^CC   -!- CATALYTIC ACTIVITY:\\s*", "", txt[idx])
      lines <- c(lines, first_line)
      
      i <- idx + 1
      while (i <= length(txt) && grepl("^CC       ", txt[i])) {
        lines <- c(lines, sub("^CC       ", "", txt[i]))
        i <- i + 1
      }
      
      block <- paste(lines, collapse = " ")
      block <- strip_evidence(block)
      catalytic <- append(catalytic, list(block))
    }
    
    return(list(
      accession = accession,
      gene = gene,
      protein = protein_name,
      fnc = fnc,               # Now a list
      catalytic = catalytic    # Also a list
    ))
  }
  
  
  
  
  
  
  observeEvent(input$register, {
    if (input$reg_password == input$reg_confirm_password) {
      success <- registerUser(input$reg_username,
                              input$reg_password,
                              input$reg_email)
      if (success) {
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
  registerUser <- function(username, password, email) {
    # if username exist
    user_exists <- dbGetQuery(
      pool,
      sprintf(
        "SELECT 1 FROM users WHERE username ='%s' OR email ='%s' LIMIT 1",
        username,
        email
      )
    )
    if (nrow((user_exists)) > 0) {
      return (FALSE)
    }
    
    password_hash <- digest::digest(password, algo = "sha256")
    
    user_data <- dbGetQuery(
      pool,
      sprintf(
        "SELECT user_id FROM users WHERE username ='%s' AND password_hash = '%s'",
        username,
        password_hash
      )
    )
    dbExecute(
      pool,
      sprintf(
        "INSERT INTO users (username, password_hash, email) VALUES ('%s', '%s', '%s')",
        username,
        password_hash,
        email
      )
    )
    return(TRUE)
  }
  
  # User login function
  verifyLogin <- function(username, password) {
    # Hash the provided password
    password_hash <- digest::digest(password, algo = "sha256")
    
    # Query the database
    user_data <- dbGetQuery(
      pool,
      sprintf(
        "SELECT user_id FROM users WHERE username = '%s' AND password_hash = '%s'",
        username,
        password_hash
      )
    )
    
    if (nrow(user_data) == 1) {
      return(user_data$user_id[1]) # Return user ID on success
    } else {
      return(NULL) # Return NULL on failure
    }
    
  }
  
  output$protein2Mapping <- renderUI({
    req(processedData())
    
    # Extract the significant proteins from the processed data
    sig_proteins <- processedData()$significant_proteins
    
    if (nrow(sig_proteins) > 0) {
      # Join the protein identifiers using "%0d"
      protein_ids <- paste(sig_proteins$Protein, collapse = "%0d")
      
      # Construct the STRING API URL
      string_url <- paste0("https://string-db.org/api/highres_image/network?identifiers=", protein_ids, "&species=9606")
      
      # Display the protein interaction network image and the legend image below it.
      tags$div(
        tags$img(src = string_url, style = "width: 100%; height: auto; max-width: 1200px;"),
        tags$br(), tags$br(),
        tags$img(
          src = "assets/images/PROTEIN_LEGENDS.PNG",
          style = "display: block; margin-left: auto; margin-right: auto; width: 50%; height: auto;"
        )
        
      )
    } else {
      p("No significant proteins found.")
    }
  })
  
  
  output$download_protein2 <- downloadHandler(
    filename = function() {
      paste("protein_interaction_map", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Ensure processed data is available
      req(processedData())
      
      # Extract the significant proteins
      sig_proteins <- processedData()$significant_proteins
      
      if (nrow(sig_proteins) > 0) {
        # Join the protein identifiers using "%0d"
        protein_ids <- paste(sig_proteins$Protein, collapse = "%0d")
        
        # Construct the STRING API URL
        string_url <- paste0("https://string-db.org/api/highres_image/network?identifiers=", protein_ids, "&species=9606")
        
        library(magick)
        
        # Read the network image from STRING
        network_img <- image_read(string_url)
        
        # Read the local legend image (ensure the file is in the working directory or accessible via correct path)
        legend_img <- image_read("images/PROTEIN_LEGENDS.PNG")
        
        # Get the width of the network image
        network_info <- image_info(network_img)
        network_width <- network_info$width
        
        # Resize the legend image to 50% of the network image's width
        target_width <- floor(network_width * 0.5)
        legend_img <- image_scale(legend_img, paste0(target_width, "x"))
        
        # Get resized legend dimensions
        legend_info <- image_info(legend_img)
        legend_height <- legend_info$height
        
        # Create a blank transparent canvas (using "none" for transparency)
        canvas <- image_blank(width = network_width, height = legend_height, color = "none")
        
        # Calculate the left offset to center the legend image on the canvas
        offset_left <- floor((network_width - target_width) / 2)
        offset <- paste0("+", offset_left, "+0")
        
        # Composite the legend image onto the transparent canvas using operator "over"
        legend_centered <- image_composite(canvas, legend_img, offset = offset, operator = "over")
        
        # Append the network image and the centered legend image vertically
        combined_img <- image_append(c(network_img, legend_centered), stack = TRUE)
        
        # Write the combined image as a PNG file, preserving transparency
        image_write(combined_img, path = file, format = "png")
      } else {
        # If no significant proteins are found, generate a blank image with a message
        png(file, width = 1000, height = 800)
        plot.new()
        text(0.5, 0.5, "No significant proteins found.", cex = 1.5)
        dev.off()
      }
    }
  )
  
  
  
  
  
  
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      png(file, width = 1000, height = 800)
      result <- processedData()$result
      
      p <- ggplot(result, aes(
        x     = logFC,
        y     = -log10(adj.P.Val),
        color = Significance
      )) +
        geom_point(alpha = 0.7, size = 2) +
        # label only the Up/Down‑regulated points
        geom_text_repel(
          data    = subset(result, Significance != "Not Significant"),
          aes(label = Protein),
          size    = 3,
          max.overlaps = 10
        ) +
        scale_color_manual(
          values = c(
            "Upregulated"     = "red",
            "Downregulated"   = "blue",
            "Not Significant" = "gray"
          )
        ) +
        theme_minimal() +
        labs(
          title = "Volcano Plot",
          x     = "Log2 Fold Change",
          y     = "-Log10 Adjusted P‑Value"
        )
      
      print(p)
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
        legend_labels = c(
          "-1.5 (Strongly Down)",
          "-1",
          "-0.5",
          "0 (No Change)",
          "0.5",
          "1",
          "1.5 (Strongly Up)"
        ),
        main = "logFC Heatmap with Clustering"
      )
      dev.off()
    }
  )
  
  output$download_table_significant <- downloadHandler(
    filename = function() {
      paste("filtered_results_significant", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      sig_results <- processedData()$result[processedData()$result$Significance != "Not Significant", ]
      write.csv(sig_results, file, row.names = FALSE)
    }
  )
  observeEvent(input$help_link, {
    updateNavbarPage(session, "main_nav", selected = "usage")
  })
  
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
      # ---- Volcano Plot PNG ----
      volcano_path <- tempfile(fileext = ".png")
      png(volcano_path, width = 1000, height = 800)
      result <- processedData()$result
      
      p_volcano <- ggplot(result, aes(
        x     = logFC,
        y     = -log10(adj.P.Val),
        color = Significance
      )) +
        geom_point(alpha = 0.7, size = 2) +
        # label significant points
        geom_text_repel(
          data    = subset(result, Significance != "Not Significant"),
          aes(label = Protein),
          size        = 3,
          max.overlaps = 10
        ) +
        scale_color_manual(
          values = c(
            "Upregulated"     = "red",
            "Downregulated"   = "blue",
            "Not Significant" = "gray"
          )
        ) +
        theme_minimal() +
        labs(
          title = "Volcano Plot",
          x     = "Log2 Fold Change",
          y     = "-Log10 Adjusted P‑Value"
        )
      
      print(p_volcano)
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
        legend_labels = c(
          "-1.5 (Strongly Down)",
          "-1",
          "-0.5",
          "0 (No Change)",
          "0.5",
          "1",
          "1.5 (Strongly Up)"
        ),
        main = "logFC Heatmap with Clustering"
      )
      dev.off()
      
      # Save Protein2 Mapping plot to PNG (combine network image & legend using magick)
      protein2_path <- tempfile(fileext = ".png")
      sig_proteins <- processedData()$significant_proteins
      if (nrow(sig_proteins) > 0) {
        # Join the protein identifiers using "%0d"
        protein_ids <- paste(sig_proteins$Protein, collapse = "%0d")
        
        # Construct the STRING API URL
        string_url <- paste0("https://string-db.org/api/highres_image/network?identifiers=", protein_ids, "&species=9606")
        
        library(magick)
        
        # Read the network image from STRING
        network_img <- image_read(string_url)
        
        # Read the local legend image (adjust path if needed)
        legend_img <- image_read("images/PROTEIN_LEGENDS.PNG")
        
        # Get the width of the network image
        network_info <- image_info(network_img)
        network_width <- network_info$width
        
        # Resize the legend image to 50% of the network image's width
        target_width <- floor(network_width * 0.5)
        legend_img <- image_scale(legend_img, paste0(target_width, "x"))
        
        # Get resized legend dimensions
        legend_info <- image_info(legend_img)
        legend_height <- legend_info$height
        
        # Create a blank transparent canvas (with no background color)
        canvas <- image_blank(width = network_width, height = legend_height, color = "none")
        
        # Calculate the left offset to center the legend image on the canvas
        offset_left <- floor((network_width - target_width) / 2)
        offset <- paste0("+", offset_left, "+0")
        
        # Composite the legend image onto the transparent canvas using operator "over"
        legend_centered <- image_composite(canvas, legend_img, offset = offset, operator = "over")
        
        # Append the network image and the centered legend image vertically
        combined_img <- image_append(c(network_img, legend_centered), stack = TRUE)
        
        # Scale the combined image up by 40% (i.e. 140% of its original size)
        combined_img <- image_scale(combined_img, "140%")
        
        # Write the combined image as a PNG file, preserving transparency
        image_write(combined_img, path = protein2_path, format = "png")
        
      } else {
        # If no significant proteins are found, generate a blank image with a message
        png(protein2_path, width = 2500, height = 2300)
        plot.new()
        text(0.5, 0.5, "No significant proteins found", cex = 1.5)
        dev.off()
      }
      
      
      # Render PDF using RMarkdown and include all three images in the report
      rmarkdown::render(
        input = "report_template.Rmd",
        output_file = file,
        params = list(
          volcano_path = volcano_path,
          heatmap_path = heatmap_path,
          protein2_path = protein2_path,
          result_table = processedData()$result[processedData()$result$Significance != "Not Significant", ],
          username = input$username  # Ensure this variable holds the current username if needed
        ),
        envir = new.env(parent = globalenv())
      )
      
      # Upload PDF to S3
      bucket_name <- "fyp-analysis-bucket"
      s3_key <- paste0("reports/", input$username, "_", Sys.Date(), ".pdf")
      
      user_id <- session$userData$user_id
      filename <- paste("Proteomic_Report_", Sys.Date(), ".pdf", sep = "")
      ci_columns <- processedData()$ci_columns
      he_columns <- processedData()$he_columns
      log2_threshold <- processedData()$log2_threshold
      pval_threshold <- processedData()$pval_threshold
      result <- dbGetQuery(
        pool,
        "
        INSERT INTO analysis_history (
          user_id, filename, s3_key,
          ci_columns, he_columns,
          log2_threshold, pval_threshold
        ) VALUES ($1, $2, $3, $4, $5, $6, $7)
        RETURNING analysis_id
        ",
        params = list(
          user_id,            # value for user_id
          filename,           # value for filename
          s3_key,        # value for s3_key
          paste(ci_columns, collapse = ","),  # convert ci_columns to a comma-separated string
          paste(he_columns, collapse = ","), # convert he_columns to a comma-separated string
          log2_threshold,     # value for log2_threshold
          pval_threshold     # value for pval_threshold
        )
      )
      analysis_id <- result$analysis_id
      s3_key <- paste0("reports/", input$username, "_", Sys.Date(), "_", analysis_id, ".pdf")
      
      put_object(
        file = file,
        object = s3_key,
        bucket = bucket_name
      )
      
      dbExecute(
        pool,
        "UPDATE analysis_history SET s3_key = $1 WHERE analysis_id = $2",
        params = list(s3_key, analysis_id)
      )
      
      output$last_s3_key <- renderText({
        s3_key
      })
      
    }
  )
  
}

shinyApp(ui = ui, server = server)
