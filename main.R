####################################
# Data Professor                   #
# http://youtube.com/dataprofessor #
# http://github.com/dataprofessor  #
####################################

# Modified from Winston Chang,
# https://shiny.rstudio.com/gallery/shiny-theme-selector.html

# Concepts about Reactive programming used by Shiny,
# https://shiny.rstudio.com/articles/reactivity-overview.html

# Load R packages
library(shiny)
library(shinythemes)


ui <- fluidPage(
  titlePanel("Text File Reader"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload a Text File", accept = c(".txt"))
    ),
    mainPanel(
      uiOutput("textOutput")  # Dynamic UI to show text lines
    )
  )
)

server <- function(input, output) {
  output$textOutput <- renderUI({
    req(input$file)  # Ensure a file is uploaded
    text_lines <- readLines(input$file$datapath, warn = FALSE)
    
    # Create a list of paragraph elements for each line
    tagList(lapply(text_lines, function(line) {
      p(line)  # Wrap each line in a <p> tag
    }))
  })
}

shinyApp(ui, server)
