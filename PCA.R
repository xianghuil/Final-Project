#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(RColorBrewer)
library(plotly)

dataset <- read.csv('gene_sample_pca.csv',header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, row.names = 1)
headerNames=colnames(dataset)

ui <- fluidPage(
  pageWithSidebar(
    
    headerPanel("PCA Explorer"),
    
    sidebarPanel(
      selectInput('x', 'X', c("None"=FALSE,headerNames),headerNames[11]),
      selectInput('y', 'Y', c("None"=FALSE,headerNames),headerNames[12]),
      selectInput('z', 'z', c("None"=FALSE,headerNames),headerNames[13])
    ),
    
    mainPanel(
      plotlyOutput('plot')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    plot_ly(dataset, x = ~get(input$x), y = ~get(input$y), z = ~get(input$z), color = dataset$icgc_donor_id, colors = c('#BF382A', '#0C4B8E')) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')))
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

