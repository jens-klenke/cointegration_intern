library(shinydashboard)
library(shiny)
library(ggplot2)
library(ggthemes)

#-----------------------------------------------------------------------------------------
# Sidebar
#-----------------------------------------------------------------------------------------
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dashboard", tabName = "Dashboard", icon = icon("dashboard")),
    menuItem("Data", tabName = "Data", icon = icon("chart-bar")),
    sliderInput("Lags", "Number of Lags included:", 
                min = 0, max = 20, 
                value = 1),
    radioButtons("Trend", "Deterministic components to be included:",
                                       choices = as.list(c("None", "Constant", "Trend")), 
                                       selected = "Constant"),
    radioButtons("Test", "Tests to aggregate:",
                 choices = as.list(c("Engle-Granger-Johansen", "all")),
                 selected = "all"))
)

#-----------------------------------------------------------------------------------------
# Body
#-----------------------------------------------------------------------------------------
body <- dashboardBody(
  tags$head(tags$style(HTML('
    .main-sidebar {background-color: #142029 !important;
    }
    .content-wrapper {background-image: linear-gradient(90deg, rgba(24,38,48,1) 0%, rgba(22,35,45,1) 35%) !important;
    }
    .box {background-color: #1B2B37 !important; border-top: #1B2B37 !important;
    }
    .box-title {color: #FFFFFF !important;
    }
    .control-label {color: #FFFFFF !important;
    }'                     
  ))),
  tabItems(
    tabItem(tabName = "Dashboard",
            fluidRow(
              box(tagList(tags$h1("bayerhanck()",
                                  style = "font-size: 29px; font-family: courier; font-weight: bold;
                                  text-align: start; color: #FFFFFF"),
                          tags$h2("A joint Test-Statistic for the Null of Non-Cointegration",
                                  style = "font-size: 18px; weight: 400; text-align: start;
                                  color: #FFFFFF")),
                  width = 12)),
            fluidRow(
              box(plotOutput("cdf_null", height = "250px"),
                  width = 8, height = "300px",
                  title = "CDF of the Null Distribution")),
            fluidRow(
              box(tableOutput("summary"), title = "Summary", width = 8, height = "300px"))),
    tabItem(tabName = "Data",
            fluidRow(
              box(title = "File input", status = "warning", width = 12,
                  fileInput(inputId = "csv_file", label = "Please add a csv file", accept = ".csv")),
              box(title = "Input", status = "warning", 
                  selectInput(inputId = "DepVar", label = "Dependent Variables", 
                              multiple = FALSE, choices = NULL,
                              selected = NULL),
                  selectInput(inputId = "IndVar", label = "Independent Variables", 
                              multiple = FALSE, choices = NULL)),
            )
    )
  )
)



ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  sidebar,
  body
)


server <- function(input, output, session) {
  output$cdf_null <- renderPlot({
    #inFile <- input$csv_file
    #if (is.null(inFile))
    #  return(NULL)
    #data <- read_csv(inFile$datapath)
    ggplot(data = iris) +
      geom_density(aes(x = Sepal.Length), col = "#18825B", linetype = 5, 
                   fill = "#18825B", alpha = 0.3) +
      theme(plot.background = element_rect(fill = "#1B2B37", colour = "#1B2B37"),
            panel.background = element_rect(fill = "#1B2B37"),
            panel.grid.major.y = element_line(size = 0.3, colour = "#546069"),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(size = 0.3, colour = "#546069"),
            panel.grid.minor.x = element_blank(),
            axis.text = element_text(colour = "#FFFFFF"),
            axis.line = element_line(size = 0.3, colour = "#546069"),
            axis.title = element_text(colour = "#FFFFFF"))
  })
  output$summary <- renderTable({
      inFile <- input$csv_file
      if (is.null(inFile))
          return(NULL)
      df <- readr::read_csv(inFile$datapath)
      bh <- bayerhanck(formula = get(input$DepVar) ~ get(input$IndVar),
                       data = df#, 
                       #lags = input$Lags,
                       #trend = input$Trend,
                       #test = input$Test
                       )
      sum_bh <- summary.bh.test(bh)[1]
      sum_bh <- tibble(sum_bh)
      sum_bh
  })
  observeEvent(input$csv_file, {
      inFile <- input$csv_file
      if (is.null(inFile))
          return(NULL)
      df <- readr::read_csv(inFile$datapath)
      updateSelectInput(session, "IndVar",
                        choices = as.list(names(select_if(df, is.numeric))))
  })
  observeEvent(input$csv_file, {
      inFile <- input$csv_file
      if (is.null(inFile))
          return(NULL)
      df <- readr::read_csv(inFile$datapath)
      updateSelectInput(session, "DepVar",
                        choices = as.list(names(select_if(df, is.numeric))))
  })
}


shinyApp(ui, server)


