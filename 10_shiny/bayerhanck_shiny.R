library(shinydashboard)
library(shiny)
library(ggthemes)
library(dplyr)
library(ggplot2)

#-----------------------------------------------------------------------------------------
# Sidebar
#-----------------------------------------------------------------------------------------
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dashboard", tabName = "Dashboard", icon = icon("dashboard")),
    menuItem("Data", tabName = "Data", icon = icon("chart-bar")),
    sliderInput("Lags", "Number of Lags included:", 
                min = 0, max = 20, 
                value = 1,
                step = 1),
    radioButtons("Trend", "Deterministic components to be included:",
                 choices = as.list(c("none", "const", "trend")), 
                 selected = "const"),
    radioButtons("Test", "Tests to aggregate:",
                 choices = as.list(c("eg-j", "all")),
                 selected = "all"),
    radioButtons("Critical", "Level for the critical value:",
                 choices = as.list(c("0.01", "0.05", "0.10")),
                 selected = "0.05")
  )
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
    }
    .shiny-text-output {background-color: #1B2B37 !important; 
    border: 1px solid #1B2B37 !important; color: #FFFFFF !important;
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
              box(verbatimTextOutput("bh_test"),
                  width = 8, height = "300px"))),
    tabItem(tabName = "Data",
            fluidRow(
              box(title = "File input", status = "warning", width = 12,
                  fileInput(inputId = "csv_file", label = "Please add a csv file", accept = ".csv")),
              box(title = "Input", status = "warning", 
                  selectInput(inputId = "DepVar", label = "Dependent Variables", 
                              multiple = FALSE, choices = NULL,
                              selected = NULL),
                  selectInput(inputId = "IndVar", label = "Independent Variables", 
                              multiple = TRUE, choices = NULL)),
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
    inFile <- input$csv_file
    if (is.null(inFile))
      return(NULL)
    dat <- readr::read_csv(inFile$datapath)
    invisible(capture.output(
        bh <- bayerhanck(reformulate(
            req(input$IndVar), req(input$DepVar)),
            data = dat, 
            lags = input$Lags,
            trend = input$Trend,
            test = input$Test,
            crit = as.numeric(input$Critical))
        ))
    plot(bh, "dark")
  })
  output$bh_test <- renderPrint({
      inFile <- input$csv_file
      if (is.null(inFile))
          return("Please upload file and select variables.")
      dat <- readr::read_csv(inFile$datapath)
      invisible(capture.output(
      bh <- bayerhanck(reformulate(
          req(input$IndVar), req(input$DepVar)),
                       data = dat, 
                       lags = input$Lags,
                       trend = input$Trend,
                       test = input$Test,
                       crit = as.numeric(input$Critical)
                       )))
      summary(bh)
  })
  observeEvent(input$csv_file, {
      inFile <- input$csv_file
      if (is.null(inFile))
          return(NULL)
      dat <- readr::read_csv(inFile$datapath)
      updateSelectInput(session, "IndVar",
                        choices = as.list(names(select_if(dat, is.numeric))))
  })
  observeEvent(input$csv_file, {
      inFile <- input$csv_file
      if (is.null(inFile))
          return(NULL)
      dat <- readr::read_csv(inFile$datapath)
      updateSelectInput(session, "DepVar",
                        choices = as.list(names(select_if(dat, is.numeric))))
  })
}


shinyApp(ui, server)


