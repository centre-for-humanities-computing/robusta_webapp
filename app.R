# guides
# https://mastering-shiny.org/action-dynamic.html

library(shiny)
library(rgdal)
library(grid)
library(gridExtra)
library(bslib)
source("R/utils_Rdata.R")
options(shiny.port = 8050)

getwd()

# loading the text for showing in the server
name_text <- paste(readLines("data/texts/name.txt"), collapse = "\n")
introduction_text <- paste(readLines("data/texts/introduction.txt"), collapse = "\n")
method_text <- paste(readLines("data/texts/method.txt"), collapse = "\n")
results_text <- paste(readLines("data/texts/results.txt"), collapse = "\n")
creators_text <- paste(readLines("data/texts/creators.txt"), collapse = "\n")
contact_text <- paste(readLines("data/texts/contact.txt"), collapse = "\n")
funding_text <- paste(readLines("data/texts/funding_information.txt"), collapse = "\n")

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  # loading css style sheet
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  sidebarLayout(
    sidebarPanel(id = "sidebar_id",
      h3("Settings", style = "text-align: center; margin: 20px 0px 40px 0px"),
      navbarPage("", id = "side_navbar",
       tabPanel("Test",
        hr(style = "margin: 0px 20px 20px 0px"),
        fluidPage(
          p("Simulations with Montecristi Dataset"),
          numericInput("mc_simulations_1", label = "Number of Monte Carlo simulations", value = 25, min = 1, max = 10000),
          br(),
          numericInput("nbRobSc_1", label = "Number of Robustness Scenarios", value = 5, min = 1, max = 10000),
          br(),
          p("Quantile Selection"),
          accordion(
            accordion_panel("Expand for options",
              checkboxInput("quantile_50", "50% quantile", FALSE),
              checkboxInput("quantile_80", "80% quantile", FALSE),
              checkboxInput("quantile_90", "90% quantile", FALSE),
              checkboxInput("quantile_95", "95% quantile", TRUE),
              checkboxInput("quantile_98", "98% quantile", FALSE),
              checkboxInput("quantile_99", "99% quantile", FALSE),
            ), open = FALSE
          ),
          div(actionButton("submit_button_1", "Submit", class = "btn-success"), id = "div_submit_button"),
          div(downloadButton("download1"), id = "div_download_button"),
        )
        ),
        tabPanel("Upload",
          hr(style = "margin: 0px 20px 20px 0px"),
          fluidPage(
            fileInput("upload1", label = HTML(
              'Upload Spatial Points
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  As a minimum you need to upload a .shp and a .shx file.
                </div>
              </span>'
            ), accept = c(".shp", ".cpg", ".dbf", ".prj", ".qmd", ".shx"), multiple = TRUE),
            br(),
            fileInput("upload2", label = "Upload Spatial Polygons", accept = c(".shp", ".cpg", ".dbf", ".prj", ".sbn", ".sbx", ".shx"), multiple = TRUE),
            br(),
            # selectInput("dataset", label = "Dataset", choices = c("Random Sampling", "Weighted Sampling")),
            # br(),
            numericInput("mc_simulations_2", label = "Number of CSR simulations", value = 25, min = 1, max = 10000),
            br(),
            numericInput("nbRobSc_2", label = "Number of Robustness Scenarios", value = 10, min = 1, max = 10000),
            br(),
            p("Quantile Selection"),
            accordion(
              accordion_panel("Expand for options",
                              checkboxInput("quantile_50", "50% quantile", FALSE),
                              checkboxInput("quantile_80", "80% quantile", FALSE),
                              checkboxInput("quantile_90", "90% quantile", FALSE),
                              checkboxInput("quantile_95", "95% quantile", TRUE),
                              checkboxInput("quantile_98", "98% quantile", FALSE),
                              checkboxInput("quantile_99", "99% quantile", FALSE),
              ), open = FALSE
            ),
            div(actionButton("submit_button_2", "Submit", class = "btn-success"), id = "div_submit_button"),
            div(downloadButton("download1"), id = "div_download_button"),
          )
        ),
      ),
      width = 2,
      style="font-size:20px;"
    ),
    mainPanel(
      navbarPage("",
        tabPanel("Home",
                 hr(style = "margin: 0px 20px 20px 0px"),
                 fluidPage(
                   h3(name_text),
                   p(introduction_text),
                   h4("Methods"),
                   p(method_text),
                   h4("How-to"),
                   p("How-to: To make it work, a range of files have to be uploaded. These are the vector files and the spatial polygon files. Vector files should be one each of the types: .cpg, .dbf, .prj, .qmd, .shp, and .shx...
                     Furthermore the number of Monte Carlo simulations, number of iterations, and number of robustness scenarios should be specified. Recommended values is in the range of XXX. More will be better, but above XXX will not change much"),
                   p("IMPORTANT: The calculations will take some time. Increasing the size of the free parameters can significantly increase processing time. A progress bar will show, how far you are. Due to the large processing time, we recommend downloading the results, if you would like to return to it."),
                   h4("Results"),
                   p(results_text),
                   h4("Creators"),
                   p(creators_text),
                   h4("Contact"),
                   p(contact_text),
                   h4("Funding Information"),
                   p(funding_text)
                   )
        ),
        tabPanel("Output",
                 hr(style = "margin: 0px 0px 20px 0px"),
                 fluidPage(
                   h3(name_text),
                   HTML("<br>"),
                   tabsetPanel(
                     tabPanel("Original PCF",
                              uiOutput("pcf_original_conditional")
                              ), 
                     tabPanel("Robustness PCF",
                              uiOutput("pcf_robust_conditional")
                              ), 
                     tabPanel("Comparison Tools",
                              uiOutput("comparison_tools_conditional")
                              )
                     ),
                   #uiOutput("dynamic_tabs")
                   # Output to display the message and click count
                   # verbatimTextOutput("message"),
                   # verbatimTextOutput("temp_func_value"),
                   # verbatimTextOutput("function_success"),
                   
                   # h4("Comparison Tools"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_1", width = "800px", height = "1100px"), plotOutput("plot_2", width = "800px", height = "1100px")),
                   # h4("Robustness Scenarios"),
                   # h5("For 50% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_50", width = "800px", height = "1100px"), plotOutput("plot_4_50", width = "800px", height = "1100px")),
                   # h5("For 60% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_60", width = "800px", height = "1100px"), plotOutput("plot_4_60", width = "800px", height = "1100px")),
                   # h5("For 70% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_70", width = "800px", height = "1100px"), plotOutput("plot_4_70", width = "800px", height = "1100px")),
                   # h5("For 80% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_80", width = "800px", height = "1100px"), plotOutput("plot_4_80", width = "800px", height = "1100px")),
                   # h5("For 90% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_90", width = "800px", height = "1100px"), plotOutput("plot_4_90", width = "800px", height = "1100px")),
                   # h5("For 100% of points"),
                   # splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_100", width = "800px", height = "1100px"), plotOutput("plot_4_100", width = "800px", height = "1100px")),
                   
                 ),
        ),
        tabPanel("More Info",
                 hr(style = "margin: 0px 0px 20px 0px"),
                 fluidPage(
                   h3(name_text),
                   h4("How to interpret the results"),
                   p("something something..."),
                   h4("Suggestions for further exploration"),
                   p("References, links, resources, whatever. Could be papers, could be youtube"),
                   h4("Original data source???"),
                   p("don't know what this is"),
                   h4("Original paper stuff??"),
                   p("also not sure"),
                 ),
        ),
      ), width = 10
    )
  ),
)

server <- function(input, output, session) {
  
  # Reactive value to track if the function has run
  function_ran <- reactiveVal(FALSE)
  
  # Observer to start function processing
  observeEvent(input$submit_button_1, {
    
    file_shp <- readOGR("data/montecristi/mc-db-95-clean.shp")
    file_poly <- readOGR("data/montecristi/nmcpoly1.shp")
    
    quantiles <- c()
    # 50 is special
    # if (input$quantile_50){
    #   quantiles <- c(quantiles, 50)
    # }
    if (input$quantile_80){
      quantiles <- c(quantiles, 0.8)
    }
    if (input$quantile_90){
      quantiles <- c(quantiles, 0.9)
    }
    if (input$quantile_95){
      quantiles <- c(quantiles, 0.95)
    }
    if (input$quantile_98){
      quantiles <- c(quantiles, 0.98)
    }
    if (input$quantile_99){
      quantiles <- c(quantiles, 0.99)
    }
    
    process_output <- big_processing_func(file_shp = file_shp,
                        file_poly = file_poly,
                        nsim = input$mc_simulations_1/5, # it should be divided by number of clusters
                        clusters = 5,
                        nbRobSc = input$nbRobSc_1,
                        quantiles = quantiles,
                        quantile_50 = input$quantile_50)
    
    output$plot_1_1 <- renderPlot({
      grid.draw(process_output[[1]][[1]])
    }, res = 96)
    
    output$plot_1_2 <- renderPlot({
      grid.draw(process_output[[1]][[2]])
    }, res = 96)
    
    output$plot_3_100_1 <- renderPlot({
      grid.draw(process_output[[6]])
    }, res = 96)
    output$plot_3_100_1_map <- renderPlot({
      grid.draw(process_output[[3]][[6]])
    }, res = 96)
    
    output$plot_4_50_1 <- renderPlot({
      grid.draw(process_output[[2]][[1]])
    }, res = 96)
    output$plot_4_50_1_map <- renderPlot({
      grid.draw(process_output[[3]][[1]])
    }, res = 96)
    output$plot_4_50_2 <- renderPlot({
      grid.draw(process_output[[4]][[1]])
    }, res = 96)
    output$plot_4_50_2_map <- renderPlot({
      grid.draw(process_output[[5]][[1]])
    }, res = 96)
    output$plot_4_60_1 <- renderPlot({
      grid.draw(process_output[[2]][[2]])
    }, res = 96)
    output$plot_4_60_1_map <- renderPlot({
      grid.draw(process_output[[3]][[2]])
    }, res = 96)
    output$plot_4_60_2 <- renderPlot({
      grid.draw(process_output[[4]][[2]])
    }, res = 96)
    output$plot_4_60_2_map <- renderPlot({
      grid.draw(process_output[[5]][[2]])
    }, res = 96)
    output$plot_4_70_1 <- renderPlot({
      grid.draw(process_output[[2]][[3]])
    }, res = 96)
    output$plot_4_70_1_map <- renderPlot({
      grid.draw(process_output[[3]][[3]])
    }, res = 96)
    output$plot_4_70_2 <- renderPlot({
      grid.draw(process_output[[4]][[3]])
    }, res = 96)
    output$plot_4_70_2_map <- renderPlot({
      grid.draw(process_output[[5]][[3]])
    }, res = 96)
    output$plot_4_80_1 <- renderPlot({
      grid.draw(process_output[[2]][[4]])
    }, res = 96)
    output$plot_4_80_1_map <- renderPlot({
      grid.draw(process_output[[3]][[4]])
    }, res = 96)
    output$plot_4_80_2 <- renderPlot({
      grid.draw(process_output[[4]][[4]])
    }, res = 96)
    output$plot_4_80_2_map <- renderPlot({
      grid.draw(process_output[[5]][[4]])
    }, res = 96)
    output$plot_4_90_1 <- renderPlot({
      grid.draw(process_output[[2]][[5]])
    }, res = 96)
    output$plot_4_90_1_map <- renderPlot({
      grid.draw(process_output[[3]][[5]])
    }, res = 96)
    output$plot_4_90_2 <- renderPlot({
      grid.draw(process_output[[4]][[5]])
    }, res = 96)
    output$plot_4_90_2_map <- renderPlot({
      grid.draw(process_output[[5]][[5]])
    }, res = 96)
    
    function_ran(TRUE)
    })

  # Observer to start function processing
  observeEvent(input$submit_button_2, {
    
    shp_upload1 <- input$upload1
    shp_upload2 <- input$upload2
    
    tempdirname_1 <- dirname(shp_upload1$datapath[1])
    tempdirname_2 <- dirname(shp_upload2$datapath[1])
    
    # Rename files
    for (i in 1:nrow(shp_upload1)) {
      file.rename(
        shp_upload1$datapath[i],
        paste0(tempdirname_1, "/", shp_upload1$name[i])
      )
    }
    
    # Rename files
    for (i in 1:nrow(shp_upload2)) {
      file.rename(
        shp_upload2$datapath[i],
        paste0(tempdirname_2, "/", shp_upload2$name[i])
      )
    }
    
    file_shp <- readOGR(paste(tempdirname_1, shp_upload1$name[grep(pattern = "*.shp$", shp_upload1$name)], sep = "/"))
    file_poly <- readOGR(paste(tempdirname_2, shp_upload2$name[grep(pattern = "*.shp$", shp_upload2$name)], sep = "/"))
    
    process_output <- big_processing_func(file_shp = file_shp,
                         file_poly = file_poly,
                         nsim = input$mc_simulations_2,
                         n_iter = input$n_iter_2,
                         nbRobSc = input$nbRobSc_2)
    
    output$plot_1_1 <- renderPlot({
      grid.draw(process_output[[1]][[1]])
    }, res = 96)
    
    output$plot_1_2 <- renderPlot({
      grid.draw(process_output[[1]][[2]])
    }, res = 96)
    
    output$plot_3_100_1 <- renderPlot({
      grid.draw(process_output[[6]])
    }, res = 96)
    output$plot_3_100_1_map <- renderPlot({
      grid.draw(process_output[[3]][[6]])
    }, res = 96)
    
    output$plot_4_50_1 <- renderPlot({
      grid.draw(process_output[[2]][[1]])
    }, res = 96)
    output$plot_4_50_1_map <- renderPlot({
      grid.draw(process_output[[3]][[1]])
    }, res = 96)
    output$plot_4_50_2 <- renderPlot({
      grid.draw(process_output[[4]][[1]])
    }, res = 96)
    output$plot_4_50_2_map <- renderPlot({
      grid.draw(process_output[[5]][[1]])
    }, res = 96)
    output$plot_4_60_1 <- renderPlot({
      grid.draw(process_output[[2]][[2]])
    }, res = 96)
    output$plot_4_60_1_map <- renderPlot({
      grid.draw(process_output[[3]][[2]])
    }, res = 96)
    output$plot_4_60_2 <- renderPlot({
      grid.draw(process_output[[4]][[2]])
    }, res = 96)
    output$plot_4_60_2_map <- renderPlot({
      grid.draw(process_output[[5]][[2]])
    }, res = 96)
    output$plot_4_70_1 <- renderPlot({
      grid.draw(process_output[[2]][[3]])
    }, res = 96)
    output$plot_4_70_1_map <- renderPlot({
      grid.draw(process_output[[3]][[3]])
    }, res = 96)
    output$plot_4_70_2 <- renderPlot({
      grid.draw(process_output[[4]][[3]])
    }, res = 96)
    output$plot_4_70_2_map <- renderPlot({
      grid.draw(process_output[[5]][[3]])
    }, res = 96)
    output$plot_4_80_1 <- renderPlot({
      grid.draw(process_output[[2]][[4]])
    }, res = 96)
    output$plot_4_80_1_map <- renderPlot({
      grid.draw(process_output[[3]][[4]])
    }, res = 96)
    output$plot_4_80_2 <- renderPlot({
      grid.draw(process_output[[4]][[4]])
    }, res = 96)
    output$plot_4_80_2_map <- renderPlot({
      grid.draw(process_output[[5]][[4]])
    }, res = 96)
    output$plot_4_90_1 <- renderPlot({
      grid.draw(process_output[[2]][[5]])
    }, res = 96)
    output$plot_4_90_1_map <- renderPlot({
      grid.draw(process_output[[3]][[5]])
    }, res = 96)
    output$plot_4_90_2 <- renderPlot({
      grid.draw(process_output[[4]][[5]])
    }, res = 96)
    output$plot_4_90_2_map <- renderPlot({
      grid.draw(process_output[[5]][[5]])
    }, res = 96)
    
    function_ran(TRUE)
  })
  
  # Dynamically update the content of the tabPanel based on whether the function has run
  output$pcf_original_conditional <- renderUI({
    # Wait until the function has run
    if (function_ran()) {
      tagList(
        HTML("<br>"),
        h4("Original PCF"),
        HTML("<br>"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_3_100_1", width = "800px", height = "600px"), plotOutput("plot_3_100_1_map", width = "800px", height = "600px"))
      )
    } else {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "Run the function to see the content of this tab."
        )
      )
    }
    
  })

  output$pcf_robust_conditional <- renderUI({
    # Wait until the function has run
    if (function_ran()) {
      tagList(
        HTML("<br>"),
        h4("Robustness PCF"),
        HTML("<br>"),
        h5("For 50% of points"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_50_1", width = "800px", height = "600px"), plotOutput("plot_4_50_1_map", width = "800px", height = "600px")),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_50_2", width = "800px", height = "600px"), plotOutput("plot_4_50_2_map", width = "800px", height = "600px")),
        HTML("<br>"),
        h5("For 60% of points"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_60_1", width = "800px", height = "600px"), plotOutput("plot_4_60_1_map", width = "800px", height = "600px")),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_60_2", width = "800px", height = "600px"), plotOutput("plot_4_60_2_map", width = "800px", height = "600px")),
        HTML("<br>"),
        h5("For 70% of points"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_70_1", width = "800px", height = "600px"), plotOutput("plot_4_70_1_map", width = "800px", height = "600px")),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_70_2", width = "800px", height = "600px"), plotOutput("plot_4_70_2_map", width = "800px", height = "600px")),
        HTML("<br>"),
        h5("For 80% of points"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_80_1", width = "800px", height = "600px"), plotOutput("plot_4_80_1_map", width = "800px", height = "600px")),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_80_2", width = "800px", height = "600px"), plotOutput("plot_4_80_2_map", width = "800px", height = "600px")),
        HTML("<br>"),
        h5("For 90% of points"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_90_1", width = "800px", height = "600px"), plotOutput("plot_4_90_1_map", width = "800px", height = "600px")),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_4_90_2", width = "800px", height = "600px"), plotOutput("plot_4_90_2_map", width = "800px", height = "600px")),
      )
    } else {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "Run the function to see the content of this tab."
        )
      )
    }
    
  })
  
  output$comparison_tools_conditional <- renderUI({
    # Wait until the function has run
    if (function_ran()) {
      tagList(
        HTML("<br>"),
        h4("Comparison Tools"),
        # HTML("<br>"),
        # splitLayout(cellWidths = c("50%", "50%"), h5("Original PCF"), h5("Robustness PCF")),
        HTML("<br>"),
        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot_1_1", width = "800px", height = "1100px"), plotOutput("plot_1_2", width = "800px", height = "1100px")),
      )
    } else {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "Run the function to see the content of this tab."
        )
      )
    }
    
  })
  
  
  # Create the download handler for the folder
  output$download1 <- downloadHandler(
    
    # Specify the name of the file to download
    filename = function() {
      paste("folder_contents", Sys.Date(), ".zip", sep = "")
    },
    
    # Define the content that will be sent to the user
    content = function(file) {
      folder_path <- "temp_data"  # Change this to your folder path
      
      # Create a temporary zip file of the folder contents
      zip::zip(zipfile = file, files = list.files(folder_path, full.names = TRUE))
    }
  )
  
}


shinyApp(ui, server)




