# guides
# https://mastering-shiny.org/action-dynamic.html

library(shiny)
library(rgdal)
library(grid)
library(gridExtra)
library(bslib)
source("R/utils_Rdata.R")
options(shiny.r.port = 8050)
library(zip)

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
              # checkboxInput("quantile_50_test", "50% quantile", FALSE),
              # checkboxInput("quantile_80_test", "80% quantile", FALSE),
              checkboxInput("quantile_90_test", "90% quantile", FALSE),
              checkboxInput("quantile_95_test", "95% quantile", TRUE),
              checkboxInput("quantile_98_test", "98% quantile", FALSE),
              checkboxInput("quantile_99_test", "99% quantile", FALSE),
              checkboxInput("quantile_995_test", "99.5% quantile", FALSE),
            ), open = FALSE
          ),
          div(actionButton("submit_button_1", "Submit", class = "btn-success"), id = "div_submit_button"),
          div(downloadButton("download1"), id = "div_download_button"),
          div(
            style = "position: fixed; bottom: 2%; left: 2.5%; width: 4%; text-align: center;",
            a(
              href = "https://github.com/centre-for-humanities-computing/robusta_webapp",
              icon("github", "fa-2x"),
              target = "_blank",
              style = "color: #0B1215;"
            )
          ),
          div(
            style = "position: fixed; bottom: 2%; left: 8.5%; width: 4%; text-align: center;",
            a(
              href = "https://osf.io/u2gyq/",
              icon("paperclip", "fa-2x"),
              target = "_blank",
              style = "color: #0B1215;"
            )
          )
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
            numericInput("mc_simulations_2", label = "Number of CSR simulations", value = 25, min = 1, max = 10000),
            br(),
            numericInput("nbRobSc_2", label = "Number of Robustness Scenarios", value = 10, min = 1, max = 10000),
            br(),
            p("Quantile Selection"),
            accordion(
              accordion_panel("Expand for options",
                # checkboxInput("quantile_50_upload", "50% quantile", FALSE),
                # checkboxInput("quantile_80_upload", "80% quantile", FALSE),
                checkboxInput("quantile_90_upload", "90% quantile", FALSE),
                checkboxInput("quantile_95_upload", "95% quantile", TRUE),
                checkboxInput("quantile_98_upload", "98% quantile", FALSE),
                checkboxInput("quantile_99_upload", "99% quantile", FALSE),
                checkboxInput("quantile_995_upload", "99.5% quantile", FALSE),
              ), open = FALSE
            ),
            div(actionButton("submit_button_2", "Submit", class = "btn-success"), id = "div_submit_button"),
            div(downloadButton("download2"), id = "div_download_button"),
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
                   uiOutput("test"),
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

  # object to hold output from the function
  data_output <- reactiveValues(data = NULL)
  
  # Observer to start function processing
  observeEvent(input$submit_button_1, {

    file_shp <- readOGR("data/montecristi/mc-db-95-clean.shp")
    file_poly <- readOGR("data/montecristi/nmcpoly1.shp")
    
    quantiles <- c()
    if (input$quantile_90_test || input$quantile_90_upload){
      quantiles <- c(quantiles, 0.9)
    }
    if (input$quantile_95_test || input$quantile_95_upload){
      quantiles <- c(quantiles, 0.95)
    }
    if (input$quantile_98_test || input$quantile_98_upload){
      quantiles <- c(quantiles, 0.98)
    }
    if (input$quantile_99_test || input$quantile_99_upload){
      quantiles <- c(quantiles, 0.99)
    }
    if (input$quantile_995_test || input$quantile_995_upload){
      quantiles <- c(quantiles, 0.995)
    }
    
    data_output$data <- big_processing_func(file_shp = file_shp,
                        file_poly = file_poly,
                        nsim = input$mc_simulations_1/5, # it should be divided by number of clusters
                        clusters = 5,
                        nbRobSc = input$nbRobSc_1,
                        quantiles = quantiles,
                        quantile_50 = TRUE)
                        # quantile_50 = input$quantile_50_test || input$quantile_50_upload)

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
    
    quantiles <- c()
    if (input$quantile_90_test || input$quantile_90_upload){
      quantiles <- c(quantiles, 0.9)
    }
    if (input$quantile_95_test || input$quantile_95_upload){
      quantiles <- c(quantiles, 0.95)
    }
    if (input$quantile_98_test || input$quantile_98_upload){
      quantiles <- c(quantiles, 0.98)
    }
    if (input$quantile_99_test || input$quantile_99_upload){
      quantiles <- c(quantiles, 0.99)
    }
    if (input$quantile_995_test || input$quantile_995_upload){
      quantiles <- c(quantiles, 0.995)
    }

    data_output$data <- big_processing_func(file_shp = file_shp,
                         file_poly = file_poly,
                         nsim = input$mc_simulations_2/5, # it should be divided by number of clusters
                         clusters = 5,
                         nbRobSc = input$nbRobSc_2,
                         quantiles = quantiles,
                         quantile_50 = TRUE)

    function_ran(TRUE)
  })
    
output$plot_100_1 <- renderPlot({
      grid.draw(data_output$data[[1]])
    }, res = 96)
  output$plot_50_1 <- renderPlot({
      grid.draw(data_output$data[[2]][["50"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_50_2 <- renderPlot({
      grid.draw(data_output$data[[2]][["50"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_60_1 <- renderPlot({
      grid.draw(data_output$data[[2]][["60"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_60_2 <- renderPlot({
      grid.draw(data_output$data[[2]][["60"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_70_1 <- renderPlot({
      grid.draw(data_output$data[[2]][["70"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_70_2 <- renderPlot({
      grid.draw(data_output$data[[2]][["70"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_80_1 <- renderPlot({
      grid.draw(data_output$data[[2]][["80"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_80_2 <- renderPlot({
      grid.draw(data_output$data[[2]][["80"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_90_1 <- renderPlot({
      grid.draw(data_output$data[[2]][["90"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_90_2 <- renderPlot({
      grid.draw(data_output$data[[2]][["90"]][[sample(data_output$data[[4]], 1)]])
    }, res = 96)
  output$plot_compare <- renderPlot({
      grid.draw(data_output$data[[3]])
    }, res = 130)
  
  # Dynamically update the content of the tabPanel based on whether the function has run
  output$pcf_original_conditional <- renderUI({
    # Wait until the function has run
    if (function_ran()) {
      tagList(
        HTML("<br>"),
        h4("Original PCF"),
        HTML("<br>"),
        plotOutput("plot_100_1", width = "100%")
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
        tabsetPanel(
            tabPanel("50% of points",
            HTML("<br>"),
            plotOutput("plot_50_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_50_2", width = "100%")
            ),
            tabPanel("60% of points",
            HTML("<br>"),
            plotOutput("plot_60_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_60_2", width = "100%")
            ),
            tabPanel("70% of points",
            HTML("<br>"),
            plotOutput("plot_70_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_70_2", width = "100%")
            ),
            tabPanel("80% of points",
            HTML("<br>"),
            plotOutput("plot_80_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_80_2", width = "100%")
            ),
            tabPanel("90% of points",
            HTML("<br>"),
            plotOutput("plot_90_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_90_2", width = "100%")
            )
        )
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
        HTML("<br>"),
        plotOutput("plot_compare", width = "80%", height = "800px")
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
    filename = function() {
      paste("robusta_plots_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()
      
      # Initialize the list of plots to generate
      plot_list <- list(
        comparison_tool = data_output$data[[3]],
        original_pcf = data_output$data[[1]]
      )
      
      for (key in names(data_output$data[[2]])) {
        for (i in seq_along(data_output$data[[2]][[key]])) {
        plot_obj <- data_output$data[[2]][[key]][[i]]
        if (inherits(plot_obj, "grob")) {
          plot_list[[paste0(key, "_robustness_", i)]] <- plot_obj
        } else {
          warning(paste("Skipping non-plot object at key", key, "index", i))
        }
        }
      }
      
      # File paths for temporary plot files
      temp_files <- vector("character", length(plot_list))
      
      for (i in seq_along(plot_list)) {
        # Create temporary file paths
        temp_file <- file.path(temp_dir, paste0(names(plot_list)[i], ".png"))
        temp_files[i] <- temp_file
        
        # Save the ggplot as a PNG file
        png(temp_file, width = 1200, height = 600)
        grid.draw(plot_list[[i]])
        dev.off()
      }

      # Create the zip file
      zip::zipr(file, files = temp_files)
      
      # Cleanup temporary plot files after creating the zip
      on.exit(unlink(temp_files), add = TRUE)
    },
    contentType = "application/zip"
)


  # Create the download handler for the folder
  output$download2 <- downloadHandler(
    filename = function() {
      paste("robusta_plots_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()
      
      # Initialize the list of plots to generate
      plot_list <- list(
        comparison_tool = data_output$data[[3]],
        original_pcf = data_output$data[[1]]
      )
      
      for (key in names(data_output$data[[2]])) {
        for (i in seq_along(data_output$data[[2]][[key]])) {
        plot_obj <- data_output$data[[2]][[key]][[i]]
        if (inherits(plot_obj, "grob")) {
          plot_list[[paste0(key, "_robustness_", i)]] <- plot_obj
        } else {
          warning(paste("Skipping non-plot object at key", key, "index", i))
        }
        }
      }
      
      # File paths for temporary plot files
      temp_files <- vector("character", length(plot_list))
      
      for (i in seq_along(plot_list)) {
        # Create temporary file paths
        temp_file <- file.path(temp_dir, paste0(names(plot_list)[i], ".png"))
        temp_files[i] <- temp_file
        
        # Save the ggplot as a PNG file
        png(temp_file, width = 1200, height = 600)
        grid.draw(plot_list[[i]])
        dev.off()
      }

      # Create the zip file
      zip::zipr(file, files = temp_files)
      
      # Cleanup temporary plot files after creating the zip
      on.exit(unlink(temp_files), add = TRUE)
    },
    contentType = "application/zip"
)
  
}


shinyApp(ui, server)