# guides
# https://mastering-shiny.org/action-dynamic.html
library(shiny)
library(sf)
library(grid)
library(gridExtra)
library(bslib)
library(zip)
library(shinyjs)
library(markdown)

source("R/utils_Rdata.R")

# loading the text for showing in the server
name_text <- paste(readLines("data/texts/name.txt"), collapse = "\n")

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  # loading css style sheet
  tags$head(includeCSS("www/style.css")),
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
  # ),
  useShinyjs(), # for toggling download button on and off
  sidebarLayout(fluid = F, 
    sidebarPanel(width = 3,
                 HTML("<br><br>"),
                 h3("Settings", style = "text-align: center"),
                 HTML("<br>"),
                 id = "sidebar_id",
     #  h3("Settings", style = "text-align: center; margin: 20px 0px 40px 0px"),
      navbarPage("", id = "side_navbar",
       tabPanel("Test",
        hr(style = "margin: 0px 20px 20px 0px"),
        fluidPage(
          p("Simulations with Montecristi Dataset"),
          numericInput("mc_simulations_1", label = HTML(
              'Number of Monte Carlo simulations
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Number of simulationes can be anything, but we recommend between 1000-5000. Beware, values below 25 and above 10.000 can cause the processing to fail.
                </div>
              </span>'), value = 25, min = 1, max = 10000),
          br(),
          numericInput("nbRobSc_1", label = HTML(
              'Number of Robustness Scenarios
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Number of robustness scenarios can be anything, but we recommend between 100-3000.
                </div>
              </span>'), value = 5, min = 1, max = 10000),
          br(),
          HTML(
              'Quantile Selection
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  The most comprehensive comparison is achieved by using all options, but the processing is faster with fewer options.
                </div>
              </span>'),
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
          HTML("<br><br>"),
          splitLayout(cellWidths = c("50%", "50%"), 
            div(actionButton("submit_button_1", "Submit", class = "btn-success"), id = "div_submit_button"),
            div(disabled(downloadButton("download1")), id = "div_download_button")),
          HTML("<br><br><br>"),
        )
        ),
        tabPanel("Upload",
          hr(style = "margin: 0px 20px 20px 0px"),
          fluidPage(
            fileInput("upload1", label = HTML(
              'Upload Spatial Points
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Upload a .shp <i>and</i> a .shx file. Files must be uploaded simultaneously, so hold CTRL down while selecting files. 
                </div>
              </span>'
            ), accept = c(".shp", ".shx"), multiple = TRUE),
            br(),
            fileInput("upload2", label = HTML(
              'Upload Spatial Polygons
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Upload a .shp <i>and</i> a .shx file. Files must be uploaded simultaneously, so hold CTRL down while selecting files. 
                </div>
              </span>'), accept = c(".shp", ".shx"), multiple = TRUE),
            br(),
            numericInput("mc_simulations_2", label = HTML(
              'Number of Monte Carlo simulations
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Number of simulationes can be anything, but we recommend between 1000-5000. Beware, values below 25 can cause the processing to fail.
                </div>
              </span>'), value = 25, min = 1, max = 10000),
            br(),
            numericInput("nbRobSc_2", label = HTML(
              'Number of Robustness Scenarios
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  Number of robustness scenarios can be anything, but we recommend between 100-3000.
                </div>
              </span>'), value = 10, min = 1, max = 10000),
            br(),
            HTML(
              'Quantile Selection
              <span class="info-icon">&#9432;
                <div class="hoverbox">
                  The most comprehensive comparison is achieved by using all options, but the processing is faster with fewer options.
                </div>
              </span>'),
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
            HTML("<br><br>"),
            splitLayout(cellWidths = c("50%", "50%"), 
              div(actionButton("submit_button_2", "Submit", class = "btn-success"), id = "div_submit_button"),
              div(disabled(downloadButton("download2")), id = "div_download_button")),
            HTML("<br><br><br>"),
            )
        ),
      ),
    ),
    mainPanel(width = 9, 
              HTML("<br><br>"),
              h2(name_text),
              HTML("<br>"),
      navbarPage(title = "", id = "main_tab",
        tabPanel("Home",
                 hr(style = "margin: 0px 20px 20px 0px"),
                 fluidPage(
                   uiOutput("test"),
                   HTML("<i> For the best experience on the web app, we recommend zooming out the browser window to 75% or 80%. </i>"),
                   HTML("<br><br>"),
                   h4("Introduction"),
                   includeMarkdown("data/texts/introduction.md"),
                   HTML("<br>"),
                   actionButton("toggle_how_to", "▼ How to use", class = "header-btn"),
                   div(id = "how_to_content", style = "display: none;",
                       includeMarkdown("data/texts/how_to.md")),
                   HTML("<br><br>"),
                   h4("Results"),
                   includeMarkdown("data/texts/results.md"),
                   HTML("<br>"),
                   h4("Creators and Contact"),
                   includeMarkdown("data/texts/creators.md"),
                  #  HTML("<br><br>"),
                   includeMarkdown("data/texts/contact.md"),
                   HTML("<br>"),
                   h4("Funding Information"),
                   includeMarkdown("data/texts/funding_information.md"),
                   HTML("<br><br>"),
                   splitLayout(cellWidths = c("33%", "33%", "34%"),
                    div(a(
                      href = "https://research-and-innovation.ec.europa.eu/funding/funding-opportunities/funding-programmes-and-open-calls/horizon-europe_en",
                      img(src = "images/Funded-by-the-European-Union.png", width = "70%"),
                      target = "_blank"
                      )
                    ),
                    div(a(
                      href = "https://marie-sklodowska-curie-actions.ec.europa.eu/",
                      img(src = "images/logo_marie-curie.jpg", width = "30%", style = "margin-left: 20%"),
                      target = "_blank"
                      )
                    ),
                    div(a(
                      href = "https://international.au.dk/",
                      img(src = "images/aarhus-university-au-3-logo.png", width = "70%"),
                      target = "_blank"
                      )
                    )
                  ),
                  #  HTML("<br><br>"),
                   id = "home_page")
        ),
        tabPanel("Output",
                 hr(style = "margin: 0px 0px 20px 0px"),
                 fluidPage(
                   HTML("<br>"),
                   tabsetPanel(
                     tabPanel("PCF based on 100% of sites",
                              uiOutput("pcf_original_conditional")
                              ), 
                     tabPanel("Robustness PCF",
                              uiOutput("pcf_robust_conditional")
                              ), 
                     tabPanel("Comparison Tools",
                              uiOutput("comparison_tools_conditional")
                              )
                     ),
                 id = "home_page"),
        ),
        tabPanel("Methods",
                 hr(style = "margin: 0px 0px 20px 0px"),
                 fluidPage(
                   HTML("<br>"),
                   h4("Framework and background"),
                   includeMarkdown("data/texts/method.md"),
                   id = "home_page"),
        ),
        tabPanel("More Info",
                 hr(style = "margin: 0px 0px 20px 0px"),
                 fluidPage(
                   HTML("<br>"),
                   h4("Original data"),
                   includeMarkdown("data/texts/original_data.md"),
                   h4("Original paper"),
                   includeMarkdown("data/texts/original_paper.md"),
                   h4("Code repository"),
                   includeMarkdown("data/texts/github.md"),
                   h4("Data privacy"),
                   includeMarkdown("data/texts/data_privacy.md"),
                 id = "home_page"),
        ), selected = "Home"
      ),  
    )
  ),
)

server <- function(input, output, session) {
  
  # Reactive value to track if the function has run
  function_ran <- reactiveVal(FALSE)
  function_error <- reactiveVal(FALSE)

  # object to hold output from the function
  data_output <- reactiveValues(data = NULL)
  
  observeEvent(input$toggle_how_to, {
    toggle("how_to_content")  # Shows/hides the how-to section
  })
  
  # Observer to start function processing
  observeEvent(input$submit_button_1, {
    shinyjs::disable("submit_button_1")
    shinyjs::disable("submit_button_2")
    shinyjs::disable("download1")
    shinyjs::disable("download2")
    
    function_error(FALSE)

    file_shp <- st_read("data/montecristi/mc-db-95-clean.shp")
    file_poly <- st_read("data/montecristi/nmcpoly1.shp")
    
    quantiles <- c()
    if (input$quantile_90_test){
      quantiles <- c(quantiles, 0.9)
    }
    if (input$quantile_95_test){
      quantiles <- c(quantiles, 0.95)
    }
    if (input$quantile_98_test){
      quantiles <- c(quantiles, 0.98)
    }
    if (input$quantile_99_test){
      quantiles <- c(quantiles, 0.99)
    }
    if (input$quantile_995_test){
      quantiles <- c(quantiles, 0.995)
    }
    
    data_output$data <- tryCatch({big_processing_func(file_shp = file_shp,
                        file_poly = file_poly,
                        nsim = as.integer(input$mc_simulations_1/5), # it should be divided by number of clusters
                        clusters = 5,
                        nbRobSc = input$nbRobSc_1,
                        quantiles = quantiles,
                        quantile_50 = TRUE)}, error = function(e) {return(NULL)})

    if(is.null(data_output$data)){
      function_error(TRUE)
    } else {
       function_ran(TRUE)
    }

    updateTabsetPanel(session, "main_tab",
      selected = "Output"
    )
    
    shinyjs::enable("submit_button_1")
    shinyjs::enable("submit_button_2")
    shinyjs::enable("download1")

    })
    
  

  # Observer to start function processing
  observeEvent(input$submit_button_2, {

    shinyjs::disable("submit_button_1")
    shinyjs::disable("submit_button_2")
    shinyjs::disable("download1")
    shinyjs::disable("download2")
    
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
    
    file_shp <- st_read(paste(tempdirname_1, shp_upload1$name[grep(pattern = "*.shp$", shp_upload1$name)], sep = "/"))
    file_poly <- st_read(paste(tempdirname_2, shp_upload2$name[grep(pattern = "*.shp$", shp_upload2$name)], sep = "/"))
    
    quantiles <- c()
    if (input$quantile_90_upload){
      quantiles <- c(quantiles, 0.9)
    }
    if (input$quantile_95_upload){
      quantiles <- c(quantiles, 0.95)
    }
    if (input$quantile_98_upload){
      quantiles <- c(quantiles, 0.98)
    }
    if (input$quantile_99_upload){
      quantiles <- c(quantiles, 0.99)
    }
    if (input$quantile_995_upload){
      quantiles <- c(quantiles, 0.995)
    }

    data_output$data <- big_processing_func(file_shp = file_shp,
                         file_poly = file_poly,
                         nsim = as.integer(input$mc_simulations_2/5), # it should be divided by number of clusters
                         clusters = 5,
                         nbRobSc = input$nbRobSc_2,
                         quantiles = quantiles,
                         quantile_50 = TRUE)

    if(is.null(data_output$data)){
      function_error(TRUE)
    } else {
       function_ran(TRUE)
    }

    updateTabsetPanel(session, "main_tab",
      selected = "Output"
    )

    shinyjs::enable("submit_button_1")
    shinyjs::enable("submit_button_2")
    shinyjs::enable("download2")
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
        h4("Pair Correlation Function based on 100% of sites"),
        HTML("<br>"),
        plotOutput("plot_100_1", width = "100%")
      )
    } else if (function_error()) {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "An error occured, and is likely due either a too low or too high number of Monte Carlo simulations. Try running it again with a number of simluations between 25 and 10.000"
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

  output$pcf_robust_conditional <- renderUI({
    # Wait until the function has run
    if (function_ran()) {
      tagList(
        HTML("<br>"),
        h4("Robustness PCF"),
        HTML("<br>"),
        p("Only two randonmly sampled robustness scenarios are shown here. The rest can be downloaded by clicking the download button."),
        tabsetPanel(
            tabPanel("90% of points",
            HTML("<br>"),
            plotOutput("plot_90_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_90_2", width = "100%"),
            HTML("<br><br>")
            ),
            tabPanel("80% of points",
            HTML("<br>"),
            plotOutput("plot_80_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_80_2", width = "100%"),
            HTML("<br><br>")
            ),
            tabPanel("70% of points",
            HTML("<br>"),
            plotOutput("plot_70_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_70_2", width = "100%"),
            HTML("<br><br>")
            ),
            tabPanel("60% of points",
            HTML("<br>"),
            plotOutput("plot_60_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_60_2", width = "100%"),
            HTML("<br><br>")
            ),
            tabPanel("50% of points",
            HTML("<br>"),
            plotOutput("plot_50_1", width = "100%"),
            HTML("<br><br>"),
            plotOutput("plot_50_2", width = "100%"),
            HTML("<br><br>")
            )
        )
      )
    } else if (function_error()) {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "An error occured, and is likely due to low number of Monte Carlo simulations. Try running it again or increase number of simulations to at least 25."
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
    } else if (function_error()) {
      div(
        class = "center-content",  # Wrapper to center content
        div(
          class = "error-message",  # Custom error styling
          "An error occured, and is likely due to low number of Monte Carlo simulations. Try running it again or increase number of simulations to at least 25."
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
  
  observe({
      # toggleState("download1", !is.null(data_output$data))
      toggleClass("download1", "btn-success", !is.null(data_output$data))
    })
  
  observe({
      # toggleState("download2", !is.null(data_output$data))
      toggleClass("download2", "btn-success", !is.null(data_output$data))
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