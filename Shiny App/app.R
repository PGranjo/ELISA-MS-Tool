#Required packages.
packages <- c("shiny", "shinyjs", "DT", "markdown", "roxygen2","shinyWidgets","svglite")

#Packages to install.
to_install <- setdiff(packages, rownames(x = installed.packages()))




#Install packages.
if(length(to_install) > 0){
  
  install.packages(to_install, repos = "http://cran.us.r-project.org")
  
}

#Load packages.
library(shiny)
library(shinyjs)
library(DT)
library(markdown)
library(roxygen2)
library(shinyWidgets)
library(svglite)

#Load Functions
source("functions.R")

#Load Database
#R object which ontains the necessary data for the Shiny application apply the predictive coverage
load("database_beta_1407.RData")

#Define the UI.
ui <- shiny::navbarPage(title = shiny::tags$div(
                          #Title of the apllication and logo in navbar
                          "ELISA-MS Tool", shiny::tags$img(src = "alphalyse_logo_full.svg",
                                                                style = "position: absolute; height: 20px; right: 15px; top: 15px;")),
                        id = "main_navigation_page",
                        #Import the custom CSS style.
                        header = shiny::includeCSS(path = "www/ELISA_MS_CSS_style.css"),
                        #Content tabPanel 1 - Analysis
                        shiny::tabPanel(title = "Analysis",
                                        value = "analysis_tab_panel",
                                        #Enable shinyjs package functionality.
                                        shinyjs::useShinyjs(),
                                        #Add Alphalyse logo at browser tab.
                                        shiny::tags$head(shiny::tags$link(rel = "icon", type = "image/png", sizes = "12x12", href = "alphalyse_logo.png")),
                                        #Analysis tabPanel contains a sidebarLayout.
                                        shiny::sidebarLayout(
                                          #Siderbar panel content.
                                          sidebarPanel = shiny::sidebarPanel(width = 3,
                                                                             shiny::uiOutput(outputId = "sidebar_panel_ui"),
                                                                             uiOutput("sliders_ui")),
                                          #Main panel (body) content.
                                          mainPanel = shiny::mainPanel(width = 9,
                                                                       shiny::uiOutput(outputId = "main_panel_ui"))
                                        )
                        ),
                        #Content tabPanel 2 - Documentation
                        shiny::tabPanel(title = "Documentation",
                                        value = "documentation_tab_panel",
                                        shiny::fluidRow(
                                          shiny::column(width = 6, offset = 3, shiny::includeMarkdown(path = "ELISA-MS Tool Documentation.md"))
                                        )
                        )
)

#The server side of the tool.
server <- function(input, output, session) {
  
  #Increase maximum file upload size to 30MB.
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  ####################################################################################################################################################
  ################################################################ SIDEBAR PANEL UI ##################################################################
  ####################################################################################################################################################
  
  #UI for the sidebar panel of the Analysis tabPanel.
  output$sidebar_panel_ui <- shiny::renderUI({
    
    #To render more than one UI elements, they are wrapped within a div.
    shiny::tags$div(
      
      # Widget for selecting the type of input data (FASTA, Excel, List of Accession Numbers).
      shiny::radioButtons(inputId = "input_type",
                          label = shiny::tags$div(shiny::icon(name = "file-alt", class = "fa-solid", style = "margin-right: 5px;"),
                                                  "Input Type:",
                                                  style = "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;"),
                          choices = c("FASTA File" = "fasta", "Excel File" = "excel", "List of Accession Numbers" = "list"),
                          selected = "fasta",
                          inline = TRUE,
                          width = "100%"),
      shiny::tags$hr(),
      
      # Conditional UI elements for different input types.
      shiny::conditionalPanel(
        condition = "input.input_type == 'fasta'",
        shiny::fileInput(inputId = "fasta_file_input",
                         label = shiny::tags$div(shiny::icon(name = "file", class = "fa-solid", style = "margin-right: 5px;"),
                                                 "Upload FASTA File:",
                                                 style = "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;"),
                         multiple = FALSE,
                         accept = c(".fasta"),
                         width = "100%")
      ),
      shiny::conditionalPanel(
        condition = "input.input_type == 'excel'",
        shiny::fileInput(inputId = "excel_file_input",
                         label = shiny::tags$div(shiny::icon(name = "file-excel", class = "fa-solid", style = "margin-right: 5px;"),
                                                 "Upload Excel File:",
                                                 style = "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;"),
                         multiple = FALSE,
                         accept = c(".xlsx", ".xlsm"),
                         width = "100%")
      ),
      shiny::conditionalPanel(
        condition = "input.input_type == 'list'",
        shiny::textAreaInput(inputId = "accession_list_input",
                             label = shiny::tags$div(shiny::icon(name = "list", class = "fa-solid", style = "margin-right: 5px;"),
                                                     "Enter Accession Numbers (comma or line separated):",
                                                     style = "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;"),
                             placeholder = "Enter accession numbers separated by commas or newlines",
                             width = "100%",
                             height = "200px")
      ),
      shiny::tags$hr(),
      
      # Widget for selecting the organism type currently available in the database
      shiny::radioButtons(inputId = "organism_input",
                          label = shiny::tags$div(shiny::icon(name = "dna", class = "fa-solid", style = "margin-right: 5px;"),
                                                  "Organism:",
                                                  style = "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;"),
                          choices = c("Escherichia Coli" = "Escherichia Coli", 
                                      "CHO" = "CHO", 
                                      "Spodoptera frugiperda" = "Spodoptera frugiperda"),
                          selected = "Escherichia Coli",
                          inline = TRUE,
                          width = "100%"),
      shiny::tags$hr(),
      # Action button to start processing
      shiny::fluidRow(
        shiny::column(width = 4, offset = 4,
                      shiny::actionButton(inputId = "process_input", label = "RUN", width = "100%")
        )
      ),
      shiny::tags$hr()
    )
  })

  
  
  
  ####################################################################################################################################################
  ################################################################## MAIN BODY UI ####################################################################
  ####################################################################################################################################################
  
  #UI for the main panel in the Analysis tab.
  output$main_panel_ui <- shiny::renderUI(expr = {
  
    shiny::tabsetPanel(id = "body_tabset_panel",
                       selected = "dataframe_view",
                       type = "tabs",
                       
                       # Tab panel for viewing the data frame
                       shiny::tabPanel(title = "Coverage Table",
                                       value = "dataframe_view",
                                       div(class = "center", DT::dataTableOutput(outputId = "dataframe_view")),
                                       uiOutput("download_button_excel")),
                       
                       # Tab panel for visualizing data coverage plot 
                      shiny::tabPanel(title = "Coverage Plot",
                                       value = "upset_plot",
                                       div(class = "center", plotOutput(outputId = "upset_plot", height = "800px")),
                                       uiOutput("download_button_coverage")
                       ),
                      
                      # Tab panel for generating and displaying Venn Diagrams
                     shiny::tabPanel(
                       title = "Venn Diagrams",
                       value = "venn_diagram",
                       uiOutput("comparison_select"),
                       div(class = "center", plotOutput("venn_diagram", height = "800px", width = "800px")),
                       uiOutput("download_button_veen")),
                     
                     # Tab panel for plotting PI vs Mw data
                      shiny::tabPanel(title = "PI vs Mw",
                                       value = "pi_vs_mw_plot",
                                       shiny::sidebarLayout(
                                         sidebarPanel = shiny::sidebarPanel(
                                           selectInput("kit1", "Select Kit 1", choices = NULL),
                                           selectInput("kit2", "Select Kit 2 (optional)", choices = c("None"))
                                         ),
                                         mainPanel = shiny::mainPanel(
                                           div(class = "center", plotOutput(outputId = "pi_vs_mw_plot", height = "800px")),
                                           uiOutput("download_button_pi_mw")
                                           )
                                         )
                                       
                       )

    )
  })
  
  
  ####################################################################################################################################################
  ############################################################## REACTIVE DATA #######################################################################
  ####################################################################################################################################################
  
  
  # Read and process Excel data reactively.
  excel_data <- reactive({
    req(input$excel_file_input)
    data <- as.data.frame(read_xlsx(input$excel_file_input$datapath, sheet = "HCP", .name_repair = "minimal"))
    #Order the data in a decreasing order based on the ppm expression
    data <- data[order(data[,"ppm"], decreasing = TRUE), ] 
    return(data)
  })
  
  
  # Process input data based on selected input type (FASTA, Excel, List) and organism.
  reactive_data_frame <- reactive({
    req(input$process_input)
    input$process_input
    
    accession_info <- switch(input$input_type,
                             "fasta" = {
                               req(input$fasta_file_input)
                               input$fasta_file_input$datapath
                             },
                             "excel" = {
                               req(input$excel_file_input)
                               gsub("\\\\", "/", input$excel_file_input$datapath)
                             },
                             "list" = {
                               req(input$accession_list_input)
                               input$accession_list_input
                             },
                             NULL)
      req(accession_info)

      organism <- input$organism_input
      protein_data <- coverage_prot(database_beta, accession_info, organism)
      table <- ic_hcp_table(database_beta, protein_data, organism)

      
      if (nrow(table) == 0 || ncol(table) == 0) {
        data.frame("Message" = "No data available to display.")
      } else {
        return(table)
      }
    })
  
  
  # Calculate the total number of HCP entries.
  total_hcp <- reactive({
    req(input$excel_file_input)
    df <- reactive_data_frame()
    if (is.null(df) || nrow(df) == 0) {
      return(0)
    } else {
      return(nrow(df))
    }
  })
  
  
  # Filter the data based on user inputs. The user input must be on excel format in order to proceed with the filtration
  filtered_data <- reactive({
    df <- reactive_data_frame()
    req(!is.null(df) && nrow(df) > 0)
    
    if (input$input_type == "excel") {
      data <- excel_data()
      req(data)
      
      req(input$numeric_slider,input$percentage_slider)
      selected_rows <- if (input$numeric_slider > 0) {
        hcp_filtration_number(data, input$numeric_slider, "numeric")
      } else if (input$numeric_slider == 0) {
        hcp_filtration_number(data, 1, "numeric")
      } else {
        NULL
      }
      
      if (is.null(selected_rows) || length(selected_rows) == 0) {
        return(df[hcp_filtration_number(data, 1, "numeric"), , drop = FALSE])
      }
      req(!is.null(selected_rows) && length(selected_rows) > 0)
      return(df[rownames(df) %in% selected_rows, , drop = FALSE])
    } else {
      return(df)
    }
  })
  
  
  # Object used in the bellow reactive objects or events which can be a filtered data frame based on user input or unfiltered one
  final_data <- reactive({
    df <- if (input$input_type == "excel") {
      req(filtered_data())
      filtered_data()
    } else {
      req(reactive_data_frame())
      reactive_data_frame()
    }
    req(!is.null(df) && nrow(df) > 0)
    return(df)
  })
  
  # Creates a reactive set list to be used to create the Venn Diagram
  set_list <- reactive({
    req(final_data())
    create_set_list(final_data())
  })
  
  #Data Frame with additional information which is the PI and MW for the PI and Mw plot
  comprehensive_table <- reactive({
    req(final_data())
    retrieve_comprehensive_protein_info(final_data())
    
  })

#########################################################################################################################################
###############################################Widgets Filtering HCP Excel part##########################################################
#########################################################################################################################################
  
  
  # Create the numeric input and sliders for the data frame filtration
  output$sliders_ui <- renderUI({
    req(input$input_type == "excel")
    req(total_hcp(), excel_data())
    total <- total_hcp()
    max_ppm <- ceiling(max(excel_data()[, "ppm"], na.rm = TRUE))
    req(total > 0)
    tagList(
      numericInput("numeric_input", "Top Highest Expressing HCP (Number)", value = total, min = 0, max = total, step = total / 5),
      sliderInput("numeric_slider", "", min = 0, max = total, step = total / 5, value = total),
      numericInput("percentage_input", "Top Highest Expressing HCP (%)", value = 100, min = 0, max = 100, step = 1),
      sliderInput("percentage_slider", "", min = 0, max = 100, value = 100, step = 1),
      numericInput("ppm_input", "Lowest IC HCP expression ppm", value = 0, min = 0, max = max_ppm, step = max_ppm / 5),
      sliderInput("ppm_slider", "", min = 0, max = max_ppm, value = 0, step = max_ppm / 5)
    )
  })


  original_values <- reactiveValues(
    numeric_slider = NULL,
    numeric_input = NULL,
    percentage_slider = NULL,
    percentage_input = NULL,
    ppm_slider = NULL,
    ppm_input = NULL
  )
  
  # Reactive values to manage update states
  updating <- reactiveValues(
    numeric_slider = FALSE,
    numeric_input = FALSE,
    percentage_slider = FALSE,
    percentage_input = FALSE,
    ppm_slider = FALSE,
    ppm_input = FALSE
  )
  
  
  observe({
    req(total_hcp(), excel_data())
    if (is.null(original_values$numeric_slider)) {
      total <- total_hcp()
      max_ppm <- ceiling(max(excel_data()$ppm, na.rm = TRUE))
      
      original_values$numeric_slider <- total
      original_values$numeric_input <- total
      original_values$percentage_slider <- 100
      original_values$percentage_input <- 100
      original_values$ppm_slider <- 0
      original_values$ppm_input <- 0
    }
  })

  
  # Function to update inputs based on changes
  update_inputs <- function(source) {
    req(total_hcp() > 0, source)
    data <- excel_data()
    total <- total_hcp()
    
    if (source == "numeric_slider" || source == "numeric_input") {
      new_numeric <- if (source == "numeric_slider") input$numeric_slider else input$numeric_input
      if (is.null(new_numeric) || is.na(new_numeric)) {
        new_numeric <- total
      } else {
        new_percentage <- round((new_numeric / total) * 100,5)
      }
    } else if (source == "percentage_slider" || source == "percentage_input") {
      new_percentage <- if (source == "percentage_slider") input$percentage_slider else input$percentage_input
      if (is.null(new_percentage) || is.na(new_percentage)) {
        new_percentage <- 100
      } else{
        new_numeric <- round((new_percentage / 100) * total,0)
      }
    } else if (source == "ppm_input" || source == "ppm_slider") {
      new_ppm <- if (source == "ppm_input") input$ppm_input else input$ppm_slider

      if (is.null(new_ppm) || is.na(new_ppm) || new_ppm < 0) {
        new_ppm <- 0
      }
      new_numeric <- nrow(data[data$ppm >= new_ppm, ])
      new_percentage <- round((new_numeric / total) * 100,5)
    }
    
    if (!is.na(new_percentage) && new_percentage < 1) {
      new_percentage <- 1
    }
    if (!is.na(new_numeric) && new_numeric < 1) {
      new_percentage <- 1
    }

    
    isolate({
      if (!is.na(new_numeric) && new_numeric > 0 && new_numeric <= total) {
        if (source != "numeric_slider") updateSliderInput(session, "numeric_slider", value = new_numeric)
        if (source != "numeric_input") updateNumericInput(session, "numeric_input", value = new_numeric)
        original_values$numeric_slider <- new_numeric
        original_values$numeric_input <- new_numeric
      }
      if (!is.na(new_percentage) && new_percentage > 0 && new_percentage <= 100) {
        if (source != "percentage_slider") updateSliderInput(session, "percentage_slider", value = new_percentage)
        if (source != "percentage_input") updateNumericInput(session, "percentage_input", value = new_percentage)
        original_values$percentage_slider <- new_percentage
        original_values$percentage_input <- new_percentage
      }
      
      if (source != "ppm_input" && source != "ppm_slider" && !is.null(data) && nrow(data) >= new_numeric) {
        new_ppm <- data[order(data$ppm, decreasing = TRUE), ][new_numeric, "ppm"]
        if (!is.na(new_ppm)) {
          if (source != "ppm_slider") updateSliderInput(session, "ppm_slider", value = new_ppm)
          if (source != "ppm_input") updateNumericInput(session, "ppm_input", value = new_ppm)
          original_values$ppm_slider <- new_ppm
          original_values$ppm_input <- new_ppm
        }
      }
    })
  }
  
  
  # Observe the numeric_slider input
  observeEvent(input$numeric_slider, {
    req(total_hcp() > 0,input$numeric_slider,)
    if (!updating$numeric_slider && original_values$numeric_slider != input$numeric_slider) {
      updating$numeric_slider <- TRUE
      update_inputs("numeric_slider")
      updating$numeric_slider <- FALSE
    }
  })
  
  # Observe the numeric_input input
  observeEvent(input$numeric_input, {
    req(total_hcp() > 0,input$numeric_input)
    if (!updating$numeric_input && original_values$numeric_input != input$numeric_input) {
      updating$numeric_input <- TRUE
      update_inputs("numeric_input")
      updating$numeric_input <- FALSE
    }
  })
  
  # Observe the percentage_slider input
  observeEvent(input$percentage_slider, {
    req(total_hcp() > 0,input$percentage_slider)
    if (!updating$percentage_slider && original_values$percentage_slider != input$percentage_slider) {
      updating$percentage_slider <- TRUE
      update_inputs("percentage_slider")
      updating$percentage_slider <- FALSE
    }
  })
  
  # Observe the percentage_input input
  observeEvent(input$percentage_input, {
    req(total_hcp() > 0,input$percentage_input)
    if (!updating$percentage_input && original_values$percentage_input != input$percentage_input) {
      updating$percentage_input <- TRUE
      update_inputs("percentage_input")
      updating$percentage_input <- FALSE
    }
  })
  
  # Observe the ppm_input input
  observeEvent(input$ppm_input, {
    req(total_hcp() > 0, input$ppm_input)
    if (!updating$ppm_input && original_values$ppm_input != input$ppm_input) {
      updating$ppm_input <- TRUE
      update_inputs("ppm_input")
      updating$ppm_input <- FALSE
    }
  })
  
  # Observe the ppm_slider input
  observeEvent(input$ppm_slider, {
    req(total_hcp() > 0, input$ppm_slider)
    if (!updating$ppm_slider && original_values$ppm_slider != input$ppm_slider) {
      updating$ppm_slider <- TRUE
      update_inputs("ppm_slider")
      updating$ppm_slider <- FALSE
    }
  })
  
  
  #################################################################################################################################################
  ############################################################## Buttons Apperances ###############################################################
  #################################################################################################################################################
  
  #Rendering of the download button from the Coverage Table tab
  output$download_button_excel <- renderUI({
    req(final_data(), nrow(final_data()) > 0 )
    shiny::downloadButton(outputId = "download_combined_excel", label = "Download Table")
  })
  
  #Download Button for the Upset Plot giving the possibility to customize Width, Height, and choose which saving format the user wants
  output$download_button_coverage <- renderUI({
    req(final_data(), nrow(final_data()) > 0 )
  fluidRow(
    column(2, numericInput("coverage_width", "Width:", value = 10, min = 5, max = 30, width = "50%")),
    column(2, numericInput("coverage_height", "Height:", value = 8, min = 5, max = 20, width = "50%")),
    column(4, shinyWidgets::radioGroupButtons(
      inputId = "coverage_format_selection",
      label = "Select Format:",
      choices = c("png", "pdf", "svg"),
      justified = TRUE,
      checkIcon = list(
        yes = icon("check-square")
      ),
      width = "100%"
    )),
    column(2, downloadButton("download_upset_plot", label = "Download", class = "btn-download download-button", width = "100%"))
  )})
  
  #Download Button for the Venn diagram giving the possibility to customize Width, Height, and choose which saving format the user wants
  output$download_button_veen <- renderUI({
    req(final_data(), nrow(final_data()) > 0, set_list())
    fluidRow(
      column(2, numericInput("venn_width", "Width:", value = 10, min = 5, max = 30, width = "50%")),
      column(2, numericInput("venn_height", "Height:", value = 8, min = 5, max = 20, width = "50%")),
      column(4, shinyWidgets::radioGroupButtons(
        inputId = "veen_format_selection",
        label = "Select Format:",
        choices = c("png", "pdf", "svg"),
        justified = TRUE,
        checkIcon = list(
          yes = icon("check-square")
        ),
        width = "100%"
      )),
      column(2, downloadButton("venn_diagram_download", label = "Download", class = "btn-download download-button", width = "100%"))
    )})
  
  #Download Button for the PI Vs Mw Plot giving the possibility to customize Width, Height, and choose which saving format the user wants
  output$download_button_pi_mw <- renderUI({
    req(comprehensive_table())
    fluidRow(
      column(2, numericInput("pi_mw_width", "Width:", value = 10, min = 5, max = 30, width = "50%")),
      column(2, numericInput("pi_mw_height", "Height:", value = 8, min = 5, max = 20, width = "50%")),
      column(4, shinyWidgets::radioGroupButtons(
        inputId = "pi_mw_format_selection",
        label = "Select Format:",
        choices = c("png", "pdf", "svg"),
        justified = TRUE,
        checkIcon = list(
          yes = icon("check-square")
        ),
        width = "100%"
      )),
      column(2, downloadButton("download_pi_mw", label = "Download", class = "btn-download download-button", width = "100%"))
    )
  })
    
  #Widget for the use to choose the kits to compare from
  output$comparison_select <- renderUI({
    req(set_list())
    selectInput("comparisons", "Select comparisons:",
                choices = names(set_list()),
                selected = names(set_list())[1:min(2, length(names(set_list())))],
                multiple = TRUE)
  })
  

  # Render the Coverage Table
   output$dataframe_view <- DT::renderDataTable({
    req(final_data())
    df <- final_data()
    if (nrow(df) == 0) {
      return(DT::datatable(data.frame("Message" = "No data available to display.")))
    } else {
    DT::datatable(data = df, options = list(scrollX = TRUE, pageLength = 20), class = "nowrap display")
  }})
   
  # Render the Venn Diagram
  output$venn_diagram <- renderPlot({
    req(input$comparisons)
    validate(
      need(length(input$comparisons) >= 2, "Please select at least two comparisons"),
      need(length(input$comparisons) <= 5, "Please select no more than five comparisons")
    )
    req(set_list())
    create_venn_diagram(set_list(), input$comparisons)
  })
  
  # Render the UpSet Plot
  output$upset_plot <- renderPlot({
    req(final_data())
    req(set_list())
    
    # Calculate input percentages
    input_percent <- calculate_input_percent(set_list(), final_data())
    
    # Generate and render the UpSet plot
    upset <- create_upset_diagram(input_percent)
    grid.draw(upset)
  })
  


  # Reactive expression to get column names (kit choices)
  kit_choices <- reactive({
    df <- final_data()
    if (is.null(df) || nrow(df) == 0) {
      return(character(0)) 
    } else {
      return(colnames(df))  # Return the column names as a character vector
    }
  })
  
  # Reactive expression for available kit2 choices
  available_kit2_choices <- reactive({
    req(input$kit1)
    df <- final_data()
    if (is.null(df) || nrow(df) == 0) {
      return(character(0)) 
    } else {
      return(setdiff(colnames(df), input$kit1))  # Return all columns except the selected kit1
    }
  })
  
  # Reactive expression for available kit1 choices
  available_kit1_choices <- reactive({
    req(input$kit2)
    df <- final_data()
    if (is.null(df) || nrow(df) == 0) {
      return(character(0))
    } else {
      return(setdiff(colnames(df), input$kit2))  # Return all columns except the selected kit2
    }
  })
  
  
  observe({
    kits <- available_kit1_choices()
    if (length(kits) > 0) {
      selected_kit1 <- if (input$kit1 %in% kits) input$kit1 else kits[1]
      updateSelectInput(session, "kit1", choices = kits, selected = selected_kit1)
    } else {
      updateSelectInput(session, "kit1", choices = character(0))  # Update with an empty vector if no kits are available
    }
  })
  
  # Observer for updating kit2 choices based on selected kit1
  observe({
    kits <- available_kit2_choices()
    if (length(kits) > 0) {
      selected_kit2 <- if (input$kit2 %in% kits) input$kit2 else "None"
      updateSelectInput(session, "kit2", choices = c("None", kits), selected = selected_kit2)
    } else {
      updateSelectInput(session, "kit2", choices = character(0))  # Update with an empty vector if no kits are available
    }
  })
  
  # Combine the two observers for updating kit1 choices
  observe({
    req(input$process_input)
    kits <- kit_choices()
    if (length(kits) > 0) {
      selected_kit1 <- if (input$kit1 %in% kits) input$kit1 else kits[1]
      updateSelectInput(session, "kit1", choices = kits, selected = selected_kit1)  # Set the first kit as the default selected value
    } else {
      updateSelectInput(session, "kit1", choices = character(0))  # Update with an empty vector if no kits are available
    }
  })
  
  # Render the PI vs Mw plot reactively based on selected kits
  output$pi_vs_mw_plot <- renderPlot({
    req(comprehensive_table())
    df <- comprehensive_table()
    if (nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data available to display", cex = 1.5)
    } else {
      if(input$kit2 == "None"){
        visualize_data(comprehensive_table(), kit1_name = input$kit1)
      } else {
        visualize_data(comprehensive_table(), kit1_name = input$kit1, kit2_name = input$kit2)
      }
      
    }
  })

  
  # Reactive value for selected format
  selected_format_venn <- reactiveVal("png")
  selected_format_coverage <- reactiveVal("png")
  selected_format_pi_mw <- reactiveVal("png")
  
  # Update selected format based on radio button input
  observeEvent(input$coverage_format_selection, {
    selected_format_coverage(input$coverage_format_selection)
  })
  
  # Update selected format based on radio button input
  observeEvent(input$veen_format_selection, {
    selected_format_venn(input$veen_format_selection)
  })
  
  # Update selected format based on radio button input
  observeEvent(input$pi_mw_format_selection, {
    selected_format_pi_mw(input$pi_mw_format_selection)
  })

  output$download_combined_excel <- downloadHandler(
    filename = function() {
      "coverage_information.xlsx"
    },
    content = function(file) {
      req(final_data())
      table <- final_data()
      # Assuming write_combined_excel returns a list of data frames
      combined_tables <- write_combined_excel(table)
      # Write to a single Excel file with multiple sheets
      writexl::write_xlsx(combined_tables, path = file)
    }
  )
  

  # Define download handlers for Coverage Plot
  output$download_upset_plot <- downloadHandler(
    filename = function() {
      file_name <- paste("Coverage_Plot", input$coverage_format_selection, sep = ".")
      file_name
    },
    content = function(file) {
      shiny::req(final_data(), set_list())
      
      dimensions <- validate.plot.dimensions(type = "rectangular", width = input$coverage_width, height = input$coverage_height)
      dimensions <- as.list(dimensions)
      
      input_percent <- calculate_input_percent(set_list(), final_data())
      
      upset_plot <- create_upset_diagram(input_percent)
      
      ggsave(file, plot = upset_plot, width = dimensions$width, height = dimensions$height, units = "in", dpi = 300, device = selected_format_coverage())
    }
  )
  
  
  
  output$venn_diagram_download <- shiny::downloadHandler(
    filename = function() {
      file_name <- paste("Venn_Diagram", input$veen_format_selection, sep = ".")
      file_name
    },
    content = function(file) {
      shiny::req(input$comparisons)
      validate(
        need(length(input$comparisons) >= 2, "Please select at least two comparisons"),
        need(length(input$comparisons) <= 5, "Please select no more than five comparisons")
      )
      shiny::req(set_list())
      # Validate the selected figure dimensions
      dimensions <- validate.plot.dimensions(type = "rectangular", width = input$venn_width, height = input$venn_height)
      dimensions <- as.list(dimensions)
      # Create the Venn Diagram
      venn_plot <- create_venn_diagram(set_list(), input$comparisons)
      # Export the Venn Diagram
      ggsave(file, plot = venn_plot, width = dimensions$width, height = dimensions$height, units = "in", dpi = 300, device = selected_format_venn())
    }
  )
  
  # Define download handlers for PI vs Mw Plot
  output$download_pi_mw <- downloadHandler(
    filename = function() {
      paste("PI_vs_Mw_Plot", selected_format_pi_mw(), sep = ".")
    },
    content = function(file) {
      req(data())
      # Validate the selected figure dimensions
      dimensions <- validate.plot.dimensions(type = "rectangular", width = input$pi_mw_width, height = input$pi_mw_height)
      dimensions <- as.list(dimensions)
 
      # Export the PI vs Mw Plot
      ggsave(file, width = dimensions$width, height = dimensions$height, units = "in", dpi = 300, device = selected_format_pi_mw())
    }
  )
  

# Define the UI
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      uiOutput("sidebar_panel_ui")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("DataFrame View", DT::dataTableOutput("dataframe_view")),
        tabPanel("Venn Diagram", plotOutput("venn_diagram"), 
                 downloadButton("save_venn", "Save Venn Diagram")),
        tabPanel("PI vs Mw Plot", plotOutput("pi_vs_mw_plot"),
                 downloadButton("save_pi_mw_plot", "Save PI vs Mw Plot")),
        tabPanel("UpSet Plot", plotOutput("upset_plot"),
                 downloadButton("save_upset_plot", "Save UpSet Plot"))
      )
    )
  )
)

shinyApp(ui = ui, server = server)



  ####################################################################################################################################################
  
  #Once the browser's window is closed... 
  session$onSessionEnded(function() {
    #...stop the app!
    shiny::stopApp()
  })
  
}

#Run the application.
shiny::shinyApp(ui = ui, server = server)

