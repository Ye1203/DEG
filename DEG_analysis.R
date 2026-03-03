library(shiny)

# ------------------- Parameters -------------------
DEBUG <- TRUE  # Set TRUE to show debug button

# ============================================================
# -------------------- Utility Functions ---------------------
# ============================================================

load_data_object <- function(path) {
  validate(need(file.exists(path), "File does not exist"))
  readRDS(path)
}

extract_meta_columns <- function(data) {
  if (!"meta.data" %in% slotNames(data)) return(NULL)
  colnames(data@meta.data)
}

extract_unique_values <- function(data, column_name) {
  if (is.null(data) || is.null(column_name)) return(NULL)
  unique(as.character(data@meta.data[[column_name]]))
}

generate_combinations <- function(values) {
  if (length(values) < 2) return(data.frame())
  comb <- expand.grid(numerator = values, denominator = values, stringsAsFactors = FALSE)
  comb <- comb[comb$numerator != comb$denominator, ]
  rownames(comb) <- NULL
  comb
}

# ============================================================
# -------------------- Comparison Module ---------------------
# ============================================================

comparisonTableUI <- function(id, title_prefix) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("table_ui"))
  )
}

comparisonTableServer <- function(id, column_values, prefix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(rows = list())
    
    add_row <- function(num = NULL, den = NULL) {
      rv$rows[[length(rv$rows) + 1]] <- list(numerator = num, denominator = den)
    }
    
    observeEvent(input$add_all, {
      rv$rows <- list()
      comb <- generate_combinations(column_values())
      if (nrow(comb) > 0) apply(comb, 1, function(x) add_row(x[1], x[2]))
    })
    
    observeEvent(input$remove_all, { rv$rows <- list() })
    observeEvent(input$add_row, { add_row() })
    
    output$table_ui <- renderUI({
      if (length(rv$rows) == 0) {
        actionButton(ns("add_row"), "Add Row")
      } else {
        tagList(
          lapply(seq_along(rv$rows), function(i) {
            fluidRow(
              column(2, strong(paste0(prefix, i))),
              column(3, selectInput(ns(paste0("num_", i)), "Treatment (numerator)",
                                    choices = column_values(), selected = rv$rows[[i]]$numerator)),
              column(3, selectInput(ns(paste0("den_", i)), "Control (denominator)",
                                    choices = column_values(), selected = rv$rows[[i]]$denominator)),
              column(2, actionButton(ns(paste0("delete_", i)), "Delete"))
            )
          }),
          br(),
          actionButton(ns("add_row"), "Add Row")
        )
      }
    })
    
    observe({
      lapply(seq_along(rv$rows), function(i) {
        observeEvent(input[[paste0("delete_", i)]], {
          rv$rows <- rv$rows[-i]
        }, ignoreInit = TRUE)
      })
    })
    
    return(reactive(rv$rows))
  })
}

# ============================================================
# ------------------------- UI -------------------------------
# ============================================================

ui <- fluidPage(
  titlePanel("Analysis Application"),
  
  # Data path input
  fluidRow(
    column(9, textInput("data_path", "Input data path")),
    column(3, br(), actionButton("load_data", "Load Data"))
  ),
  
  tags$hr(),
  
  # Debug button if DEBUG = TRUE
  conditionalPanel(
    condition = paste0(DEBUG),
    absolutePanel(bottom = 10, left = 10,
                  actionButton("debug_btn", "DEBUG"))
  ),
  
  # Hidden step state
  tags$script(HTML("
    Shiny.setInputValue('step_state', 1);
    Shiny.addCustomMessageHandler('setStep', function(value) {
      Shiny.setInputValue('step_state', value, {priority: 'event'});
    });
  ")),
  
  # Tabset
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel(
      "Data processing",
      
      # Show UI only if data_obj() is not NULL
      uiOutput("data_processing_ui")
    ),
    
    tabPanel("Visualization", h3("Visualization content will be defined later."))
  )
)

# ============================================================
# ------------------------ SERVER ----------------------------
# ============================================================

server <- function(input, output, session) {
  
  data_obj <- reactiveVal(NULL)
  
  # ----------------- DEBUG button -----------------
  if (DEBUG) {
    observeEvent(input$debug_btn, { browser() })
  }
  
  # ----------------- Load Data -----------------
  observeEvent(input$load_data, {
    showModal(modalDialog(title = "Loading Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      withProgress(message = "Reading file...", value = 0, {
        incProgress(0.3)
        obj <- load_data_object(input$data_path)
        data_obj(obj)
        incProgress(0.6)
        cols <- extract_meta_columns(obj)
        updateSelectInput(session, "sample_column", choices = cols)
        updateSelectInput(session, "cluster_column", choices = cols)
        # Initialize modules with one row
        session$sendInputMessage("sample_table-add_row", list(value = 1))
        session$sendInputMessage("cluster_table-add_row", list(value = 1))
        incProgress(1)
      })
      removeModal()
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # ----------------- Reactive unique values -----------------
  sample_values <- reactive({ req(data_obj(), input$sample_column)
    extract_unique_values(data_obj(), input$sample_column) })
  cluster_values <- reactive({ req(data_obj(), input$cluster_column)
    extract_unique_values(data_obj(), input$cluster_column) })
  
  # ----------------- Modules -----------------
  sample_rows <- comparisonTableServer("sample_table", sample_values, prefix = "S")
  cluster_rows <- comparisonTableServer("cluster_table", cluster_values, prefix = "C")
  
  # ----------------- Top buttons -----------------
  observeEvent(input$sample_add_all, { session$sendInputMessage("sample_table-add_all", list(value = 1)) })
  observeEvent(input$sample_remove_all, { session$sendInputMessage("sample_table-remove_all", list(value = 1)) })
  observeEvent(input$cluster_add_all, { session$sendInputMessage("cluster_table-add_all", list(value = 1)) })
  observeEvent(input$cluster_remove_all, { session$sendInputMessage("cluster_table-remove_all", list(value = 1)) })
  
  # ----------------- Step navigation -----------------
  observeEvent(input$next_step, { session$sendCustomMessage("setStep", 2) })
  observeEvent(input$back_step, { session$sendCustomMessage("setStep", 1) })
  
  # ----------------- Data Processing UI -----------------
  output$data_processing_ui <- renderUI({
    req(data_obj())  # Only show after data loaded
    
    # Step 1 UI
    conditionalPanel(
      condition = "input.step_state == 1",
      fluidRow(
        column(6, selectInput("sample_column", "meta.data sample column", choices = NULL)),
        column(6, selectInput("cluster_column", "meta.data cluster column", choices = NULL))
      ),
      hr(),
      h3("Sample Analysis"),
      fluidRow(
        column(6, actionButton("sample_add_all", "Add All Comparison")),
        column(6, actionButton("sample_remove_all", "Remove All Comparison"))
      ),
      comparisonTableUI("sample_table", "S"),
      hr(),
      h3("Cluster Analysis"),
      fluidRow(
        column(6, actionButton("cluster_add_all", "Add All Comparison")),
        column(6, actionButton("cluster_remove_all", "Remove All Comparison"))
      ),
      comparisonTableUI("cluster_table", "C"),
      hr(),
      actionButton("next_step", "NEXT STEP")
    )
  })
}

shinyApp(ui, server)