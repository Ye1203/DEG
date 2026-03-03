library(shiny)
library(shinyjs)

# ============================================================
# -------------------- Utility Functions ---------------------
# ============================================================

load_data_object <- function(path) {
  validate(need(file.exists(path), "File does not exist"))
  readRDS(path)
}

extract_meta_columns <- function(data) {
  if (!"meta.data" %in% slotNames(data)) return(NULL)
  c("",colnames(data@meta.data))
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
  uiOutput(ns("table_ui"))
}

comparisonTableServer <- function(id, column_values, prefix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(rows = list(list(numerator = NULL, denominator = NULL)))
    
    add_row <- function(num = NULL, den = NULL) {
      rv$rows[[length(rv$rows) + 1]] <- list(numerator = num, denominator = den)
    }
    
    observeEvent(input$add_all, {
      values <- column_values()
      if (length(values) > 1) {
        comb <- t(combn(values, 2))
        rv$rows <- lapply(1:nrow(comb), function(i) {
          list(numerator = comb[i, 1], denominator = comb[i, 2])
        })
      } else {
        rv$rows <- list(list(numerator = NULL, denominator = NULL))
      }
    })
    
    observeEvent(input$remove_all, {
      rv$rows <- list(list(numerator = NULL, denominator = NULL))
    })
    
    observeEvent(input$add_row, {
      add_row()
    })
    
    delete_observers <- list()
    exchange_observers <- list()
    num_observers <- list()
    den_observers <- list()
    
    observe({
      current_rows <- length(rv$rows)
      existing_observers <- length(delete_observers)
      
      if (current_rows < existing_observers) {
        for (i in (current_rows + 1):existing_observers) {
          if (!is.null(delete_observers[[i]])) delete_observers[[i]]$destroy()
          if (!is.null(exchange_observers[[i]])) exchange_observers[[i]]$destroy()
          if (!is.null(num_observers[[i]])) num_observers[[i]]$destroy()
          if (!is.null(den_observers[[i]])) den_observers[[i]]$destroy()
        }
        delete_observers <<- delete_observers[1:current_rows]
        exchange_observers <<- exchange_observers[1:current_rows]
        num_observers <<- num_observers[1:current_rows]
        den_observers <<- den_observers[1:current_rows]
      }
      
      if (current_rows > length(delete_observers)) {
        for (i in (length(delete_observers) + 1):current_rows) {
          delete_observers[[i]] <<- observeEvent(input[[paste0("delete_", i)]], {
            current_rows <- rv$rows
            if (i <= length(current_rows)) {
              new_rows <- current_rows[-i]
              if (length(new_rows) == 0) {
                new_rows <- list(list(numerator = NULL, denominator = NULL))
              }
              rv$rows <- new_rows
            }
          }, ignoreInit = TRUE, autoDestroy = FALSE)
          
          exchange_observers[[i]] <<- observeEvent(input[[paste0("exchange_", i)]], {
            if (i <= length(rv$rows)) {
              current_num <- rv$rows[[i]]$numerator
              current_den <- rv$rows[[i]]$denominator
              rv$rows[[i]]$numerator <- current_den
              rv$rows[[i]]$denominator <- current_num
            }
          }, ignoreInit = TRUE, autoDestroy = FALSE)
          
          num_observers[[i]] <<- observeEvent(input[[paste0("num_", i)]], {
            if (i <= length(rv$rows)) {
              rv$rows[[i]]$numerator <- input[[paste0("num_", i)]]
            }
          }, ignoreInit = TRUE, autoDestroy = FALSE)
          
          den_observers[[i]] <<- observeEvent(input[[paste0("den_", i)]], {
            if (i <= length(rv$rows)) {
              rv$rows[[i]]$denominator <- input[[paste0("den_", i)]]
            }
          }, ignoreInit = TRUE, autoDestroy = FALSE)
        }
      }
    })
    
    output$table_ui <- renderUI({
      req(column_values())
      tagList(
        fluidRow(
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:20px;",
            actionButton(ns("add_all"), "Add All Comparison", class = "btn-primary", width = "100%")
          ),
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:5px;",
            actionButton(ns("remove_all"), "Remove All", class = "btn-warning", width = "100%")
          ),
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:5px;",
            actionButton(ns("add_row"), "Add One Row", class = "btn-success", width = "100%")
          )),
        br(),
        fluidRow(
          column(1, strong("Index")),
          column(3, strong("Treatment (numerator)"), style = "text-align:center;"),
          column(3, strong("Control (denominator)"), style = "text-align:center;"),
          column(2, strong("Swap N/E")),
          column(2, strong("Action"))
        ),
        lapply(seq_along(rv$rows), function(i) {
          fluidRow(
            style = "margin-top: 5px;",
            column(1, div(style = "margin-top: 2.5px;",
                          strong(paste0(prefix, i)))),
            column(3, selectInput(ns(paste0("num_", i)), NULL,
                                  choices = column_values(), selected = rv$rows[[i]]$numerator,
                                  width = "100%")),
            column(3, selectInput(ns(paste0("den_", i)), NULL,
                                  choices = column_values(), selected = rv$rows[[i]]$denominator,
                                  width = "100%")),
            column(2, actionButton(ns(paste0("exchange_", i)), "⇄ Swap",
                                   style = "color:#fff;background-color:#337ab7;border-color:#2e6da4;",
                                   width = "100%")),
            column(2, actionButton(ns(paste0("delete_", i)), "✕ Delete",
                                   style = "color:#fff;background-color:#d9534f;border-color:#d43f3a;",
                                   width = "100%"))
          )
        })
      )
    })
    
    return(reactive(rv$rows))
  })
}

# ============================================================
# ------------------------- UI -------------------------------
# ============================================================

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Analysis Application"),
  
  h4("Input Data Path"),
  fluidRow(
    column(9, textInput("data_path", value = "/projectnb/wax-es/00_shinyapp/Clustering/file/G193_Male.rds", label = NULL, width = "100%")),
    column(3, actionButton("load_data", "Load Data", class = "btn-success"))
  ),
  
  tabsetPanel(
    id = "main_tabs",
    tabPanel("Data processing", uiOutput("data_processing_ui")),
    tabPanel("Visualization", uiOutput("visualization_ui"))
  ),
  br(),br(),br(),br(),br(),br()
)

# ============================================================
# ------------------------ SERVER ----------------------------
# ============================================================

server <- function(input, output, session) {
  
  data_obj <- reactiveVal(NULL)
  
  # ----------------- Load Data -----------------
  observeEvent(input$load_data, {
    showModal(modalDialog(title = "Loading Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$data_path)
      data_obj(obj)
      removeModal()
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # ----------------- Reactive unique values -----------------
  sample_values <- reactive({
    req(data_obj(), input$sample_column)
    extract_unique_values(data_obj(), input$sample_column)
  })
  
  cluster_values <- reactive({
    req(data_obj(), input$cluster_column)
    extract_unique_values(data_obj(), input$cluster_column)
  })
  
  # ----------------- Modules -----------------
  sample_rows <- comparisonTableServer("sample_table", sample_values, prefix = "S")
  cluster_rows <- comparisonTableServer("cluster_table", cluster_values, prefix = "C")
  
  # ----------------- Step navigation -----------------
  observeEvent(input$next_step, {
    sample_comparisons <- sample_rows()
    cluster_comparisons <- cluster_rows()
    
    # Check if at least one comparison is selected (not NULL)
    has_sample <- any(sapply(sample_comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator)))
    has_cluster <- any(sapply(cluster_comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator)))
    
    if (!has_sample && !has_cluster) {
      showModal(modalDialog(
        title = "No Comparisons Selected",
        "Please add at least one comparison for Sample or Cluster analysis before proceeding.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    } else {
      shinyjs::hide("step1_ui")
      shinyjs::show("step2_ui")
    }
  })
  
  observeEvent(input$back_step, {
    shinyjs::show("step1_ui")
    shinyjs::hide("step2_ui")
  })
  
  # ----------------- Data Processing UI -----------------
  output$data_processing_ui <- renderUI({
    if (is.null(data_obj())) {
      h4("Please load data first.")
    } else {
      tagList(
        # Step 1 (original content wrapped in div)
        div(id = "step1_ui",
            fluidRow(
              column(6, selectInput("sample_column", "meta.data sample column", choices = extract_meta_columns(data_obj()))),
              column(6, selectInput("cluster_column", "meta.data cluster column", choices = extract_meta_columns(data_obj())))
            ),
            hr(),
            h3("Sample Analysis"),
            comparisonTableUI("sample_table", "S"),
            hr(),
            h3("Cluster Analysis"),
            comparisonTableUI("cluster_table", "C"),
            hr(),
            actionButton("next_step", "NEXT STEP", class = "btn-primary", style = "width:200px;")
        ),
        
        # Step 2 (hidden initially)
        div(id = "step2_ui", style = "display: none;",
            h4("Sample Comparisons:"),
            uiOutput("sample_comparisons_summary"),
            hr(),
            h4("Cluster Comparisons:"),
            uiOutput("cluster_comparisons_summary"),
            hr(),
            actionButton("back_step", "Back to Step 1", class = "btn-warning", style = "width:200px;")
        )
      )
    }
  })
  
  # ----------------- Generate comparison summaries -----------------
  output$sample_comparisons_summary <- renderUI({
    comparisons <- sample_rows()
    # Filter out rows with NULL values
    valid_comparisons <- comparisons[sapply(comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
    
    if (length(valid_comparisons) == 0) {
      return(p("No valid sample comparisons selected."))
    }
    
    tags$ul(
      lapply(seq_along(valid_comparisons), function(i) {
        comp <- valid_comparisons[[i]]
        tags$li(sprintf("Comparison %d: %s vs %s", i, comp$numerator, comp$denominator))
      })
    )
  })
  
  output$cluster_comparisons_summary <- renderUI({
    comparisons <- cluster_rows()
    # Filter out rows with NULL values
    valid_comparisons <- comparisons[sapply(comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
    
    if (length(valid_comparisons) == 0) {
      return(p("No valid cluster comparisons selected."))
    }
    
    tags$ul(
      lapply(seq_along(valid_comparisons), function(i) {
        comp <- valid_comparisons[[i]]
        tags$li(sprintf("Comparison %d: %s vs %s", i, comp$numerator, comp$denominator))
      })
    )
  })
  
  # ----------------- Visualization UI -----------------
  output$visualization_ui <- renderUI({
    if (is.null(data_obj())) {
      h4("Please load data first.")
    } else {
      tagList(
        h3("Visualization will be implemented here")
      )
    }
  })
}

shinyApp(ui, server)