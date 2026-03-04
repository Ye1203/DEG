
library(shiny)
library(shinyjs)
library(Seurat)
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
    rv <- reactiveValues(rows = list())
    row_count <- reactiveVal(0)
    
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
          )
        ),
        br(),
        fluidRow(
          column(3, strong("Treatment (numerator)"), style = "text-align:center;"),
          column(3, strong("Control (denominator)"), style = "text-align:center;"),
          column(2, strong("Swap N/E")),
          column(2, strong("Action"))
        ),
        br(),
        div(id = ns("rows_container"))
      )
    })
    
    add_row_ui <- function(num = NULL, den = NULL, index = NULL) {
      if (is.null(index)) {
        index <- row_count() + 1
        row_count(index)
      }
      
      row_id <- paste0("row_", index)
      
      insertUI(
        selector = paste0("#", ns("rows_container")),
        where = "beforeEnd",
        ui = div(
          id = ns(row_id),
          fluidRow(
            style = "margin-top: 5px;",
            column(3, selectInput(ns(paste0("num_", index)), NULL,
                                  choices = column_values(), selected = num,
                                  width = "100%")),
            column(3, selectInput(ns(paste0("den_", index)), NULL,
                                  choices = column_values(), selected = den,
                                  width = "100%")),
            column(2, actionButton(ns(paste0("exchange_", index)), "⇄ Swap",
                                   style = "color:#fff;background-color:#337ab7;border-color:#2e6da4;",
                                   width = "100%")),
            column(2, actionButton(ns(paste0("delete_", index)), "✕ Delete",
                                   style = "color:#fff;background-color:#d9534f;border-color:#d43f3a;",
                                   width = "100%"))
          )
        )
      )
      
      new_rows <- rv$rows
      new_rows[[index]] <- list(numerator = num, denominator = den)
      rv$rows <- new_rows
      
      local({
        row_idx <- index
        
        observeEvent(input[[paste0("delete_", row_idx)]], {
          removeUI(selector = paste0("#", ns(paste0("row_", row_idx))))
          
          new_rows <- rv$rows
          new_rows[[row_idx]] <- NULL
          rv$rows <- new_rows[!sapply(new_rows, is.null)]
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("exchange_", row_idx)]], {
          current_rows <- isolate(rv$rows)
          if (!is.null(current_rows[[row_idx]])) {
            current_num <- current_rows[[row_idx]]$numerator
            current_den <- current_rows[[row_idx]]$denominator
            
            updateSelectInput(session, paste0("num_", row_idx), selected = current_den)
            updateSelectInput(session, paste0("den_", row_idx), selected = current_num)
            
            new_rows <- rv$rows
            new_rows[[row_idx]]$numerator <- current_den
            new_rows[[row_idx]]$denominator <- current_num
            rv$rows <- new_rows
          }
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("num_", row_idx)]], {
          new_rows <- rv$rows
          if (!is.null(new_rows[[row_idx]])) {
            new_rows[[row_idx]]$numerator <- input[[paste0("num_", row_idx)]]
            rv$rows <- new_rows
          }
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("den_", row_idx)]], {
          new_rows <- rv$rows
          if (!is.null(new_rows[[row_idx]])) {
            new_rows[[row_idx]]$denominator <- input[[paste0("den_", row_idx)]]
            rv$rows <- new_rows
          }
        }, ignoreInit = TRUE)
      })
    }
    
    clear_all_rows <- function() {
      for (i in seq_along(rv$rows)) {
        if (!is.null(rv$rows[[i]])) {
          removeUI(selector = paste0("#", ns(paste0("row_", i))), immediate = TRUE)
        }
      }
      rv$rows <- list()
      row_count(0)
    }
    
    observeEvent(input$add_row, {
      add_row_ui()
    })
    
    observeEvent(input$add_all, {
      values <- column_values()
      if (length(values) > 1) {
        clear_all_rows()
        comb <- t(combn(values, 2))
        for (i in 1:nrow(comb)) {
          add_row_ui(num = comb[i, 1], den = comb[i, 2], index = i)
        }
        row_count(nrow(comb))
      }
    })
    
    observeEvent(input$remove_all, {
      clear_all_rows()
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