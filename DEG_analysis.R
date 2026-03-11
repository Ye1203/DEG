
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
    
    add_row_ui <- function(num = NULL, den = NULL) {
      row_idx <- as.character(Sys.time())
      row_idx <- gsub("[^0-9]", "", row_idx)  
      row_id <- paste0("row_", row_idx)
      
      insertUI(
        selector = paste0("#", ns("rows_container")),
        where = "beforeEnd",
        ui = div(
          id = ns(row_id),
          fluidRow(
            style = "margin-top: 5px;",
            column(3, selectInput(ns(paste0("num_", row_id)), NULL,
                                  choices = column_values(), selected = num,
                                  width = "100%")),
            column(3, selectInput(ns(paste0("den_", row_id)), NULL,
                                  choices = column_values(), selected = den,
                                  width = "100%")),
            column(2, actionButton(ns(paste0("exchange_", row_id)), "⇄ Swap",
                                   style = "color:#fff;background-color:#337ab7;border-color:#2e6da4;",
                                   width = "100%")),
            column(2, actionButton(ns(paste0("delete_", row_id)), "✕ Delete",
                                   style = "color:#fff;background-color:#d9534f;border-color:#d43f3a;",
                                   width = "100%"))
          )
        )
      )
      
      rv$rows[[row_id]] <- list(numerator = num, denominator = den)
      
      # DELETE
      observeEvent(input[[paste0("delete_", row_id)]], {
        removeUI(selector = paste0("#", ns(row_id)))
        rv$rows[[row_id]] <- NULL
      }, ignoreInit = TRUE)
      
      # SWAP
      observeEvent(input[[paste0("exchange_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        
        updateSelectInput(session, paste0("num_", row_id), selected = row$denominator)
        updateSelectInput(session, paste0("den_", row_id), selected = row$numerator)
        
        rv$rows[[row_id]]$numerator <- row$denominator
        rv$rows[[row_id]]$denominator <- row$numerator
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("num_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        rv$rows[[row_id]]$numerator <- input[[paste0("num_", row_id)]]
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("den_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        rv$rows[[row_id]]$denominator <- input[[paste0("den_", row_id)]]
      }, ignoreInit = TRUE)
      
      row_count(length(rv$rows))
    }
    
    clear_all_rows <- function() {
      lapply(names(rv$rows), function(rid) {
        removeUI(selector = paste0("#", ns(rid)), immediate = TRUE)
      })
      rv$rows <- list()
      row_count(0)
    }
    
    observeEvent(input$add_row, { add_row_ui() })
    
    observeEvent(input$add_all, {
      vals <- column_values()
      if (length(vals) > 1) {
        clear_all_rows()
        comb <- t(combn(vals, 2))
        apply(comb, 1, function(row) add_row_ui(num = row[1], den = row[2]))
      }
    })
    
    observeEvent(input$remove_all, { clear_all_rows() })
    
    return(reactive(rv$rows))
  })
}

step2UI <- function(id) {
  ns <- NS(id)
  
  tagList(
      h4("Sample Comparisons:"),
      uiOutput(ns("add_sample_comp_ui")),
      br(), br(),
      div(id = ns("sample_cards_container"), class = "row"),
      
      fluidRow(column(12, hr())),
      
      h4("Cluster Comparisons:"),
      uiOutput(ns("add_cluster_comp_ui")),
      br(), br(),
      div(id = ns("cluster_cards_container"), class = "row"),
      
    fluidRow(
      column(12,hr()),
      column(
        12,
        actionButton(
          "back_step",
          "Back to Step 1",
          class = "btn-warning",
          style = "width:200px;"
        )
      ))
  )
}

# ----------------- Step2 Server Module -----------------
step2Server <- function(id, sample_rows, cluster_rows, sample_meta, cluster_meta) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$add_sample_comp_ui <- renderUI({
      rows <- sample_rows()
      if (is.null(rows) || length(rows) == 0){
        div(
          style = "color:#555; font-style:italic;",
          "No sample comparisons found. You can add comparisons in Step 1 using the 'Back to Step 1' button."
        )
      }else{
        actionButton(
          ns("add_sample_comp"),
          "Add Comparison",
          class = "btn-primary"
        )
      }
    })
    
    output$add_cluster_comp_ui <- renderUI({
      rows <- cluster_rows()
      if (is.null(rows) || length(rows) == 0){
        div(
          style = "color:#555; font-style:italic;",
          "No celltype cluster comparisons found. You can add comparisons in Step 1 using the 'Back to Step 1' button."
        )
      }else{
        actionButton(
          ns("add_cluster_comp"),
          "Add Comparison",
          class = "btn-primary"
        )
      }
    })
    
    # Helper: Create card UI
    createCard <- function(card_title, available_comps, meta_choices, card_id, is_sample = TRUE) {
      # Create checkbox choices from all available comparisons
      comp_choices <- list()
      for (comp in available_comps) {
        comp_choices[[comp$name]] <- comp$id
      }
      
      div(
        id = ns(card_id),  # Add ID to the main div for easy removal
        class = "col-sm-4",  
        style = "margin-bottom: 15px;",
        wellPanel(
          h5(card_title, style = "margin-top: 0; font-weight: bold;"),
          
          selectInput(
            ns(paste0(card_id, "_meta_selector")),  # Use underscore instead of hyphen for easier referencing
            label = ifelse(is_sample, "Select clusters to compare in:", "Select samples to compare in:"),
            choices = meta_choices,
            selected = NULL,
            multiple = TRUE,
            width = "100%"
          ),
          
          checkboxGroupInput(
            ns(paste0(card_id, "_comparison_checkboxes")),
            label = "Select comparisons to perform:",
            choices = comp_choices,
            selected = NULL,
            width = "100%"
          ),
          
          fluidRow(
            column(6, actionButton(ns(paste0(card_id, "_select_all")), "Select All", 
                                   class = "btn-sm btn-primary", width = "100%")),
            column(6, actionButton(ns(paste0(card_id, "_clear_all")), "Clear", 
                                   class = "btn-sm btn-warning", width = "100%"))
          ),
          br(),
          actionButton(ns(paste0(card_id, "_delete")), "Delete", 
                       class = "btn-danger", style = "width:100%; margin-top:5px;")
        )
      )
    }
    
    # Get all valid sample comparisons from Step1
    all_sample_comps <- reactive({
      comps <- sample_rows()
      comps <- comps[sapply(comps, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
      
      comp_list <- list()
      for (i in seq_along(comps)) {
        comp <- comps[[i]]
        comp_id <- paste0("S", i, "_", comp$numerator, "_vs_", comp$denominator)
        comp_name <- paste0("S", i, ": ", comp$numerator, " vs ", comp$denominator)
        comp_list[[i]] <- list(id = comp_id, name = comp_name, 
                               numerator = comp$numerator, denominator = comp$denominator)
      }
      comp_list
    })
    
    # Get all valid cluster comparisons from Step1
    all_cluster_comps <- reactive({
      comps <- cluster_rows()
      comps <- comps[sapply(comps, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
      
      comp_list <- list()
      for (i in seq_along(comps)) {
        comp <- comps[[i]]
        comp_id <- paste0("C", i, "_", comp$numerator, "_vs_", comp$denominator)
        comp_name <- paste0("C", i, ": ", comp$numerator, " vs ", comp$denominator)
        comp_list[[i]] <- list(id = comp_id, name = comp_name, 
                               numerator = comp$numerator, denominator = comp$denominator)
      }
      comp_list
    })
    
    # ----------------- Sample Cards Management -----------------
    sample_cards <- reactiveValues(ids = character(0), data = list())
    
    # Add new sample comparison card
    observeEvent(input$add_sample_comp, {
      sample_comps <- all_sample_comps()
      if (length(sample_comps) == 0) {
        showModal(modalDialog(
          title = "No Sample comparisons", 
          "Please add at least one sample comparison in Step1."
        ))
        return()
      }
      
      card_id <- paste0("sample_card_", length(sample_cards$ids) + 1)
      sample_cards$ids <- c(sample_cards$ids, card_id)
      
      sample_cards$data[[card_id]] <- list(
        type = "sample",
        selected_comps = character(0),
        selected_meta = character(0)
      )
      
      # Insert the card UI
      insertUI(
        selector = paste0("#", ns("sample_cards_container")),
        ui = createCard(
          card_title = "Sample Comparisons",
          available_comps = sample_comps,
          meta_choices = sample_meta(),
          card_id = card_id,
          is_sample = TRUE
        ),
        immediate = TRUE
      )
      
      # Set up observers for this card
      setupCardObservers(card_id, "sample", sample_cards, all_sample_comps)
    })
    
    # ----------------- Cluster Cards Management -----------------
    cluster_cards <- reactiveValues(ids = character(0), data = list())
    
    # Add new cluster comparison card
    observeEvent(input$add_cluster_comp, {
      cluster_comps <- all_cluster_comps()
      if (length(cluster_comps) == 0) {
        showModal(modalDialog(
          title = "No Cluster comparisons", 
          "Please add at least one cluster comparison in Step1."
        ))
        return()
      }
      
      card_id <- paste0("cluster_card_", length(cluster_cards$ids) + 1)
      cluster_cards$ids <- c(cluster_cards$ids, card_id)
      
      cluster_cards$data[[card_id]] <- list(
        type = "cluster",
        selected_comps = character(0),
        selected_meta = character(0)
      )
      
      # Insert the card UI
      insertUI(
        selector = paste0("#", ns("cluster_cards_container")),
        ui = createCard(
          card_title = "Cluster Comparisons",
          available_comps = cluster_comps,
          meta_choices = cluster_meta(),
          card_id = card_id,
          is_sample = FALSE
        ),
        immediate = TRUE
      )
      
      # Set up observers for this card
      setupCardObservers(card_id, "cluster", cluster_cards, all_cluster_comps)
    })
    
    # Helper function to set up card observers
    setupCardObservers <- function(card_id, card_type, cards_reactive, all_comps_func) {
      
      # Create input ID patterns
      delete_id <- paste0(card_id, "_delete")
      select_id <- paste0(card_id, "_select_all")
      clear_id <- paste0(card_id, "_clear_all")
      checkbox_id <- paste0(card_id, "_comparison_checkboxes")
      meta_id <- paste0(card_id, "_meta_selector")
      
      # Delete observer
      observeEvent(input[[delete_id]], {
        removeUI(selector = paste0("#", ns(card_id)))
        cards_reactive$ids <- setdiff(cards_reactive$ids, card_id)
        cards_reactive$data[[card_id]] <- NULL
      }, ignoreInit = TRUE)
      
      # Select All observer
      observeEvent(input[[select_id]], {
        all_comp_choices <- sapply(all_comps_func(), function(x) x$id)
        
        updateCheckboxGroupInput(
          session,
          checkbox_id,
          selected = all_comp_choices
        )
        
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- all_comp_choices
        }
      }, ignoreInit = TRUE)
      
      # Clear All observer
      observeEvent(input[[clear_id]], {
        updateCheckboxGroupInput(
          session,
          checkbox_id,
          selected = character(0)
        )
        
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- character(0)
        }
      }, ignoreInit = TRUE)
      
      # Monitor checkbox changes
      observeEvent(input[[checkbox_id]], {
        selected <- input[[checkbox_id]]
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- selected
        }
      }, ignoreNULL = FALSE)
      
      # Monitor meta selector changes
      observeEvent(input[[meta_id]], {
        selected <- input[[meta_id]]
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_meta <- selected
        }
      }, ignoreNULL = FALSE)
    }
    
    # Return all selected data
    return(
      list(
        sample_cards_data = reactive(sample_cards$data),
        cluster_cards_data = reactive(cluster_cards$data),
        all_sample_comps = all_sample_comps,
        all_cluster_comps = all_cluster_comps
      )
    )
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
    }else if(input$sample_column == "" | input$cluster_column == ""){
      showModal(modalDialog(
        title = "Sample meta.data or Cluster meta.data was selected.",
        "Please make sure you have select meta.data for both sample and cluster.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    }
    else {
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
        div(id = "step2_ui", style="display:none;", step2UI("step2_module"))
      )
    }
  })
  
  step2Server(
    "step2_module",
    sample_rows = sample_rows,
    cluster_rows = cluster_rows,
    sample_meta = cluster_values, 
    cluster_meta = sample_values  
  )

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