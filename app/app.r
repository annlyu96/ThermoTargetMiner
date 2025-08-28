## app.R

#setwd("I:/TTM/MyShinyApp/app")

library(shiny)
library(shinydashboard)
library(plotly)
library(tidyverse)
library(readxl)
library(DT)
library(shinycssloaders)
library(ComplexUpset)
library(clusterProfiler)
library(org.Hs.eg.db)
library(future)
library(promises)
plan(multisession)
library(ggplot2)
library(patchwork)

ui <- dashboardPage(
  dashboardHeader(
    title = tags$div(
      tags$img(
        src = "RZlab_ThermoTargetMiner.jpg",
        height = 80, 
        align = "left"
      )
    ),
    titleWidth = 250,
    tags$li(
      class = "dropdown",
      tags$span(
        "ThermoTargetMiner",
        style = "color: white; font-size: 30px; font-weight: bold; 
               line-height: 60px; margin-left: 15px;"
      )
    )
  ),
  
  dashboardSidebar(
    width = 250,
    selectInput("dataset", label = h3("Select your proteomics data set"),
                choices = list(
                  "A549 cell lysate",
                  "A549 intact cell",
                  "H82 cell lysate",
                  "H82 intact cell"),
                selected = NULL),   
    
    selectInput("drugname", label = h3("Select Drug Name"),
                choices = NULL, selected = NULL),
    
    sliderInput("threshold", 
                label = h3("Set Threshold"), 
                min = 0, max = 10, value = 4.5, step = 0.1)  # set at 4.5
  ),
  
  dashboardBody(
    tags$style(HTML('
    .main-header {
      height: 80px !important;
    }
    .main-header .logo > img {
      object-fit: contain;
      height: 100%;
      width: auto;
      display: block;
    }
    .main-header .logo {
      height: 80px !important;
      padding-left: 0 !important;
      margin-left: 0 !important;
    }
    .main-header .navbar {
      height: 80px !important;
      min-height: 80px !important;
      display: flex;
      align-items: center;
      padding-left: 10px;
    }
    body, .content-wrapper, .main-header, .main-sidebar {
      background-color: #FFFFFF !important;
      color: #4F0433 !important;
    }
    .skin-blue .main-header .logo,
    .skin-blue .main-header .navbar {
      background-color: #4F0433 !important;
    }
    .skin-blue .main-sidebar {
    background-color: #EDF4F4 !important; 
    padding-top: 80px !important;
    box-sizing: border-box;
    }
    .main-sidebar h3 {
      color: #4F0433 !important;  
    }
    .sidebar-menu > li.active > a {
      background-color: #FF876F !important;
      color: #4F0433 !important;
    }
    .btn-primary, .nav-tabs-custom > .nav-tabs > li.active > a {
      background-color: #FF876F !important;
      color: #4F0433 !important;
      border-color: #870052 !important;
    }
  ')),
    tabsetPanel(type = "tabs",
                tabPanel("ThermoTargetMiner",
                         h3("Selected OPLS-DA Table"),
                         DTOutput("opls_table") %>% withSpinner(),
                         downloadButton("download_opls", "Download Full OPLS-DA Table"),
                         
                         h3("OPLS-DA & Transformed Data Scatter Plots and Bar Chart"),
                         fluidRow(
                           column(5, plotlyOutput("scatter1") %>% withSpinner()),
                           column(5, plotlyOutput("scatter2") %>% withSpinner()),
                           column(2, plotlyOutput("bar_chart") %>% withSpinner())
                         ),
                         
                         h3("GO Enrichment Analysis (based on genes exceeding the threshold in Transformed OPLS-DA)"),
                         uiOutput("go_status"),
                         
                         fluidRow(
                           column(4, plotOutput("plot_go_bp") %>% withSpinner()),
                           column(4, plotOutput("plot_go_cc") %>% withSpinner()),
                           column(4, plotOutput("plot_go_mf") %>% withSpinner())
                         ),
                         
                         h3("UpSet Plot of Drug Effect Overlaps"),
                         plotOutput("upset_plot") %>% withSpinner(),
                         
                         h3("Pro-target Table"),
                         dataTableOutput("protarget_table") %>% withSpinner(),
                         downloadButton("download_protarget", "Download Full Pro-target Table")
                )
    )
  )
)

server <- function(input, output, session) {
  
  # Populate drugname choices dynamically based on folder
  observe({
    drugs_dir <- "data/drugs"
    if (dir.exists(drugs_dir)) {
      files <- list.files(drugs_dir, pattern = "\\.xlsx$", full.names = FALSE)
      drugnames <- tools::file_path_sans_ext(files)
      updateSelectInput(session, "drugname", choices = drugnames)
    }
  })
  
  # OPLS-DA file reading 
  opls_data <- reactive({
    req(input$dataset, input$drugname)
    file_path <- file.path("data", "OPLS-DA", input$dataset, paste0(input$drugname, ".xlsx"))
    if (file.exists(file_path)) {
      df <- read_excel(file_path)
      colnames(df)[1:3] <- c("Gene names", "pq", "poso")
      n <- nrow(df)
      if (n >= 2) {
        df$`Gene names`[n-1] <- "All other drugs"
        df$`Gene names`[n] <- input$drugname
      }
      return(df)
    } else {
      return(NULL)
    }
  })
  
  output$opls_table <- DT::renderDT({
    req(opls_data())
    datatable(opls_data())
  })
  
  output$download_opls <- downloadHandler(
    filename = function() {
      paste0(input$drugname, "_", gsub(" ", "_", input$dataset), "_OPLSDA.xlsx")
    },
    content = function(file) {
      df <- opls_data()
      if (!is.null(df)) {
        writexl::write_xlsx(df, path = file)
      }
    }
  )
  
  # Transformed OPLS-DA data（change to 0 if < threshold）
  transformed_data <- reactive({
    req(input$dataset, input$drugname, input$threshold)
    file_path <- file.path("data", "Transformed OPLS-DA", paste0(input$dataset, ".xlsx"))
    if (file.exists(file_path)) {
      df <- read_excel(file_path)
      if (all(c("Gene names", input$drugname) %in% colnames(df))) {
        df <- df %>%
          dplyr::select(`Gene names`, all_of(input$drugname)) %>%
          mutate(y = 0) %>%
          mutate(across(where(is.numeric), ~ ifelse(. < input$threshold, 0, .)))
        return(df)
      } else {
        NULL
      }
    } else {
      NULL
    }
  })
  
  # ----- Scatter 1: OPLS-DA（threshold and color decided by transformed OPLS-DA） -----
  output$scatter1 <- renderPlotly({
    df <- opls_data(); req(df)
    trans_df <- transformed_data(); req(trans_df)
    
    colnames(df)[1:3] <- c("Gene names", "pq", "poso")
    colnames(trans_df)[1] <- "Gene names"
    
    n <- nrow(df)
    df$`Gene names`[n-1] <- "All other drugs"
    df$`Gene names`[n] <- input$drugname
    special_genes <- df$`Gene names`[(n-1):n]
    
    value_col <- input$drugname
    
    top10_genes <- trans_df %>%
      filter(!is.na(.data[[value_col]])) %>%
      arrange(desc(.data[[value_col]])) %>%
      slice_head(n = 10) %>%
      pull(`Gene names`)
    
    df_joined <- df %>%
      left_join(
        trans_df %>% dplyr::select(`Gene names`, all_of(value_col)),
        by = "Gene names"
      )
    
    # apply threshold
    genes_over_thr <- df_joined %>%
      filter(!is.na(.data[[value_col]]), .data[[value_col]] > input$threshold) %>%
      pull(`Gene names`)
    
    df_joined <- df_joined %>%
      mutate(
        gene_type = case_when(
          `Gene names` %in% special_genes ~ "special",
          `Gene names` %in% top10_genes & pq < 0 ~ "top10_neg",
          `Gene names` %in% top10_genes & pq >= 0 ~ "top10_pos",
          `Gene names` %in% genes_over_thr ~ "over_thr",
          TRUE ~ "other"
        )
      )
    
    colors <- c(
      special = "#29ADB2",
      top10_neg = "#65B741",
      top10_pos = "#C51605",
      over_thr = "#637A9F",
      other = "#B0B0B0"
    )
    
    sizes <- c(
      special = 15,
      top10_neg = 12,
      top10_pos = 12,
      over_thr = 12,
      other = 8
    )
    
    df_joined$color <- colors[df_joined$gene_type]
    df_joined$size <- sizes[df_joined$gene_type]
    df_joined$color[is.na(df_joined$color)] <- "#B0B0B0"
    df_joined$size[is.na(df_joined$size)] <- 8
    
    plot_ly(
      df_joined,
      x = ~pq, y = ~poso,
      type = 'scatter', mode = 'markers',
      text = ~`Gene names`, hoverinfo = 'text',
      marker = list(
        color = df_joined$color,
        size = df_joined$size,
        line = list(width = 1, color = '#333333')
      ),
      source = "scatter1"
    ) %>%
      layout(
        title = "OPLS-DA Loading Plot",
        plot_bgcolor = "#EDF4F4",
        paper_bgcolor = "#EDF4F4",
        font = list(color = '#4F0433')
      )
  })
  
  # ----- Scatter 2: Transformed（color decided by threshold） -----
  output$scatter2 <- renderPlotly({
    df <- transformed_data(); req(df)
    colnames(df)[1] <- "Gene names"
    
    opls_df <- opls_data(); req(opls_df)
    colnames(opls_df)[1:3] <- c("Gene names", "pq", "poso")
    
    value_col <- input$drugname
    trans_sub <- df %>% dplyr::select(`Gene names`, dplyr::all_of(value_col))
    
    opls_joined <- opls_df %>% dplyr::left_join(trans_sub, by = "Gene names")
    
    # 用阈值而不是 4.5
    genes_over_thr <- opls_joined %>%
      dplyr::filter(!is.na(.data[[value_col]]), .data[[value_col]] > input$threshold) %>%
      dplyr::pull(`Gene names`)
    
    top10_genes <- opls_joined %>%
      dplyr::filter(!is.na(.data[[value_col]])) %>%
      dplyr::arrange(dplyr::desc(.data[[value_col]])) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::pull(`Gene names`)
    
    opls_joined <- opls_joined %>%
      dplyr::mutate(
        gene_type = dplyr::case_when(
          `Gene names` %in% top10_genes & pq < 0 ~ "top10_neg",
          `Gene names` %in% top10_genes & pq >= 0 ~ "top10_pos",
          `Gene names` %in% genes_over_thr ~ "over_thr",
          TRUE ~ "other"
        )
      )
    
    colors <- c(
      top10_neg = "#65B741",
      top10_pos = "#C51605",
      over_thr = "#637A9F",
      other = "#B0B0B0"
    )
    
    sizes <- c(
      top10_neg = 12,
      top10_pos = 12,
      over_thr = 8,
      other = 6
    )
    
    opls_joined$color <- colors[opls_joined$gene_type]
    opls_joined$size <- sizes[opls_joined$gene_type]
    
    df_joined <- df %>%
      dplyr::left_join(dplyr::select(opls_joined, `Gene names`, color, size), by = "Gene names")
    
    df_joined$color[is.na(df_joined$color)] <- "#B0B0B0"
    df_joined$size[is.na(df_joined$size)] <- 6
    
    plot_ly(
      df_joined,
      x = ~get(value_col),
      y = ~y,
      type = 'scatter',
      mode = 'markers',
      text = ~`Gene names`,
      hoverinfo = 'text',
      marker = list(
        color = df_joined$color,
        size = df_joined$size,
        line = list(width = 1, color = '#333333')
      ),
      source = "scatter2"
    ) %>%
      layout(
        title = paste("Transformed Data –", input$drugname),
        xaxis = list(title = "Transformed OPLS x-coordinates"),
        plot_bgcolor = "#EDF4F4",
        paper_bgcolor = "#EDF4F4",
        font = list(color = '#4F0433')
      )
  })
  
  # proteins for enrichment（significant proteins !=0 is >= threshold）
  selected_genes <- reactive({
    df <- transformed_data(); req(df)
    value_col <- input$drugname
    df %>%
      filter(.data[[value_col]] != 0) %>%
      pull(`Gene names`) %>%
      unique()
  })
  
  go_result_bp <- reactiveVal(NULL)
  go_result_cc <- reactiveVal(NULL)
  go_result_mf <- reactiveVal(NULL)
  
  # Listen to changes in selected_genes() and trigger asynchronous computation
  observeEvent(selected_genes(), {
    genes <- selected_genes()
    if (length(genes) < 10) {
      go_result_bp(NULL); go_result_cc(NULL); go_result_mf(NULL)
      return()
    }
    future({
      enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    }) %...>% go_result_bp()
    future({
      enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    }) %...>% go_result_cc()
    future({
      enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    }) %...>% go_result_mf()
  })
  
  go_ready <- reactive({
    !is.null(go_result_bp()) && !is.null(go_result_cc()) && !is.null(go_result_mf())
  })
  
  output$go_status <- renderUI({
    if (!go_ready()) {
      tags$div(style = "color: #777; font-style: italic;", "GO enrichment is running, please wait...")
    } else {
      NULL
    }
  })
  
  output$plot_go_bp <- renderPlot({
    req(go_ready())
    res <- go_result_bp()
    if (nrow(res) == 0) { plot.new(); title("No enriched GO Biological Process terms") }
    else { dotplot(res, showCategory = 10, title = "GO Biological Process") }
  })
  
  output$plot_go_cc <- renderPlot({
    req(go_ready())
    res <- go_result_cc()
    if (nrow(res) == 0) { plot.new(); title("No enriched GO Cellular Component terms") }
    else { dotplot(res, showCategory = 10, title = "GO Cellular Component") }
  })
  
  output$plot_go_mf <- renderPlot({
    req(go_ready())
    res <- go_result_mf()
    if (nrow(res) == 0) { plot.new(); title("No enriched GO Molecular Function terms") }
    else { dotplot(res, showCategory = 10, title = "GO Molecular Function") }
  })
  
  selected_gene <- reactiveVal(NULL)
  
  observeEvent(event_data("plotly_click", source = "scatter1"), {
    df <- opls_data()
    click <- event_data("plotly_click", source = "scatter1")
    if (!is.null(click)) selected_gene(df$`Gene names`[click$pointNumber + 1])
  })
  
  observeEvent(event_data("plotly_click", source = "scatter2"), {
    df <- transformed_data()
    click <- event_data("plotly_click", source = "scatter2")
    if (!is.null(click)) selected_gene(df$`Gene names`[click$pointNumber + 1])
  })
  
  drug_gene_data <- reactive({
    req(input$drugname, selected_gene())
    file_path <- file.path("data", "drugs", paste0(input$drugname, ".xlsx"))
    if (!file.exists(file_path)) return(NULL)
    df <- read_excel(file_path)
    df %>% filter(`Gene names` == selected_gene())
  })
  
  # Bar chart（bar color and line follows threshold selection）
  output$bar_chart <- renderPlotly({
    df <- drug_gene_data(); req(df)
    
    values <- as.numeric(as.character(df[1, -1]))
    names <- colnames(df)[-1]
    valid_idx <- !is.na(values)
    values <- values[valid_idx]
    names <- names[valid_idx]
    
    colors <- ifelse(values >= input$threshold, "#637A9F", "#B0B0B0")
    colors[is.na(colors)] <- "#B0B0B0"
    
    plot_ly(
      x = names, y = values,
      type = "bar",
      marker = list(color = colors)
    ) %>%
      layout(
        title = paste("Bar chart for", selected_gene()),
        plot_bgcolor = "#EDF4F4",
        paper_bgcolor = "#EDF4F4",
        yaxis = list(
          title = "Transformed OPLS-DA x coordinates",
          zeroline = FALSE, showline = TRUE, linecolor = '#333',
          rangemode = "tozero", zerolinecolor = "#999",
          zerolinewidth = 2, gridcolor = "#ccc"
        ),
        xaxis = list(tickangle = -45),
        shapes = list(
          list(
            type = "line", x0 = 0, x1 = 1, xref = "paper",
            y0 = input$threshold, y1 = input$threshold,
            line = list(color = "red", width = 2, dash = "dash")
          )
        )
      )
  })
  
  # ---- UpSet & Pro-target（threshold linkage + dynamic column name） ----
  drug_data <- reactive({
    req(input$drugname)
    file_path <- file.path("data", "drugs", paste0(input$drugname, ".xlsx"))
    if (!file.exists(file_path)) return(NULL)
    readxl::read_excel(file_path)
  })
  
  upset_data <- reactive({
    df <- drug_data(); req(df)
    if (ncol(df) < 5) return(NULL)
    
    group_cols <- colnames(df)[2:5]
    
    df_raw <- df %>%
      rename_with(~ gsub(" ", "_", .), all_of(group_cols)) %>%
      mutate(Gene_name = df[[1]])
    
    group_cols_new <- colnames(df_raw)[2:5]
    
    # threshold linkage
    count_col <- paste0("count_above_", format(round(input$threshold, 1), nsmall = 1))
    
    df_raw[[count_col]] <- rowSums(df_raw[, group_cols_new] > input$threshold, na.rm = TRUE)
    
    # same logic and threshold for UpSet
    df_logic <- df_raw %>%
      mutate(across(all_of(group_cols_new), ~ . > input$threshold)) %>%
      mutate(across(all_of(group_cols_new), ~ tidyr::replace_na(., FALSE)))
    
    list(
      raw = df_raw,       # for table
      logic = df_logic,   # for upset plot
      groups = group_cols_new,
      count_col = count_col
    )
  })
  
  output$upset_plot <- renderPlot({
    res <- upset_data(); req(res)
    df_logic <- res$logic
    groups <- res$groups
    df_logic <- df_logic[rowSums(df_logic[groups]) > 0, ]
    p <- ComplexUpset::upset(
      df_logic,
      intersect = groups,   
      name = "Condition",
      matrix = intersection_matrix(geom = geom_point(size = 5)),
      queries = list(
        upset_query(set = "A549_cell_lysate", color = "#870052", fill = "#870052"),
        upset_query(set = "A549_intact_cell", color = "#4DB5BC", fill = "#4DB5BC"),
        upset_query(set = "H82_cell_lysate", color = "#54B986", fill = "#54B986"),
        upset_query(set = "H82_intact_cell", color = "#FFC66D", fill = "#FFC66D")
      )
    ) + 
      labs(x = "combination") +
      theme(axis.text.y = element_text(size = 14))
    p
  })
  
  output$protarget_table <- DT::renderDT({
    res <- upset_data(); req(res)
    df_raw <- res$raw
    count_col <- res$count_col
    
    df_filtered <- df_raw %>%
      dplyr::select(-Gene_name) %>%
      dplyr::filter(.data[[count_col]] > 1)
    
    DT::datatable(df_filtered, options = list(pageLength = 10))
  })
  
  output$download_protarget <- downloadHandler(
    filename = function() {
      paste0(input$drugname, "_Protarget.xlsx")
    },
    content = function(file) {
      res <- upset_data()
      if (!is.null(res)) {
        df_raw <- res$raw
        count_col <- res$count_col
        df_filtered <- df_raw %>%
          dplyr::select(-Gene_name) %>%
          dplyr::filter(.data[[count_col]] > 1)
        writexl::write_xlsx(df_filtered, path = file)
      }
    }
  )
}

shinyApp(ui, server)

#required files
#renv::snapshot()


