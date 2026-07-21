#Load Packages
library(shiny)
library(shinyFiles)
library(data.table)
library(damr)
library(bslib)
library(sleepr)
library(colourpicker)
library(scales)
library(survival)
library(survminer) 
library(ggplot2)
library(multcomp)

ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),
  
  tags$style(HTML("
    .center-wrap { max-width: 860px; margin: 0 auto; }
    .cardish {
      background: rgba(255,255,255,.96);
      border: 1px solid rgba(0,0,0,.08);
      border-radius: 14px;
      padding: 16px;
      box-shadow: 0 6px 18px rgba(0,0,0,.06);
      margin-bottom: 14px;
    }
    .small-help { font-size: 0.95rem; opacity: 0.85; }

    .filebox {
      max-height: 260px;
      overflow-y: auto;
      background: #f8fafc;
      color: #0f172a;
      border: 1px solid rgba(0,0,0,.10);
      border-radius: 12px;
      padding: 10px 12px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;
      font-size: 0.85rem;
      white-space: pre;
    }

    .kpi {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
      margin-bottom: 8px;
    }
    .kpiCard {
      border-radius: 14px;
      padding: 12px 14px;
      border: 1px solid rgba(0,0,0,.08);
      background: rgba(255,255,255,.96);
    }
    .kpiLabel { font-size: 0.85rem; opacity: .75; margin-bottom: 4px; }
    .kpiValue { font-size: 1.6rem; font-weight: 800; }

    .btn { border-radius: 12px; }
    .btn-primary { border-radius: 12px; }

    .errorBox {
      border: 1px solid rgba(220,38,38,.35);
      background: rgba(220,38,38,.08);
      color: #7f1d1d;
      border-radius: 14px;
      padding: 12px 14px;
      margin-top: 12px;
    }
    .infoBox {
      border: 1px solid rgba(2,132,199,.25);
      background: rgba(2,132,199,.08);
      color: #0c4a6e;
      border-radius: 14px;
      padding: 12px 14px;
      margin-top: 12px;
    }
  ")),
  
  uiOutput("main_ui")
)

server <- function(input, output, session) {
  
  # ============================================================================
  # Plot helper (NO survminer::surv_summary(); avoids NA/Inf terminal CI rows)
  # ============================================================================
  
  make_palette <- function(levels_chr) {
    levels_chr <- as.character(levels_chr)
    cols <- scales::hue_pal()(length(levels_chr))
    setNames(cols, levels_chr)
  }
  
  app_primary_colour <- "#1f77b4"
  
  custom_theme <- function(base_size = 14) {
    theme_classic(base_size = base_size) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size * 0.85),
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = base_size * 0.85),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
      )
  }
  
  km_plot_df <- function(fit, title = NULL, palette = NULL) {
    
    s <- summary(fit)
    
    df <- data.frame(
      time  = s$time,
      surv  = s$surv,
      lower = s$lower,
      upper = s$upper,
      strata = if (!is.null(s$strata)) as.character(s$strata) else NA_character_
    )
    
    if (!all(is.na(df$strata))) {
      df$strata <- gsub("^.*=", "", df$strata)
    }
    
    # robust filtering (prevents "blank" plots due to NA/Inf CI rows)
    rib <- df[is.finite(df$lower) & is.finite(df$upper) & is.finite(df$time), , drop = FALSE]
    stp <- df[is.finite(df$surv)  & is.finite(df$time), , drop = FALSE]
    
    # overall (no strata)
    if (all(is.na(stp$strata))) {
      p <- ggplot() +
        geom_step(
          data = stp,
          aes(x = time, y = surv),
          linewidth = 1.1,
          colour = app_primary_colour
        )
      
      if (nrow(rib) > 0) {
        p <- p +
          geom_ribbon(
            data = rib,
            aes(x = time, ymin = lower, ymax = upper),
            alpha = 0.15,
            linewidth = 0,
            fill = app_primary_colour
          )
      }
      
      return(
        p +
          labs(title = title, x = "Time (Days)", y = "Survival Probability") +
          custom_theme(base_size = 12)
      )
    }
    
    # stratified
    p <- ggplot() +
      geom_step(
        data = stp,
        aes(x = time, y = surv, colour = strata, group = strata),
        linewidth = 1.1
      )
    
    if (nrow(rib) > 0) {
      p <- p +
        geom_ribbon(
          data = rib,
          aes(x = time, ymin = lower, ymax = upper, fill = strata, group = strata),
          alpha = 0.15,
          linewidth = 0
        )
    }
    
    p <- p +
      labs(
        title = title,
        x = "Time (Days)",
        y = "Survival Probability",
        colour = "Group",
        fill   = "Group"
      ) +
      custom_theme(base_size = 12)
    
    if (!is.null(palette) && length(palette) > 0) {
      p <- p +
        scale_color_manual(values = palette, drop = FALSE) +
        scale_fill_manual(values = palette, drop = FALSE)
    }
    
    p
  }
  
  # ============================================================================
  # Page controller
  # ============================================================================
  page <- reactiveVal("welcome")
  observeEvent(input$go_to_load,       { page("load") },       ignoreInit = TRUE)
  observeEvent(input$back_to_welcome,  { page("welcome") },    ignoreInit = TRUE)
  observeEvent(input$back_to_load,     { page("load") },       ignoreInit = TRUE)
  observeEvent(input$go_to_lifespan,   { page("lifespan") },   ignoreInit = TRUE)
  observeEvent(input$back_to_curation, { page("curation") },   ignoreInit = TRUE)
  
  # ============================================================================
  # Directory chooser
  # ============================================================================
  # Use absolute, normalized roots on every platform. On Windows, also expose
  # all currently mounted drive letters (C:, D:, network-mounted drives, etc.).
  home_dir <- normalizePath(
    path.expand("~"),
    winslash = "/",
    mustWork = FALSE
  )
  
  roots <- c(Home = home_dir)
  if (.Platform$OS.type == "windows") {
    windows_volumes <- shinyFiles::getVolumes()()
    roots <- c(roots, windows_volumes)
    
    # Avoid showing the same location twice if Home is a drive root.
    root_keys <- tolower(gsub("[\\\\/]+$", "", unname(roots)))
    roots <- roots[!duplicated(root_keys)]
  }
  
  shinyDirChoose(input, "data_dir", roots = roots, session = session)
  
  data_dir <- reactive({
    req(input$data_dir)
    selected_dir <- parseDirPath(roots, input$data_dir)
    
    validate(need(
      length(selected_dir) == 1L && dir.exists(selected_dir),
      "Please choose an existing data folder."
    ))
    
    # Forward slashes are accepted by R on Windows and avoid backslash escaping
    # problems while remaining valid on macOS and Linux.
    normalizePath(selected_dir, winslash = "/", mustWork = TRUE)
  })
  
  output$files_preview <- renderText({
    req(data_dir())
    paste(list.files(data_dir()), collapse = "\n")
  })
  
  results_rv <- reactiveVal(NULL)
  err_rv <- reactiveVal(NULL)
  
  # ============================================================================
  # Run curation
  # ============================================================================
  observeEvent(input$run_step2, {
    req(data_dir())
    err_rv(NULL)
    
    res <- tryCatch({
      
      if (isTRUE(input$set_wd)) setwd(data_dir())
      
      meta_path <- file.path(data_dir(), input$metadata_name)
      validate(need(file.exists(meta_path),
                    paste0("Could not find metadata file: ", meta_path)))
      
      metadata <- fread(meta_path)
      validate(need("start_datetime" %in% names(metadata),
                    "metadata must contain a 'start_datetime' column."))
      
      metadata[, batch := as.character(start_datetime)]
      metadata <- link_dam_metadata(metadata, result_dir = data_dir())
      
      dt <- load_dam(metadata)
      
      # curate dead animals
      dt_curated <- sleepr::curate_dead_animals(dt, moving_var = activity)
      
      all_ids <- dt[, id, meta = TRUE]
      kept_ids <- dt_curated[, id, meta = TRUE]
      removed_ids <- setdiff(all_ids, kept_ids)
      
      list(
        metadata = metadata,
        dt = dt,
        dt_curated = dt_curated,
        removed_ids = removed_ids
      )
      
    }, error = function(e) {
      err_rv(conditionMessage(e))
      NULL
    })
    
    if (!is.null(res)) {
      results_rv(res)
      page("curation")
    }
  }, ignoreInit = TRUE)
  
  output$load_error <- renderUI({
    msg <- err_rv()
    if (is.null(msg)) return(NULL)
    div(class = "errorBox", tags$b("ERROR! "), tags$span(msg))
  })
  
  # ============================================================================
  # KPI boxes
  # ============================================================================
  output$curation_kpis <- renderUI({
    res <- results_rv()
    if (is.null(res)) {
      return(div(class = "infoBox",
                 tags$b("No results yet. "),
                 "Go back to the data loading page and click Run curation."))
    }
    
    total <- length(res$dt[, id, meta = TRUE])
    kept  <- length(res$dt_curated[, id, meta = TRUE])
    removed <- length(res$removed_ids)
    
    div(class = "kpi",
        div(class = "kpiCard",
            div(class = "kpiLabel", "Total animals"),
            div(class = "kpiValue", total)
        ),
        div(class = "kpiCard",
            div(class = "kpiLabel", "Removed (dead/excluded)"),
            div(class = "kpiValue", removed)
        ),
        div(class = "kpiCard",
            div(class = "kpiLabel", "Kept for analysis"),
            div(class = "kpiValue", kept)
        )
    )
  })
  
  # ============================================================================
  # Multifactor lifespan analysis (metadata-driven)
  # ============================================================================
  meta_cols_r <- reactive({
    res <- results_rv()
    req(res)
    md <- copy(res$metadata)
    if ("file_info" %in% names(md)) md[, file_info := NULL]
    
    exclude <- c("id", "region_id", "experiment_id",
                 "start_datetime", "stop_datetime", "batch")
    
    setdiff(names(md), exclude)
  })
  
  output$groupvar_selector <- renderUI({
    choices <- meta_cols_r()
    if (length(choices) == 0) {
      return(div(class = "infoBox",
                 tags$b("No grouping variables found. "),
                 "Check your metadata columns."))
    }
    
    default_sel <- intersect(c("treatment", "sex"), choices)
    if (length(default_sel) == 0) default_sel <- choices[1]
    
    selectizeInput(
      "group_vars",
      "Choose grouping variables (e.g., treatment, sex)",
      choices = choices,
      selected = default_sel,
      multiple = TRUE
    )
  })
  
  lifespan_dt_r <- reactive({
    res <- results_rv()
    req(res)
    
    dt_lifespan <- res$dt_curated
    md <- copy(res$metadata)
    
    lifespan_by_id <- dt_lifespan[, .(lifespan = max(t)), by = id]
    lifespan_by_id[, deathdays := lifespan / (24 * 60 * 60)]
    
    md <- unique(md, by = "id")
    if ("file_info" %in% names(md)) md[, file_info := NULL]
    
    merge(lifespan_by_id, md, by = "id", all.x = TRUE)
  })
  
  output$level_selectors <- renderUI({
    df <- lifespan_dt_r()
    req(df)
    gvars <- input$group_vars
    req(gvars)
    
    ui_list <- lapply(gvars, function(g) {
      if (!g %in% names(df)) return(NULL)
      
      levs <- sort(unique(as.character(df[[g]])))
      levs <- levs[!is.na(levs) & nzchar(levs)]
      
      tagList(
        tags$hr(),
        tags$h5(sprintf("Filter levels: %s", g)),
        selectizeInput(
          inputId = paste0("levels__", g),
          label = "Levels to include",
          choices = levs,
          selected = levs,
          multiple = TRUE
        ),
        selectInput(
          inputId = paste0("ref__", g),
          label = "Reference level (for Cox model)",
          choices = c("(no change)" = "", levs),
          selected = if (g == "treatment" && "Control" %in% levs) "Control" else ""
        )
      )
    })
    
    tagList(
      div(class = "infoBox",
          tags$b("Tip: "),
          "Use these selectors to include/exclude levels for each factor. Reference level affects Cox model interpretation."
      ),
      ui_list
    )
  })
  
  lifespan_filt_r <- reactive({
    df <- lifespan_dt_r()
    req(df)
    gvars <- input$group_vars
    req(gvars)
    
    out <- copy(df)
    
    for (g in gvars) {
      lev_id <- paste0("levels__", g)
      if (g %in% names(out) && !is.null(input[[lev_id]])) {
        out <- out[get(g) %in% input[[lev_id]]]
      }
    }
    
    out <- out[!is.na(deathdays)]
    out <- out[is.finite(deathdays)]
    
    for (g in gvars) {
      if (g %in% names(out)) out[, (g) := droplevels(as.factor(get(g)))]
    }
    
    for (g in gvars) {
      ref_id <- paste0("ref__", g)
      if (g %in% names(out)) {
        ref <- input[[ref_id]]
        if (!is.null(ref) && nzchar(ref) && ref %in% levels(out[[g]])) {
          out[, (g) := relevel(get(g), ref = ref)]
        }
      }
    }
    
    out
  })
  
  # ============================================================================
  # Overall KM
  # ============================================================================
  output$km_overall <- renderPlot({
    df <- lifespan_filt_r()
    validate(need(nrow(df) > 0, "No animals available after filtering."))
    
    dff <- as.data.frame(df)
    fit <- survfit(Surv(deathdays) ~ 1, data = dff)
    
    km_plot_df(fit, title = "Kaplan–Meier (overall)")
  })
  
  # ============================================================================
  # Choose one factor for single grouped display
  # ============================================================================
  output$which_factor_plot_ui <- renderUI({
    req(input$group_vars)
    selectInput(
      "plot_factor",
      "Factor for grouped KM plot",
      choices = input$group_vars,
      selected = input$group_vars[[1]]
    )
  })
  
  output$colour_selectors <- renderUI({
    df <- lifespan_filt_r()
    req(input$plot_factor)
    g <- as.character(input$plot_factor)
    
    dff <- as.data.frame(df)
    validate(need(g %in% names(dff), "Selected factor not in data."))
    
    levs <- sort(unique(as.character(dff[[g]])))
    levs <- levs[!is.na(levs) & nzchar(levs)]
    validate(need(length(levs) >= 1, "No levels available for colour selection."))
    
    default_cols <- setNames(scales::hue_pal()(length(levs)), levs)
    
    tagList(
      div(class = "infoBox",
          tags$b("Colours: "),
          "Pick a colour for each level of the selected factor (applies to the grouped KM plot below)."
      ),
      lapply(levs, function(lv) {
        colourInput(
          inputId = paste0("col__", g, "__", lv),
          label   = lv,
          value   = default_cols[[lv]],
          showColour = "background"
        )
      })
    )
  })
  
  output$km_by_factor <- renderPlot({
    df <- lifespan_filt_r()
    req(input$plot_factor)
    g <- as.character(input$plot_factor)
    
    dff <- as.data.frame(df)
    validate(need(nrow(dff) > 0, "No animals available after filtering."))
    validate(need(g %in% names(dff), "Selected factor not in data."))
    
    dff <- dff[!is.na(dff[[g]]), , drop = FALSE]
    dff[[g]] <- droplevels(as.factor(dff[[g]]))
    validate(need(nlevels(dff[[g]]) >= 1, "No levels available for the selected factor."))
    
    fml <- as.formula(sprintf("Surv(deathdays) ~ `%s`", g))
    fit <- survfit(fml, data = dff)
    
    levs <- levels(dff[[g]])
    pal <- setNames(rep(NA_character_, length(levs)), levs)
    for (lv in levs) {
      id <- paste0("col__", g, "__", lv)
      v  <- input[[id]]
      if (!is.null(v) && nzchar(v)) pal[[lv]] <- v
    }
    pal <- pal[!is.na(pal)]
    
    km_plot_df(fit, title = paste("Kaplan–Meier by", g), palette = pal)
  })
  
  output$cox_pairwise <- renderTable({
    df <- lifespan_filt_r()
    req(input$plot_factor)
    g <- as.character(input$plot_factor)
    
    dff <- as.data.frame(df)
    validate(need(nrow(dff) > 0, "No animals available after filtering."))
    validate(need(g %in% names(dff), "Selected factor not in data."))
    
    dff <- dff[!is.na(dff[[g]]), , drop = FALSE]
    dff[[g]] <- droplevels(as.factor(dff[[g]]))
    validate(need(nlevels(dff[[g]]) >= 2, "Need at least 2 levels for Cox model comparisons."))
    
    fml <- as.formula(sprintf("Surv(deathdays) ~ `%s`", g))
    cox_model <- coxph(fml, data = dff)
    
    if (nlevels(dff[[g]]) <= 2) {
      s <- summary(cox_model)
      return(data.table(
        term = rownames(s$coefficients),
        coef = s$coefficients[, "coef"],
        HR = s$coefficients[, "exp(coef)"],
        se = s$coefficients[, "se(coef)"],
        z = s$coefficients[, "z"],
        p = s$coefficients[, "Pr(>|z|)"]
      ))
    }
    
    pw <- glht(cox_model, linfct = do.call(mcp, setNames(list("Tukey"), g)))
    s <- summary(pw)
    
    data.table(
      comparison = names(s$test$coefficients),
      estimate = as.numeric(s$test$coefficients),
      se = as.numeric(s$test$sigma),
      z = as.numeric(s$test$tstat),
      p = as.numeric(s$test$pvalues)
    )
  }, striped = TRUE, hover = TRUE, spacing = "s")
  
  # ============================================================================
  # Auto KM + Cox stats for ALL selected factors
  # ============================================================================
  output$all_factor_panels <- renderUI({
    req(input$group_vars)
    gvars <- input$group_vars
    validate(need(length(gvars) >= 1, "Select at least one grouping variable."))
    
    tagList(lapply(gvars, function(g) {
      tagList(
        tags$hr(),
        tags$h4(sprintf("Factor: %s", g)),
        plotOutput(outputId = paste0("km_all__", g), height = "360px"),
        tags$h5("Cox model / pairwise comparisons"),
        tableOutput(outputId = paste0("cox_all__", g))
      )
    }))
  })
  
  observeEvent(input$group_vars, {
    gvars <- input$group_vars
    if (is.null(gvars) || length(gvars) == 0) return(NULL)
    
    for (g in gvars) {
      local({
        gg <- as.character(g)
        
        output[[paste0("km_all__", gg)]] <- renderPlot({
          df <- lifespan_filt_r()
          dff <- as.data.frame(df)
          validate(need(nrow(dff) > 0, "No animals available after filtering."))
          validate(need(gg %in% names(dff), sprintf("Factor '%s' not in data.", gg)))
          
          dff <- dff[!is.na(dff[[gg]]), , drop = FALSE]
          dff[[gg]] <- droplevels(as.factor(dff[[gg]]))
          validate(need(nlevels(dff[[gg]]) >= 1, sprintf("No levels available for '%s'.", gg)))
          
          fml <- as.formula(sprintf("Surv(deathdays) ~ `%s`", gg))
          fit <- survfit(fml, data = dff)
          
          pal <- make_palette(levels(dff[[gg]]))
          km_plot_df(fit, title = paste("Kaplan–Meier by", gg), palette = pal)
        })
        
        output[[paste0("cox_all__", gg)]] <- renderTable({
          df <- lifespan_filt_r()
          dff <- as.data.frame(df)
          validate(need(nrow(dff) > 0, "No animals available after filtering."))
          validate(need(gg %in% names(dff), sprintf("Factor '%s' not in data.", gg)))
          
          dff <- dff[!is.na(dff[[gg]]), , drop = FALSE]
          dff[[gg]] <- droplevels(as.factor(dff[[gg]]))
          validate(need(nlevels(dff[[gg]]) >= 2, sprintf("Need >=2 levels for '%s'.", gg)))
          
          fml <- as.formula(sprintf("Surv(deathdays) ~ `%s`", gg))
          cox_model <- coxph(fml, data = dff)
          
          if (nlevels(dff[[gg]]) <= 2) {
            s <- summary(cox_model)
            return(data.table(
              term = rownames(s$coefficients),
              coef = s$coefficients[, "coef"],
              HR = s$coefficients[, "exp(coef)"],
              se = s$coefficients[, "se(coef)"],
              z = s$coefficients[, "z"],
              p = s$coefficients[, "Pr(>|z|)"]
            ))
          }
          
          pw <- glht(cox_model, linfct = do.call(mcp, setNames(list("Tukey"), gg)))
          s <- summary(pw)
          
          data.table(
            comparison = names(s$test$coefficients),
            estimate = as.numeric(s$test$coefficients),
            se = as.numeric(s$test$sigma),
            z = as.numeric(s$test$tstat),
            p = as.numeric(s$test$pvalues)
          )
        }, striped = TRUE, hover = TRUE, spacing = "s")
        
      })
    }
  }, ignoreInit = TRUE)
  
  # ============================================================================
  # Interaction section (choose any two metadata columns)
  # ============================================================================
  output$interaction_selector <- renderUI({
    df <- lifespan_dt_r()
    req(df)
    
    exclude <- c("id", "lifespan", "deathdays", "start_datetime", "stop_datetime", "batch")
    choices <- setdiff(names(df), exclude)
    validate(need(length(choices) >= 2, "Not enough metadata columns available for interaction."))
    
    default_a <- if ("treatment" %in% choices) "treatment" else choices[[1]]
    default_b <- if ("sex" %in% choices) "sex" else choices[[2]]
    
    tagList(
      div(class = "infoBox",
          tags$b("Interaction: "),
          "Choose two factors. The app will plot combined groups (A:B) and test A*B in a Cox model."
      ),
      fluidRow(
        column(6, selectInput("int_a", "Factor A", choices = choices, selected = default_a)),
        column(6, selectInput("int_b", "Factor B", choices = choices, selected = default_b))
      )
    )
  })
  
  output$km_interaction <- renderPlot({
    df <- lifespan_filt_r()
    req(input$int_a, input$int_b)
    
    a <- as.character(input$int_a)
    b <- as.character(input$int_b)
    validate(need(a != b, "Choose two different factors."))
    
    dff <- as.data.frame(df)
    validate(need(a %in% names(dff) && b %in% names(dff), "Interaction factors not in data."))
    
    dff <- dff[!is.na(dff[[a]]) & !is.na(dff[[b]]), , drop = FALSE]
    validate(need(nrow(dff) > 0, "No rows for interaction after filtering."))
    
    dff[[a]] <- droplevels(as.factor(dff[[a]]))
    dff[[b]] <- droplevels(as.factor(dff[[b]]))
    dff$int_group <- interaction(dff[[a]], dff[[b]], sep = " : ", drop = TRUE)
    
    validate(need(nlevels(dff$int_group) >= 2, "Not enough interaction groups after filtering."))
    
    fit <- survfit(Surv(deathdays) ~ int_group, data = dff)
    pal <- make_palette(levels(dff$int_group))
    
    km_plot_df(fit, title = paste("Kaplan–Meier:", a, "×", b), palette = pal)
  })
  
  output$cox_interaction_text <- renderPrint({
    df <- lifespan_filt_r()
    req(input$int_a, input$int_b)
    
    a <- as.character(input$int_a)
    b <- as.character(input$int_b)
    validate(need(a != b, "Choose two different factors."))
    
    dff <- as.data.frame(df)
    validate(need(a %in% names(dff) && b %in% names(dff), "Interaction factors not in data."))
    
    dff <- dff[!is.na(dff[[a]]) & !is.na(dff[[b]]), , drop = FALSE]
    validate(need(nrow(dff) > 0, "No rows for interaction after filtering."))
    
    dff[[a]] <- droplevels(as.factor(dff[[a]]))
    dff[[b]] <- droplevels(as.factor(dff[[b]]))
    
    validate(need(nlevels(dff[[a]]) >= 2 && nlevels(dff[[b]]) >= 2,
                  "Need at least 2 levels in each factor for an interaction test."))
    
    f_full <- as.formula(sprintf("Surv(deathdays) ~ `%s` * `%s`", a, b))
    f_main <- as.formula(sprintf("Surv(deathdays) ~ `%s` + `%s`", a, b))
    
    cox_full <- coxph(f_full, data = dff)
    cox_main <- coxph(f_main, data = dff)
    
    cat("Cox model with interaction:\n")
    print(summary(cox_full))
    
    cat("\nLikelihood ratio test for adding interaction (main effects vs interaction):\n")
    print(anova(cox_main, cox_full, test = "LRT"))
  })
  
  output$cox_interaction_pairwise <- renderTable({
    df <- lifespan_filt_r()
    req(input$int_a, input$int_b)
    
    a <- as.character(input$int_a)
    b <- as.character(input$int_b)
    validate(need(a != b, "Choose two different factors."))
    
    dff <- as.data.frame(df)
    validate(need(a %in% names(dff) && b %in% names(dff), "Interaction factors not in data."))
    
    dff <- dff[!is.na(dff[[a]]) & !is.na(dff[[b]]), , drop = FALSE]
    validate(need(nrow(dff) > 0, "No rows for interaction after filtering."))
    
    dff[[a]] <- droplevels(as.factor(dff[[a]]))
    dff[[b]] <- droplevels(as.factor(dff[[b]]))
    dff$int_group <- interaction(dff[[a]], dff[[b]], sep = " : ", drop = TRUE)
    
    validate(need(nlevels(dff$int_group) >= 2, "Not enough interaction groups for pairwise tests."))
    
    cox_grp <- coxph(Surv(deathdays) ~ int_group, data = dff)
    
    if (nlevels(dff$int_group) <= 2) {
      s <- summary(cox_grp)
      return(data.table(
        term = rownames(s$coefficients),
        coef = s$coefficients[, "coef"],
        HR = s$coefficients[, "exp(coef)"],
        se = s$coefficients[, "se(coef)"],
        z = s$coefficients[, "z"],
        p = s$coefficients[, "Pr(>|z|)"]
      ))
    }
    
    pw <- glht(cox_grp, linfct = mcp(int_group = "Tukey"))
    s <- summary(pw)
    
    data.table(
      comparison = names(s$test$coefficients),
      estimate = as.numeric(s$test$coefficients),
      se = as.numeric(s$test$sigma),
      z = as.numeric(s$test$tstat),
      p = as.numeric(s$test$pvalues)
    )
  }, striped = TRUE, hover = TRUE, spacing = "s")
  
  # ============================================================================
  # Main UI pages
  # ============================================================================
  output$main_ui <- renderUI({
    
    if (page() == "welcome") {
      
      fluidPage(
        br(), br(),
        div(class = "center-wrap",
            div(class = "cardish",
                h2("Welcome to Lab 519's Trikinetics App"),
                p("For more detailed instructions on how to use this app, see the manual at the link below:"),
                tags$a(
                  href = "https://github.com/C-L-Thomas/Nasonia_Trikinetics/blob/main/Manual__Interactive_Trikinetics_App.pdf",
                  target = "_blank",
                  "Open the Interactive Trikinetics App Manual (PDF)"
                ),
                br(), br(),
                p("To proceed, click Next to go to the data loading page."),
                actionButton("go_to_load", "Next", class = "btn-primary btn-lg")
            )
        )
      )
      
    } else if (page() == "load") {
      
      fluidPage(
        br(),
        div(class = "center-wrap",
            div(class = "cardish",
                h3("Step 1: Load data"),
                p(class = "small-help",
                  "Select the folder that contains your metadata CSV and DAM results, then click Run curation."
                ),
                shinyDirButton("data_dir", "Choose data folder", "Browse…"),
                checkboxInput("set_wd", "Set as working directory (setwd)", value = TRUE),
                textInput("metadata_name", "Metadata filename", value = "metadata.csv"),
                br(),
                actionButton("run_step2", "Run curation", class = "btn-primary"),
                uiOutput("load_error"),
                br(),
                actionButton("back_to_welcome", "← Back", class = "btn-secondary")
            ),
            
            div(class = "cardish",
                h4("Files in selected folder"),
                div(class = "filebox", verbatimTextOutput("files_preview"))
            )
        )
      )
      
    } else if (page() == "curation") {
      
      fluidPage(
        br(),
        div(class = "center-wrap",
            div(class = "cardish",
                h3("Data Curation / Removal of Dead"),
                uiOutput("curation_kpis"),
                br(),
                actionButton("go_to_lifespan", "Next: Lifespan analysis →", class = "btn-primary"),
                br(), br(),
                actionButton("back_to_load", "← Back to data loading", class = "btn-secondary")
            )
        )
      )
      
    } else {
      
      fluidPage(
        br(),
        div(class = "center-wrap",
            div(class = "cardish",
                h3("Lifespan Analysis (Multifactor)"),
                
                uiOutput("groupvar_selector"),
                uiOutput("level_selectors"),
                
                br(),
                h4("Kaplan–Meier (overall)"),
                plotOutput("km_overall", height = "340px"),
                
                br(),
                h4("Kaplan–Meier (by selected factor)"),
                uiOutput("which_factor_plot_ui"),
                
                uiOutput("colour_selectors"),
                plotOutput("km_by_factor", height = "380px"),
                
                br(),
                h4("Cox model: pairwise comparisons (selected factor)"),
                tableOutput("cox_pairwise"),
                
                br(),
                h4("All selected factors (KM + Cox stats)"),
                uiOutput("all_factor_panels"),
                
                br(),
                h4("Interaction (combined plot + stats)"),
                uiOutput("interaction_selector"),
                plotOutput("km_interaction", height = "420px"),
                br(),
                h5("Interaction model output"),
                verbatimTextOutput("cox_interaction_text"),
                br(),
                h5("Pairwise comparisons across interaction groups"),
                tableOutput("cox_interaction_pairwise"),
                
                br(),
                actionButton("back_to_curation", "← Back to curation", class = "btn-secondary")
            )
        )
      )
    }
  })
}

shinyApp(ui, server)
