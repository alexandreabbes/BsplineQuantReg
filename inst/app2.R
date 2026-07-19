# BsplineQuantReg Shiny Interface
# Author: Alexandre Abbes
# Inspired by the tkinter interface in Quant_reg_tk.py
#
# Run with:
# shiny::runApp("inst/app.R")

library(shiny)
library(shinythemes)
library(shinyjs)
library(BsplineQuantReg)
library(ggplot2)
library(DT)
library(plotly)
library(colourpicker)

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),

  # Custom CSS pour les notifications
  tags$style(HTML("
    .shiny-notification-success {
      background-color: #d4edda;
      border-color: #c3e6cb;
      color: #155724;
    }
    .region-box {
      border: 2px solid #ff9800;
      background-color: rgba(255, 152, 0, 0.15);
      border-radius: 5px;
      padding: 8px;
      margin: 4px 0;
    }
  ")),

  titlePanel(
    h1("BsplineQuantReg - Régression Quantile avec Splines Contraintes",
       align = "center", style = "color: #2c3e50;"),
    windowTitle = "BsplineQuantReg - Spline Quantile Regression"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "background-color: #f8f9fa; border-radius: 5px;",

      # Onglet principal : Données
      h4("1. Données", class = "text-primary"),

      fluidRow(
        column(6, actionButton("load_csv", "📂 CSV", icon = icon("file-csv"),
                               class = "btn-sm btn-info", style = "width:100%;")),
        column(6, actionButton("load_excel", "📊 Excel", icon = icon("file-excel"),
                               class = "btn-sm btn-info", style = "width:100%;"))
      ),
      br(),
      fluidRow(
        column(6, actionButton("test_data", "🧪 Test", icon = icon("flask"),
                               class = "btn-sm btn-success", style = "width:100%;")),
        column(6, actionButton("temp_data", "🌡️ Temp", icon = icon("thermometer-half"),
                               class = "btn-sm btn-warning", style = "width:100%;"))
      ),
      br(),
      p("Fonction personnalisée:"),
      fluidRow(
        column(8, textInput("custom_func", NULL,
                           value = "2*x + 0.5*sin(6*pi*x) + 0.05*rnorm(n)",
                           placeholder = "f(x)")),
        column(4, numericInput("n_points", "n", value = 100, min = 10, max = 1000, step = 10))
      ),
      actionButton("generate_custom", "Générer", icon = icon("play"),
                   class = "btn-sm btn-primary", style = "width:100%;"),

      hr(),

      # Paramètres spline
      h4("2. Spline", class = "text-primary"),

      fluidRow(
        column(6, numericInput("degree", "Degré:", value = 3, min = 1, max = 4, step = 1)),
        column(6, numericInput("knots_count", "Nœuds:", value = 10, min = 4, max = 30, step = 1))
      ),

      fluidRow(
        column(6, sliderInput("tau", "τ (quantile):", min = 0.05, max = 0.95, value = 0.5, step = 0.05)),
        column(6, selectInput("solver", "Solveur:",
                             choices = c("CLARABEL", "OSQP", "ECOS", "SCS")))
      ),

      hr(),

      # Contraintes - Mode
      h4("3. Contraintes", class = "text-primary"),

      radioButtons("constraint_mode", "Mode:",
                   choices = c("Uniformes" = "uniform", "Par région" = "region"),
                   selected = "uniform", inline = TRUE),

      # Contraintes uniformes
      conditionalPanel(
        condition = "input.constraint_mode == 'uniform'",
        radioButtons("monot", "Monotonie:",
                     choices = c("✗" = "0", "↗" = "1", "↘" = "-1"),
                     selected = "0", inline = TRUE),

        radioButtons("conv", "Convexité:",
                     choices = c("✗" = "0", "∪" = "1", "∩" = "-1"),
                     selected = "0", inline = TRUE),

        conditionalPanel(
          condition = "input.degree >= 3",
          radioButtons("der3", "Dérivée 3e:",
                       choices = c("✗" = "0", "+" = "1", "-" = "-1"),
                       selected = "0", inline = TRUE)
        )
      ),

      # Contraintes par région
      conditionalPanel(
        condition = "input.constraint_mode == 'region'",
        div(style = "font-size: 13px; color: #555; margin-bottom: 10px;",
            "Cliquez sur le graphique pour définir une région, puis ajustez les paramètres ci-dessous."),

        fluidRow(
          column(6, numericInput("region_xmin", "X min:", value = 0.3, step = 0.05)),
          column(6, numericInput("region_xmax", "X max:", value = 0.6, step = 0.05))
        ),

        radioButtons("region_monot", "Monotonie:",
                     choices = c("✗" = "0", "↗" = "1", "↘" = "-1"),
                     selected = "0", inline = TRUE),

        radioButtons("region_conv", "Convexité:",
                     choices = c("✗" = "0", "∪" = "1", "∩" = "-1"),
                     selected = "0", inline = TRUE),

        conditionalPanel(
          condition = "input.degree >= 3",
          radioButtons("region_der3", "Dérivée 3e:",
                       choices = c("✗" = "0", "+" = "1", "-" = "-1"),
                       selected = "0", inline = TRUE)
        ),

        fluidRow(
          column(6, actionButton("add_region", "➕ Ajouter région",
                                 class = "btn-sm btn-primary", style = "width:100%;")),
          column(6, actionButton("clear_regions", "🗑️ Effacer tout",
                                 class = "btn-sm btn-danger", style = "width:100%;"))
        ),
        br(),
        div(id = "regions_list", style = "max-height: 120px; overflow-y: auto;")
      ),

      hr(),

      # Boutons d'exécution
      actionButton("run", "▶ Lancer", icon = icon("play"),
                   class = "btn-success btn-lg", style = "width:100%; margin-bottom: 5px;"),

      fluidRow(
        column(6, actionButton("clear_all", "Effacer tout", icon = icon("trash"),
                               class = "btn-sm btn-danger", style = "width:100%;")),
        column(6, actionButton("clear_curves", "Effacer courbes", icon = icon("eraser"),
                               class = "btn-sm btn-warning", style = "width:100%;"))
      ),

      hr(),

      # Démos
      h4("4. Démos", class = "text-primary"),
      div(style = "display: flex; flex-wrap: wrap; gap: 5px;",
          actionButton("demo_comp", "Comprehensive", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_log", "Logistic", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_temp", "Température", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_conv", "Convexité", class = "btn-sm btn-info", style = "flex:1;")
      ),

      hr(),

      # Info package
      div(style = "font-size: 11px; color: #666; text-align: center;",
          p("BsplineQuantReg v0.1.0"),
          p("Basé sur Karlin-Studden (1966)"),
          a("GitHub", href = "https://github.com/alexandreabbes/BsplineQuantReg", target = "_blank")
      )
    ),

    mainPanel(
      width = 9,

      tabsetPanel(
        id = "tabs",

        # Onglet Plot
        tabPanel("📈 Visualisation", value = "plot",
                 br(),
                 fluidRow(
                   column(10, plotlyOutput("spline_plot", height = "500px")),
                   column(2,
                          h5("Couleur:"),
                          fluidRow(
                            column(6, colourpicker::colourInput("curve_color", NULL, value = "blue")),
                            column(6, actionButton("apply_color", "Appliquer", class = "btn-sm"))
                          ),
                          br(),
                          p("Courbes:", textOutput("curve_count", inline = TRUE)),
                          br(),
                          actionButton("export_png", "📊 Exporter PNG", class = "btn-sm btn-success", width = "100%")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(6,
                          h5("Information"),
                          verbatimTextOutput("fit_info", placeholder = TRUE)
                   ),
                   column(6,
                          h5("Coefficients"),
                          verbatimTextOutput("coef_info", placeholder = TRUE)
                   )
                 )
        ),

        # Onglet Données
        tabPanel("📊 Données", value = "data",
                 br(),
                 fluidRow(
                   column(6,
                          h4("Résumé des données"),
                          verbatimTextOutput("data_summary")
                   ),
                   column(6,
                          h4("Nœuds"),
                          verbatimTextOutput("knots_info")
                   )
                 ),
                 br(),
                 h4("Tableau des données"),
                 DTOutput("data_table")
        ),

        # Onglet Code
        tabPanel("📝 Code R", value = "code",
                 br(),
                 h4("Code pour reproduire l'analyse:"),
                 verbatimTextOutput("r_code"),
                 br(),
                 actionButton("copy_code", "📋 Copier", icon = icon("copy"), class = "btn-sm btn-info"),
                 actionButton("export_code", "💾 Exporter", icon = icon("save"), class = "btn-sm btn-success")
        ),

        # Onglet Régions
        tabPanel("🎯 Régions", value = "regions",
                 br(),
                 h4("Régions définies"),
                 verbatimTextOutput("regions_info", placeholder = TRUE),
                 br(),
                 fluidRow(
                   column(6,
                          h5("Régions actives"),
                          uiOutput("regions_list_ui")
                   ),
                   column(6,
                          h5("Instructions"),
                          p("1. Passez en mode 'Par région' dans le panneau de gauche"),
                          p("2. Cliquez sur le graphique pour sélectionner une zone"),
                          p("3. Ajustez les paramètres dans le panneau"),
                          p("4. Cliquez sur 'Ajouter région' pour valider")
                   )
                 )
        )
      )
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # Reactive values
  values <- reactiveValues(
    xtab = NULL,
    ytab = NULL,
    knots = NULL,
    fit = NULL,
    x_eval = NULL,
    y_eval = NULL,
    curve_lines = list(),
    regions = list(),
    data_name = "Aucune donnée",
    region_id = 0,
    click_data = NULL
  )

  # ============ GENERATION DES DONNEES ============

  # Données test
  observeEvent(input$test_data, {
    withProgress(message = "Génération des données...", {
      set.seed(42)
      n <- 200
      x <- seq(0, 1, length.out = n)
      y <- 2*x + 0.2*sin(10*pi*x) + 0.05*rnorm(n)
      values$xtab <- x
      values$ytab <- y
      values$data_name <- "Test (sinusoïde)"
      values$fit <- NULL
      values$curve_lines <- list()
      values$regions <- list()
      showNotification("Données test générées", type = "message")
    })
  })

  # Données température
  observeEvent(input$temp_data, {
    withProgress(message = "Chargement des données température...", {
      temp_data <- c(
        -0.32, -0.32, -0.40, -0.39, -0.65, -0.43, -0.40, -0.52, -0.30, -0.12,
        -0.40, -0.42, -0.39, -0.45, -0.35, -0.36, -0.19, -0.14, -0.37, -0.22,
         0.00, -0.08, -0.24, -0.36, -0.49, -0.27, -0.19, -0.43, -0.29, -0.30,
        -0.29, -0.29, -0.28, -0.23, -0.04, -0.02, -0.24, -0.42, -0.35, -0.16,
        -0.17, -0.09, -0.13, -0.16, -0.14, -0.14,  0.10, -0.03,  0.03, -0.18,
        -0.06,  0.04,  0.02, -0.13,  0.03, -0.06,  0.02,  0.13,  0.13, -0.03,
         0.15,  0.12,  0.10,  0.04,  0.11, -0.04,  0.01,  0.13, -0.01, -0.06,
        -0.14, -0.02,  0.04,  0.14, -0.07, -0.06, -0.17,  0.10,  0.10,  0.05,
        -0.01,  0.08,  0.02,  0.02, -0.26, -0.16, -0.09, -0.02, -0.12,  0.03,
         0.04, -0.11, -0.07,  0.19, -0.07, -0.05, -0.22,  0.16,  0.09,  0.14,
         0.28,  0.39,  0.07,  0.29,  0.11,  0.11,  0.16,  0.32,  0.35,  0.25,
         0.47,  0.41,  0.13
      )
      years <- 1880:1992
      x <- (years - 1880) / (1992 - 1880)
      y <- temp_data
      values$xtab <- x
      values$ytab <- y
      values$data_name <- "Température globale (1880-1992)"
      values$fit <- NULL
      values$curve_lines <- list()
      values$regions <- list()
      year_knots <- c(1880, 1889, 1900, 1910, 1930, 1940, 1965, 1992)
      knots <- (year_knots - 1880) / (1992 - 1880)
      values$knots <- knots
      showNotification("Données température chargées", type = "message")
    })
  })

  # Données personnalisées
  observeEvent(input$generate_custom, {
    tryCatch({
      n <- input$n_points
      x <- seq(0, 1, length.out = n)
      func_str <- input$custom_func
      func_str <- gsub("sin\\(", "sin(", func_str)
      func_str <- gsub("cos\\(", "cos(", func_str)
      func_str <- gsub("pi", "pi", func_str)
      func_str <- gsub("randn\\(", "rnorm(", func_str)

      y <- eval(parse(text = func_str))
      values$xtab <- x
      values$ytab <- y
      values$data_name <- "Fonction personnalisée"
      values$fit <- NULL
      values$curve_lines <- list()
      values$regions <- list()
      showNotification("Données générées", type = "message")
    }, error = function(e) {
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })

  # ============ NOEUDS ============

  observe({
    if (!is.null(values$xtab)) {
      kn <- input$knots_count
      values$knots <- quantile(values$xtab, probs = seq(0, 1, length.out = kn + 1))
    }
  })

  # ============ REGIONS ============

  # Ajouter une région
  observeEvent(input$add_region, {
    req(values$xtab, values$knots)

    xmin <- input$region_xmin
    xmax <- input$region_xmax

    if (xmin >= xmax) {
      showNotification("X min doit être inférieur à X max", type = "warning")
      return()
    }

    values$region_id <- values$region_id + 1
    region <- list(
      id = values$region_id,
      xmin = xmin,
      xmax = xmax,
      monot = as.numeric(input$region_monot),
      conv = as.numeric(input$region_conv),
      der3 = as.numeric(input$region_der3)
    )

    values$regions <- c(values$regions, list(region))

    showNotification(paste("Région ajoutée:", xmin, "-", xmax), type = "message")
  })

  # Effacer toutes les régions
  observeEvent(input$clear_regions, {
    values$regions <- list()
    values$region_id <- 0
    showNotification("Régions effacées", type = "message")
  })

  # Supprimer une région spécifique
  observeEvent(input$delete_region, {
    id <- as.numeric(input$delete_region)
    values$regions <- values$regions[!sapply(values$regions, function(r) r$id == id)]
    showNotification(paste("Région", id, "supprimée"), type = "message")
  })

  # ============ CONSTRUCTION DES CONTRAINTES ============

  build_constraints <- function() {
    degree <- input$degree
    kn <- length(values$knots) - 1

    if (input$constraint_mode == "uniform") {
      monot <- rep(as.numeric(input$monot), kn)
      conv <- rep(as.numeric(input$conv), kn + 1)
      # Pour degré 3 : der3 est un vecteur de longueur kn (par intervalle)
      der3 <- rep(as.numeric(input$der3), kn)
    } else {
      monot <- rep(0, kn)
      conv <- rep(0, kn + 1)
      der3 <- rep(0, kn)  # ← CORRECTION : kn et non kn+1

      for (region in values$regions) {
        for (i in 1:kn) {
          x1 <- values$knots[i]
          x2 <- values$knots[i + 1]
          if (x2 > region$xmin && x1 < region$xmax) {
            if (region$monot != 0) monot[i] <- region$monot
            if (region$conv != 0) {
              conv[i] <- region$conv
              conv[i + 1] <- region$conv
            }
            # Pour degré 3 : der3 par intervalle
            if (region$der3 != 0 && degree >= 3) {
              der3[i] <- region$der3  # ← CORRECTION : pas der3[i+1]
            }
          }
        }
      }
    }

    list(monot = monot, conv = conv, der3 = der3)
  }
  # ============ REGRESSION ============

  observeEvent(input$run, {
    req(values$xtab, values$ytab, values$knots)

    withProgress(message = "Régression en cours...", {

      degree <- input$degree
      tau <- input$tau
      solver <- input$solver

      constraints <- build_constraints()
      monot <- constraints$monot
      conv <- constraints$conv
      der3 <- constraints$der3

      if (degree < 3) der3 <- rep(0, length(der3))

      incProgress(0.3, detail = "Construction des B-splines...")

      fit <- tryCatch({
        quantile_spline(
          values$xtab, values$ytab, values$knots, tau,
          degree = degree,
          monot = monot,
          convcons = conv,
          der3cons = der3,
          solver = solver,
          verbose = FALSE,
          callable = TRUE
        )
      }, error = function(e) {
        showNotification(paste("Erreur:", e$message), type = "error")
        NULL
      })

      if (!is.null(fit)) {
        incProgress(0.8, detail = "Évaluation...")
        x_eval <- seq(min(values$xtab), max(values$xtab), length.out = 300)
        y_eval <- fit(x_eval)

        values$fit <- fit
        values$x_eval <- x_eval
        values$y_eval <- y_eval

        color <- input$curve_color
        values$curve_lines <- c(values$curve_lines, list(list(x = x_eval, y = y_eval, color = color)))

        showNotification("Régression réussie!", type = "message")
      }

      incProgress(1)
    })
  })

  # ============ VISUALISATION ============

  output$spline_plot <- renderPlotly({
    if (is.null(values$xtab)) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = "Chargez des données"))
    }

    p <- plot_ly()

    # Données
    p <- p %>% add_trace(
      x = values$xtab, y = values$ytab,
      type = "scatter", mode = "markers",
      marker = list(color = "gray", size = 6, opacity = 0.5),
      name = "Données"
    )

    # Régions (en mode région)
    if (input$constraint_mode == "region" && length(values$regions) > 0) {
      for (region in values$regions) {
        y_range <- range(values$ytab)
        p <- p %>% add_trace(
          x = c(region$xmin, region$xmax, region$xmax, region$xmin, region$xmin),
          y = c(y_range[1], y_range[1], y_range[2], y_range[2], y_range[1]),
          type = "scatter", mode = "lines",
          fill = "toself",
          fillcolor = "rgba(255, 152, 0, 0.15)",
          line = list(color = "rgba(255, 152, 0, 0.8)", width = 1),
          name = paste0("Région [", region$xmin, ", ", region$xmax, "]"),
          hoverinfo = "text",
          text = paste0(
            "Région ", region$id, "\n",
            "Monotonie: ", c("✗", "↗", "↘")[region$monot + 2], "\n",
            "Convexité: ", c("✗", "∪", "∩")[region$conv + 2], "\n",
            "Dérivée 3e: ", c("✗", "+", "-")[region$der3 + 2]
          )
        )
      }
    }

    # Courbes sauvegardées
    for (curve in values$curve_lines) {
      p <- p %>% add_trace(
        x = curve$x, y = curve$y,
        type = "scatter", mode = "lines",
        line = list(color = curve$color, width = 2),
        name = paste0("τ=", input$tau)
      )
    }

    # Nœuds
    if (!is.null(values$knots)) {
      y_range <- range(values$ytab)
      y_pos <- y_range[2] - 0.1 * diff(y_range)
      p <- p %>% add_trace(
        x = values$knots, y = rep(y_pos, length(values$knots)),
        type = "scatter", mode = "markers",
        marker = list(color = "red", symbol = "triangle-down", size = 10),
        name = "Nœuds"
      )
    }

    p <- p %>% layout(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      hovermode = "closest",
      legend = list(orientation = "h", y = -0.1),
      dragmode = "select"
    )

    p
  })

  # ============ INFORMATION ============

  output$fit_info <- renderPrint({
    if (is.null(values$fit)) {
      cat("Aucune régression effectuée")
    } else {
      cat("Statut: Réussi\n")
      cat("Degré:", attr(values$fit, "degree"), "\n")
      cat("τ:", input$tau, "\n")
      cat("Nœuds:", length(values$knots), "\n")
      cat("Coefficients:", length(attr(values$fit, "coefficients")), "\n")
    }
  })

  output$coef_info <- renderPrint({
    if (is.null(values$fit)) {
      cat("Aucun coefficient")
    } else {
      coefs <- attr(values$fit, "coefficients")
      cat("Min:", round(min(coefs), 4), "\n")
      cat("Max:", round(max(coefs), 4), "\n")
      cat("Moyenne:", round(mean(coefs), 4), "\n")
      cat("Ecart-type:", round(sd(coefs), 4), "\n")
      cat("Premiers coefficients:\n")
      print(head(coefs, 5))
    }
  })

  output$curve_count <- renderText({
    length(values$curve_lines)
  })

  # ============ REGIONS UI ============

  output$regions_list_ui <- renderUI({
    if (length(values$regions) == 0) {
      return(p("Aucune région définie", style = "color: #999;"))
    }

    tags$div(
      lapply(values$regions, function(r) {
        tags$div(
          class = "region-box",
          tags$div(
            style = "display: flex; justify-content: space-between; align-items: center;",
            tags$span(
              style = "font-weight: bold;",
              paste0("Région ", r$id, " [", round(r$xmin, 2), ", ", round(r$xmax, 2), "]")
            ),
            actionButton(paste0("del_", r$id), "✕",
                         class = "btn-sm btn-danger",
                         style = "padding: 0px 6px;",
                         onclick = paste0("Shiny.setInputValue('delete_region', ", r$id, ")"))
          ),
          tags$div(
            style = "font-size: 12px; color: #555; margin-top: 4px;",
            paste0(
              "M: ", c("✗", "↗", "↘")[r$monot + 2], " | ",
              "C: ", c("✗", "∪", "∩")[r$conv + 2], " | ",
              "D3: ", c("✗", "+", "-")[r$der3 + 2]
            )
          )
        )
      })
    )
  })

  # ============ DONNEES ============

  output$data_summary <- renderPrint({
    if (is.null(values$xtab)) {
      cat("Aucune donnée")
    } else {
      cat("Source:", values$data_name, "\n")
      cat("Nombre de points:", length(values$xtab), "\n")
      cat("x: [", min(values$xtab), ",", max(values$xtab), "]\n")
      cat("y: [", min(values$ytab), ",", max(values$ytab), "]")
    }
  })

  output$knots_info <- renderPrint({
    if (is.null(values$knots)) {
      cat("Aucun nœud")
    } else {
      cat("Nombre de nœuds:", length(values$knots), "\n")
      cat("Nœuds:\n")
      print(round(values$knots, 4))
    }
  })

  output$data_table <- renderDT({
    if (is.null(values$xtab)) return(NULL)
    df <- data.frame(
      x = round(values$xtab, 4),
      y = round(values$ytab, 4)
    )
    datatable(df, options = list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip'
    ))
  })

  output$regions_info <- renderPrint({
    if (length(values$regions) == 0) {
      cat("Aucune région définie\n")
      cat("Passez en mode 'Par région' et ajoutez une région.")
    } else {
      cat("Régions définies:\n")
      for (r in values$regions) {
        cat(paste0(
          "  ", r$id, ": [", round(r$xmin, 2), ", ", round(r$xmax, 2), "]  ",
          "M=", c("✗", "↗", "↘")[r$monot + 2], " ",
          "C=", c("✗", "∪", "∩")[r$conv + 2], " ",
          "D3=", c("✗", "+", "-")[r$der3 + 2], "\n"
        ))
      }
    }
  })

  # ============ CODE R ============

  output$r_code <- renderText({
    if (is.null(values$fit)) {
      return("# Lancez d'abord une régression")
    }

    constraints <- build_constraints()

    paste0(
      "library(BsplineQuantReg)\n\n",
      "# Données\n",
      "x <- c(", paste(round(head(values$xtab, 5), 4), collapse = ", "), ", ...)\n",
      "y <- c(", paste(round(head(values$ytab, 5), 4), collapse = ", "), ", ...)\n\n",
      "# Nœuds\n",
      "knots <- c(", paste(round(values$knots, 4), collapse = ", "), ")\n\n",
      "# Régression\n",
      "fit <- quantile_spline(x, y, knots,\n",
      "                       tau = ", input$tau, ",\n",
      "                       degree = ", input$degree, ",\n",
      "                       monot = c(", paste(constraints$monot, collapse = ", "), "),\n",
      "                       convcons = c(", paste(constraints$conv, collapse = ", "), "),\n",
      "                       der3cons = c(", paste(constraints$der3, collapse = ", "), "),\n",
      "                       solver = '", input$solver, "',\n",
      "                       callable = TRUE)\n\n",
      "# Évaluation\n",
      "x_eval <- seq(min(x), max(x), length.out = 300)\n",
      "y_eval <- fit(x_eval)\n\n",
      "# Visualisation\n",
      "plot(x, y, pch = 16, cex = 0.5, col = 'gray')\n",
      "lines(x_eval, y_eval, col = '", input$curve_color, "', lwd = 2)"
    )
  })

  # ============ ACTIONS ============

  observeEvent(input$apply_color, {
    showNotification(paste("Couleur appliquée:", input$curve_color), type = "message")
  })

  observeEvent(input$clear_curves, {
    values$curve_lines <- list()
    showNotification("Courbes effacées", type = "message")
  })

  observeEvent(input$clear_all, {
    values$xtab <- NULL
    values$ytab <- NULL
    values$knots <- NULL
    values$fit <- NULL
    values$curve_lines <- list()
    values$regions <- list()
    values$region_id <- 0
    values$data_name <- "Aucune donnée"
    showNotification("Tout effacé", type = "message")
  })

  # ============ DEMOS ============

  observeEvent(input$demo_comp, {
    showNotification("Lancement de la démo Comprehensive...", type = "message")
    demo("comprehensive", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_log, {
    showNotification("Lancement de la démo Logistic...", type = "message")
    demo("logistic", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_temp, {
    showNotification("Lancement de la démo Température...", type = "message")
    demo("temperature", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_conv, {
    showNotification("Lancement de la démo Convexité...", type = "message")
    demo("convexity", package = "BsplineQuantReg")
  })

  # ============ EXPORT ============

  observeEvent(input$export_png, {
    showNotification("Export PNG: utilisez le bouton 'Télécharger' de Plotly", type = "message")
  })

  observeEvent(input$copy_code, {
    showNotification("Code copié dans le presse-papier", type = "message")
  })

  observeEvent(input$export_code, {
    showNotification("Export du code...", type = "message")
  })
}

# Run app
shinyApp(ui = ui, server = server)
