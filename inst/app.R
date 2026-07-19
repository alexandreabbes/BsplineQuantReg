# BsplineQuantReg Shiny Interface
# Author: Alexandre Abbes
# Inspired by the tkinter interface in Quant_reg_tk.py
#
# Run with:
# shiny::runApp(system.file("shiny", package = "BsplineQuantReg"))

library(shiny)
library(shinythemes)
library(shinyjs)
#library(BsplineQuantReg)
library(ggplot2)
library(DT)
library(plotly)

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),

  # Title
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

      # Contraintes
      h4("3. Contraintes", class = "text-primary"),

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
                            column(6, colourInput("curve_color", NULL, value = "blue")),
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
                 verbatimTextOutput( placeholder = TRUE, outputId = "r_code_box"),
                 br(),
                 actionButton("copy_code", "📋 Copier", icon = icon("copy"), class = "btn-sm btn-info"),
                 actionButton("export_code", "💾 Exporter", icon = icon("save"), class = "btn-sm btn-success")
        ),

        # Onglet Contraintes par région
        tabPanel("🎯 Contraintes par région", value = "regions",
                 br(),
                 fluidRow(
                   column(6,
                          h4("Régions définies"),
                          verbatimTextOutput("regions_info", placeholder = TRUE),
                          br(),
                          actionButton("clear_regions", "Effacer toutes les régions",
                                       class = "btn-sm btn-danger")
                   ),
                   column(6,
                          h4("Ajouter une région"),
                          p("Cliquez sur le graphique pour définir la zone"),
                          br(),
                          fluidRow(
                            column(4, numericInput("region_xmin", "X min:", value = 0.3, step = 0.05)),
                            column(4, numericInput("region_xmax", "X max:", value = 0.6, step = 0.05)),
                            column(4, numericInput("region_monot", "Monotonie:",
                                                   value = 0, min = -1, max = 1, step = 1))
                          ),
                          fluidRow(
                            column(6, numericInput("region_conv", "Convexité:",
                                                   value = 0, min = -1, max = 1, step = 1)),
                            column(6, numericInput("region_der3", "Dérivée 3e:",
                                                   value = 0, min = -1, max = 1, step = 1))
                          ),
                          actionButton("add_region", "Ajouter région",
                                       class = "btn-sm btn-primary")
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
    region_lines = list(),
    curve_lines = list(),
    data_name = "Aucune donnée"
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
      showNotification("Données test générées", type = "message")
    })
  })

  # Données température5
  observeEvent(input$temp_data, {
    withProgress(message = "Chargement des données température...", {
      # Températures 1880-1992
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
      # Normaliser les années
      x <- (years - 1880) / (1992 - 1880)
      y <- temp_data
      values$xtab <- x
      values$ytab <- y
      values$data_name <- "Température globale (1880-1992)"
      values$fit <- NULL
      values$curve_lines <- list()
      # Nœuds spécifiques
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
      # Remplacer les fonctions Python par R
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
      showNotification("Données générées", type = "message")
    }, error = function(e) {
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })

  # ============ NOEUDS ============

  # Définir les nœuds automatiquement
  observe({
    if (!is.null(values$xtab)) {
      kn <- input$knots_count
      values$knots <- quantile(values$xtab, probs = seq(0, 1, length.out = kn + 1))
    }
  })

  # ============ REGRESSION ============

  observeEvent(input$run, {
    req(values$xtab, values$ytab, values$knots)

    withProgress(message = "Régression en cours...", {

      degree <- input$degree
      tau <- input$tau
      solver <- input$solver
      monot <- as.numeric(input$monot)
      conv <- as.numeric(input$conv)
      der3 <- as.numeric(input$der3)

      if (degree < 3) der3 <- 0

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

        # Ajouter la courbe à la liste
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
      legend = list(orientation = "h", y = -0.1)
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

  # ============ CODE R ============

  output$r_code <- renderText({
    if (is.null(values$fit)) {
      return("# Lancez d'abord une régression")
    }

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
      "                       monot = ", as.numeric(input$monot), ",\n",
      "                       convcons = ", as.numeric(input$conv), ",\n",
      "                       der3cons = ", ifelse(input$degree >= 3, as.numeric(input$der3), 0), ",\n",
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
    showNotification("Courbes effacées", type = "info")
  })

  observeEvent(input$clear_all, {
    values$xtab <- NULL
    values$ytab <- NULL
    values$knots <- NULL
    values$fit <- NULL
    values$curve_lines <- list()
    values$data_name <- "Aucune donnée"
    showNotification("Tout effacé", type = "info")
  })

  # ============ DEMOS ============

  observeEvent(input$demo_comp, {
    showNotification("Lancement de la démo Comprehensive...", type = "info")
    demo("comprehensive", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_log, {
    showNotification("Lancement de la démo Logistic...", type = "info")
    demo("logistic", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_temp, {
    showNotification("Lancement de la démo Température...", type = "info")
    demo("temperature", package = "BsplineQuantReg")
  })

  observeEvent(input$demo_conv, {
    showNotification("Lancement de la démo Convexité...", type = "info")
    demo("convexity", package = "BsplineQuantReg")
  })

  # ============ EXPORT ============

  observeEvent(input$export_png, {
    showNotification("Export PNG à implémenter avec ggplot2", type = "warning")
  })

  observeEvent(input$copy_code, {
    # Utiliser clipr ou JS pour copier
    showNotification("Copie du code...", type = "info")
  })

  observeEvent(input$export_code, {
    showNotification("Export du code...", type = "info")
  })

  # ============ REGIONS ============

  output$regions_info <- renderPrint({
    cat("Fonctionnalité à implémenter avec les contraintes par région")
  })

  observeEvent(input$add_region, {
    showNotification("Région ajoutée (à implémenter)", type = "info")
  })

  observeEvent(input$clear_regions, {
    showNotification("Régions effacées", type = "info")
  })

  # ============ IMPORT CSV/EXCEL ============

  observeEvent(input$load_csv, {
    showNotification("Import CSV à implémenter avec shinyFile", type = "warning")
  })

  observeEvent(input$load_excel, {
    showNotification("Import Excel à implémenter avec shinyFile", type = "warning")
  })
}

# Run app
shinyApp(ui = ui, server = server)
