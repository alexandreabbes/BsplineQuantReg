# BsplineQuantReg Shiny Interface
# Author: Alexandre Abbes
# Version stable avec sélection graphique des régions
#
# Run with:
# shiny::runApp("inst/app.R")

library(shiny)
library(shinythemes)
library(shinyjs)
library(BsplineQuantReg)
library(DT)
library(plotly)
library(colourpicker)

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),

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
    windowTitle = "BsplineQuantReg"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "background-color: #f8f9fa; border-radius: 5px;",

      # ============ DONNEES ============
      h4("1. Données", class = "text-primary"),

      fluidRow(
        column(6, actionButton("test_data", "🧪 Test",
                               class = "btn-sm btn-success", style = "width:100%;")),
        column(6, actionButton("temp_data", "🌡️ Temp",
                               class = "btn-sm btn-warning", style = "width:100%;"))
      ),
      br(),

      h4("2. Spline", class = "text-primary"),
      fluidRow(
        column(6, numericInput("degree", "Degré:", value = 3, min = 1, max = 4)),
        column(6, numericInput("knots_count", "Nœuds:", value = 10, min = 4, max = 30))
      ),
      fluidRow(
        column(6, actionButton("add_knot_mode", "➕ Ajouter nœud",
                               class = "btn-sm btn-primary", style = "width:100%;")),
        column(6, actionButton("clear_knots", "🗑️ Effacer nœuds",
                               class = "btn-sm btn-danger", style = "width:100%;"))
      ),
      br(),

      fluidRow(
        column(6, sliderInput("tau", "τ:", min = 0.05, max = 0.95, value = 0.5)),
        column(6, selectInput("solver", "Solveur:",
                              choices = c("CLARABEL", "OSQP", "ECOS", "SCS")))
      )

      ,



      h5("Intervalle:"),
      fluidRow(
        column(6, numericInput("data_xmin", "X min:", value = 0, step = 0.05)),
        column(6, numericInput("data_xmax", "X max:", value = 1, step = 0.05))
      ),

      p("Fonction personnalisée:"),
      fluidRow(
        column(8, textInput("custom_func", NULL,
                            value = "2*x + 0.5*sin(6*pi*x) + 0.05*rnorm(n)")),
        column(4, numericInput("n_points", "n", value = 100, min = 10, max = 1000))
      ),
      actionButton("generate_custom", "Générer",
                   class = "btn-sm btn-primary", style = "width:100%;"),

      hr(),

      # ============ SPLINE ============
      h4("2. Spline", class = "text-primary"),

      fluidRow(
        column(6, numericInput("degree", "Degré:", value = 3, min = 1, max = 4)),
        column(6, numericInput("knots_count", "Nœuds:", value = 10, min = 4, max = 30))
      ),

      fluidRow(
        column(6, sliderInput("tau", "τ:", min = 0.05, max = 0.95, value = 0.5)),
        column(6, selectInput("solver", "Solveur:",
                              choices = c("CLARABEL", "OSQP", "ECOS", "SCS")))
      ),

      hr(),

      # ============ CONTRAINTES ============
      h4("3. Contraintes", class = "text-primary"),

      radioButtons("constraint_mode", "Mode:",
                   choices = c("Uniformes" = "uniform", "Par région" = "region"),
                   selected = "uniform", inline = TRUE),

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

      conditionalPanel(
        condition = "input.constraint_mode == 'region'",
        div(style = "font-size: 13px; color: #555; margin-bottom: 10px;",
            "1. Cliquez sur 'Sélectionner'"),
        div(style = "font-size: 13px; color: #555; margin-bottom: 10px;",
            "2. Sélectionnez une zone sur le graphique"),
        div(style = "font-size: 13px; color: #555; margin-bottom: 10px;",
            "3. Les champs X min/max se mettent à jour"),

        fluidRow(
          column(6, actionButton("start_selection", "🎯 Sélectionner",
                                 class = "btn-sm btn-warning", style = "width:100%;")),
          column(6, actionButton("clear_regions", "🗑️ Effacer",
                                 class = "btn-sm btn-danger", style = "width:100%;"))
        ),
        br(),

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
          column(6, actionButton("add_region", "➕ Ajouter",
                                 class = "btn-sm btn-primary", style = "width:100%;")),
          column(6, actionButton("update_region", "🔄 Mettre à jour",
                                 class = "btn-sm btn-info", style = "width:100%;"))
        ),
        br(),
        div(id = "regions_list", style = "max-height: 120px; overflow-y: auto;")
      ),

      hr(),

      # ============ EXECUTION ============
      actionButton("run", "▶ Lancer",
                   class = "btn-success btn-lg", style = "width:100%;"),

      fluidRow(
        column(6, actionButton("clear_all", "Effacer tout",
                               class = "btn-sm btn-danger", style = "width:100%;")),
        column(6, actionButton("clear_curves", "Effacer courbes",
                               class = "btn-sm btn-warning", style = "width:100%;"))
      ),

      hr(),

      h4("4. Démos", class = "text-primary"),
      div(style = "display: flex; flex-wrap: wrap; gap: 5px;",
          actionButton("demo_comp", "Comprehensive", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_log", "Logistic", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_temp", "Température", class = "btn-sm btn-info", style = "flex:1;"),
          actionButton("demo_conv", "Convexité", class = "btn-sm btn-info", style = "flex:1;")
      ),

      hr(),

      div(style = "font-size: 11px; color: #666; text-align: center;",
          p("BsplineQuantReg v0.1.0"),
          p("Basé sur Karlin-Studden (1966)"),
          a("GitHub", href = "https://github.com/alexandreabbes/BsplineQuantReg", target = "_blank")
      )
    ),

    mainPanel(
      width = 9,

      tabsetPanel(
        tabPanel("📈 Visualisation",
                 br(),
                 fluidRow(
                   column(10, plotlyOutput("spline_plot", height = "500px")),
                   column(2,
                          h5("Couleur:"),
                          colourpicker::colourInput("curve_color", NULL, value = "blue"),
                          actionButton("apply_color", "Appliquer", class = "btn-sm", style = "width:100%;"),
                          br(), br(),
                          p("Courbes:", textOutput("curve_count", inline = TRUE)),
                          br()
                          #,
                          #h5("Statut sélection:"),
                          #verbatimTextOutput("selection_status", placeholder = TRUE)
                   )
                 ),
                 br(),
                 fluidRow(
                   column(6, h5("Information"), verbatimTextOutput("fit_info")),
                   column(6, h5("Coefficients"), verbatimTextOutput("coef_info"))
                 )
        ),

        tabPanel("📊 Données",
                 br(),
                 fluidRow(
                   column(6, h4("Résumé"), verbatimTextOutput("data_summary")),
                   column(6, h4("Nœuds"), verbatimTextOutput("knots_info"))
                 ),
                 br(),
                 DTOutput("data_table")
        ),

        tabPanel("📝 Code R",
                 br(),
                 h4("Code pour reproduire l'analyse:"),
                 verbatimTextOutput("r_code"),
                 #htmlOutput("r_code_display"),
                 ),


        tabPanel("🎯 Régions",
                 br(),
                 h4("Régions définies"),
                 verbatimTextOutput("regions_info"),
                 br(),
                 fluidRow(
                   column(6, h5("Régions actives"), uiOutput("regions_list_ui")),
                   column(6, h5("Instructions"),
                          p("1. Mode 'Par région'"),
                          p("2. 'Sélectionner' puis rectangle sur le graphique"),
                          p("3. Ajuster les contraintes"),
                          p("4. 'Ajouter région'")
                   )
                 )
        )
      )
    )
  )
)

# SERVER ------------------------------------------------------------------

server <- function(input, output, session) {
  output$debug_output <- renderPrint({
    cat("Debug - Dernières valeurs:\n")
    cat("region_xmin:", input$region_xmin, "\n")
    cat("region_xmax:", input$region_xmax, "\n")
    cat("mode sélection:", values$selecting_region, "\n")
    cat("nb régions:", length(values$regions), "\n")
  })
  #
  # ============ REACTIVE VALUES ============
  values <- reactiveValues(
    xtab = NULL,
    ytab = NULL,
    knots = NULL,
    manual_knots = list(),  # ← Liste des nœuds ajoutés manuellement
    adding_knot = FALSE,    # ← État du mode ajout
    fit = NULL,
    x_eval = NULL,
    y_eval = NULL,
    curve_lines = list(),
    regions = list(),
    data_name = "Aucune donnée",
    region_id = 0,
    selected_region_id = NULL,
    selecting_region = FALSE
  )
  # ============ FONCTION DE MISE A JOUR DES CHAMPS ============

  update_region_fields <- function(xmin, xmax) {
    if (is.null(xmin) || is.null(xmax) || is.na(xmin) || is.na(xmax)) {
      return()
    }
    if (xmin >= xmax) {
      showNotification("X min doit être inférieur à X max", type = "warning")
      return()
    }
    updateNumericInput(session, "region_xmin", value = round(xmin, 3))
    updateNumericInput(session, "region_xmax", value = round(xmax, 3))
  }
  # ============ DEBUG : Récupération des valeurs sélectionnées ============

  debug_selection <- function(msg = "", xmin = NULL, xmax = NULL) {
  # Affiche également dans l'interface si vous voulez
    showNotification(
      paste("Debug:", msg, "| xmin:", round(xmin, 4), "xmax:", round(xmax, 4)),
      type = "message",
      duration = 3
    )
  }
  # ============ GENERATION DES DONNEES ============

  observeEvent(input$test_data, {
    withProgress(message = "Génération...", {
      set.seed(42)
      n <- 200
      xmin <- input$data_xmin
      xmax <- input$data_xmax
      x <- as.vector(seq(xmin, xmax, length.out = n))
      y <- as.vector(2*x + 0.2*sin(10*pi*x) + 0.05*rnorm(n))
      values$xtab <- x
      values$ytab <- y
      values$data_name <- paste("Test [", xmin, ",", xmax, "]")
      values$fit <- NULL
      values$curve_lines <- list()
      values$regions <- list()
      showNotification("Données test générées", type = "message")
    })
  })

  observeEvent(input$temp_data, {
    withProgress(message = "Chargement...", {
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
      values$data_name <- "Température (1880-1992)"
      values$fit <- NULL
      values$curve_lines <- list()
      values$regions <- list()
      year_knots <- c(1880, 1889, 1900, 1910, 1930, 1940, 1965, 1992)
      knots <- (year_knots - 1880) / (1992 - 1880)
      values$knots <- knots
      showNotification("Données température chargées", type = "message")
    })
  })

  output$r_code_display <- renderUI({
    code <- code_r()
    HTML(paste0(
      "<pre id='r_code_text' style='white-space: pre-wrap; background: #f5f5f5; padding: 10px; border-radius: 4px; font-family: monospace;'>",
      code,
      "</pre>"
    ))
  })

  observeEvent(input$generate_custom, {
    tryCatch({
      n <- input$n_points
      xmin <- input$data_xmin
      xmax <- input$data_xmax
      x <- as.vector(seq(xmin, xmax, length.out = n))
      func_str <- gsub("sin\\(", "sin(", input$custom_func)
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
  # ============ NOEUDS ============

  # Initialisation automatique des nœuds
  observe({
    if (!is.null(values$xtab) && length(values$manual_knots) == 0) {
      kn <- max(input$knots_count,1)
      values$knots <- quantile(values$xtab, probs = seq(0, 1, length.out = kn + 1))
    }
  })
  # Observer pour mettre à jour l'affichage quand les nœuds changent
  observe({
    # Forcer la mise à jour des outputs quand values$knots change
    output$knot_count_display <- renderPrint({
      if (!is.null(values$knots)) {
        cat("Nœuds:", length(values$knots))
      } else {
        cat("Aucun nœud")
      }
    })

    output$knots_list_display <- renderPrint({
      if (!is.null(values$knots) && length(values$knots) > 0) {
        knots_formatted <- round(values$knots, 4)
        cat(paste(knots_formatted, collapse = ", "))
      } else {
        cat("Aucun nœud défini")
      }
    })
  })
  # Activer/désactiver le mode ajout de nœuds
  observeEvent(input$add_knot_mode, {
    values$adding_knot <- !values$adding_knot
    if (values$adding_knot) {
      showNotification("Mode ajout de nœuds: cliquez sur le graphique", type = "message")
      updateActionButton(session, "add_knot_mode", label = "⏹ Arrêter nœuds")
    } else {
      updateActionButton(session, "add_knot_mode", label = "➕ Ajouter nœud")
    }
  })

  # Ajouter un nœud par clic
  observeEvent(event_data("plotly_click", source = "plot"), {
    if (values$adding_knot) {
      click <- event_data("plotly_click", source = "plot")
      if (!is.null(click) && !is.null(values$xtab)) {
        x <- click$x

        # Vérifier que le nœud est dans l'intervalle
        if (x > min(values$xtab) && x < max(values$xtab)) {
          # Éviter les doublons
          if (!any(abs(values$knots - x) < 1e-6)) {
            values$manual_knots <- c(values$manual_knots, x)
            values$knots <- sort(c(values$knots, x))
            showNotification(paste("Nœud ajouté à x =", round(x, 3)), type = "message")
          } else {
            showNotification("Ce nœud existe déjà", type = "warning")
          }
        } else {
          showNotification("Le nœud doit être à l'intérieur de l'intervalle", type = "warning")
        }
      }
    }
  })

  # Effacer les nœuds manuels
  observeEvent(input$clear_knots, {
    values$manual_knots <- list()
    if (!is.null(values$xtab)) {
      kn <- input$knots_count
      values$knots <- quantile(values$xtab, probs = seq(0, 1, length.out = kn + 1))
    }
    showNotification("Nœuds réinitialisés", type = "message")
  })
  observe({
    if (!is.null(values$xtab)) {
      kn <- input$knots_count
      values$knots <- quantile(values$xtab, probs = seq(0, 1, length.out = kn + 1))
    }
  })

  # ============ GESTION DE LA SELECTION ============

  observeEvent(input$start_selection, {
    values$selecting_region <- !values$selecting_region
    if (values$selecting_region) {
      showNotification("Sélectionnez une région sur le graphique (rectangle)", type = "message")
      updateActionButton(session, "start_selection", label = "⏹ Arrêter")
    } else {
      updateActionButton(session, "start_selection", label = " Sélectionner")
    }
  })

  # === FONCTION PRINCIPALE : Sélection par rectangle ===
  observeEvent(event_data("plotly_selected", source = "plot"), {
    # DEBUG : afficher l'événement brut
    selected <- event_data("plotly_selected", source = "plot")
    cat("\n=== EVENT DATA ===\n")
    print(str(selected))
    cat("==================\n")

    if (input$constraint_mode == "region" && values$selecting_region) {
      if (!is.null(selected) && nrow(selected) > 0) {
        x_vals <- selected$x
        cat("x_vals:", x_vals, "\n")
        cat("length:", length(x_vals), "\n")

        if (length(x_vals) >= 2) {
          xmin <- min(x_vals, na.rm = TRUE)
          xmax <- max(x_vals, na.rm = TRUE)

          # DEBUG : afficher les valeurs calculées
          debug_selection("Sélection rectangle", xmin, xmax)

          # Mettre à jour les champs
          updateNumericInput(session, "region_xmin", value = round(xmin, 3))
          updateNumericInput(session, "region_xmax", value = round(xmax, 3))

          values$selecting_region <- FALSE
          updateActionButton(session, "start_selection", label = "🎯 Sélectionner")
          showNotification(
            paste("Région sélectionnée: [", round(xmin, 3), ", ", round(xmax, 3), "]"),
            type = "message"
          )
        } else {
          debug_selection("Pas assez de points", NULL, NULL)
        }
      } else {
        debug_selection("selected est NULL ou vide", NULL, NULL)
      }
    }
  })

  # ============ REGIONS ============

  observeEvent(input$add_region, {
    req(values$xtab, values$knots)
    xmin <- input$region_xmin
    xmax <- input$region_xmax
    if (xmin >= xmax) {
      showNotification("X min < X max", type = "warning")
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
    values$selected_region_id <- NULL
    showNotification(paste("Région ajoutée: [", round(xmin, 3), ", ", round(xmax, 3), "]"), type = "message")
  })

  observeEvent(input$update_region, {
    if (!is.null(values$selected_region_id)) {
      idx <- which(sapply(values$regions, function(r) r$id == values$selected_region_id))
      if (length(idx) > 0) {
        values$regions[[idx]]$xmin <- input$region_xmin
        values$regions[[idx]]$xmax <- input$region_xmax
        values$regions[[idx]]$monot <- as.numeric(input$region_monot)
        values$regions[[idx]]$conv <- as.numeric(input$region_conv)
        values$regions[[idx]]$der3 <- as.numeric(input$region_der3)
        showNotification("Région mise à jour", type = "message")
      }
    } else {
      showNotification("Sélectionnez d'abord une région", type = "warning")
    }
  })

  observeEvent(input$clear_regions, {
    values$regions <- list()
    values$region_id <- 0
    values$selected_region_id <- NULL
    showNotification("Régions effacées", type = "message")
  })

  # Avec ignoreNULL = TRUE (par défaut, il ignore les NULL)
  observeEvent(input$delete_region, {
    id <- as.numeric(input$delete_region)

    # Vérifier que id est valide
    if (is.na(id)) {
      showNotification("ID invalide", type = "warning")
      return()
    }

    # Filtrer
    values$regions <- values$regions[!sapply(values$regions, function(r) r$id == id)]

    if (!is.null(values$selected_region_id) && values$selected_region_id == id) {
      values$selected_region_id <- NULL
    }

    showNotification(paste("Région", id, "supprimée"), type = "message")
  }, ignoreNULL = TRUE)  # ← Clé : ignorer les NULL

  # ============ CONSTRUCTION DES CONTRAINTES ============

  build_constraints <- function() {
    degree <- input$degree
    kn <- length(values$knots) - 1

    if (kn < 1) {
      showNotification("Pas assez de nœuds!", type = "warning")
      return(NULL)
    }

    safe_repeat <- function(val, len) {
      if (length(val) == 1) {
        return(rep(as.numeric(val), len))
      }
      v <- as.numeric(val)
      if (length(v) > len) return(v[1:len])
      if (length(v) < len) return(c(v, rep(0, len - length(v))))
      return(v)
    }

    if (input$constraint_mode == "uniform") {
      # Contraintes uniformes
      monot <- safe_repeat(input$monot, kn)
      conv <- safe_repeat(input$conv, kn + 1)
      der3 <- safe_repeat(input$der3, kn)
    } else {
      # Mode région - tout à 0 par défaut (aucune contrainte)
      monot <- rep(0, kn)
      conv <- rep(0, kn + 1)
      der3 <- rep(0, kn)

      # Appliquer les régions
      for (region in values$regions) {
        for (i in 1:kn) {
          x1 <- values$knots[i]
          x2 <- values$knots[i + 1]
          # Vérifier si l'intervalle est dans la région
          if (x2 > region$xmin && x1 < region$xmax) {
            # Appliquer les contraintes de la région
            if (region$monot != 0) {
              monot[i] <- region$monot
            }
            if (region$conv != 0) {
              conv[i] <- region$conv
              conv[i + 1] <- region$conv
            }
            if (region$der3 != 0 && degree >= 3) {
              der3[i] <- region$der3
            }
          }
        }
      }
    }

    if (degree < 3) der3 <- rep(0, kn)

    list(monot = monot, conv = conv, der3 = der3)
  }

  # ============ REGRESSION ============

  observeEvent(input$run, {

    req(values$xtab, values$ytab, values$knots)
    # Vérifier qu'il y a au moins 2 nœuds
    if (length(values$knots) < 2) {
      showNotification("Il faut au moins 2 nœuds!", type = "error")
      return()
    }

    # Vérifier que kn > 0
    if (length(values$knots) - 1 < 1) {
      showNotification("Il faut au moins 1 intervalle!", type = "error")
      return()
    }
    if (length(values$xtab) != length(values$ytab)) {
      showNotification("x et y longueurs différentes!", type = "error")
      return()
    }

    withProgress(message = "Régression...", {
      constraints <- build_constraints()
      if (is.null(constraints)) return()
      fit <- tryCatch({
        quantile_spline(
          as.vector(values$xtab),
          as.vector(values$ytab),
          as.vector(values$knots),
          tau = input$tau,
          degree = input$degree,
          monot = constraints$monot,
          convcons = constraints$conv,
          der3cons = constraints$der3,
          solver = input$solver,
          #verbose=verbose,
          callable = TRUE
        )
      }, error = function(e) {
        showNotification(paste("Erreur:", e$message), type = "error")
        NULL
      })
      if (!is.null(fit)) {
        x_eval <- seq(min(values$xtab), max(values$xtab), length.out = 300)
        y_eval <- fit(x_eval)
        values$fit <- fit
        values$x_eval <- x_eval
        values$y_eval <- y_eval
        color <- input$curve_color
        values$curve_lines <- c(values$curve_lines, list(list(x = x_eval, y = y_eval, color = color)))
        showNotification("Régression réussie!", type = "message")
      }
    })
  })


  # ============ VISUALISATION ============

  output$spline_plot <- renderPlotly({
    req(values$xtab)

    #  p <- plot_ly()
    p <- plot_ly(source = "plot")
    # Données
    p <- p %>% add_trace(
      x = values$xtab, y = values$ytab,
      type = "scatter", mode = "markers",
      marker = list(color = "gray", size = 6, opacity = 0.5),
      name = "Données"
    )
    # Dans output$spline_plot, ajoutez :

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

    # Régions
    if (input$constraint_mode == "region" && length(values$regions) > 0) {
      y_range <- range(values$ytab)
      for (region in values$regions) {
        is_selected <- !is.null(values$selected_region_id) && values$selected_region_id == region$id
        border_color <- if (is_selected) "#ff0000" else "rgba(255, 152, 0, 0.8)"
        fill_color <- if (is_selected) "rgba(255, 0, 0, 0.15)" else "rgba(255, 152, 0, 0.15)"
        p <- p %>% add_trace(
          x = c(region$xmin, region$xmax, region$xmax, region$xmin, region$xmin),
          y = c(y_range[1], y_range[1], y_range[2], y_range[2], y_range[1]),
          type = "scatter", mode = "lines",
          fill = "toself",
          fillcolor = fill_color,
          line = list(color = border_color, width = ifelse(is_selected, 3, 1)),
          name = paste0("Région ", region$id),
          hoverinfo = "text",
          text = paste0(
            "Région ", region$id, "\n",
            "[", round(region$xmin, 3), ", ", round(region$xmax, 3), "]\n",
            "M: ", c("✗", "↗", "↘")[region$monot + 2], "\n",
            "C: ", c("✗", "∪", "∩")[region$conv + 2], "\n",
            "D3: ", c("✗", "+", "-")[region$der3 + 2]
          )
        )
      }
    }

    # Courbes
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
    # Ajouter une annotation sur le plot
    p <- p %>% layout(
      annotations = list(
        x = 0.02,
        y = 0.98,
        text = paste("Nœuds:", length(values$knots)),
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
        font = list(size = 12, color = "red")
      )
    )
    p <- p %>% layout(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      hovermode = "closest",
      legend = list(orientation = "h", y = -0.1),
      dragmode = if (values$selecting_region && input$constraint_mode == "region") "select" else "zoom"
    )

    p <- p %>% config(
      scrollZoom = TRUE,
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("sendDataToCloud", "resetViews")
    )

    p
  })

  # ============ SELECTION D'UNE REGION PAR CLIC ============

  observeEvent(event_data("plotly_click", source = "plot"), {
    click <- event_data("plotly_click", source = "plot")

    # DEBUG : afficher le clic
    if (!is.null(click)) {
      cat("\n=== CLICK EVENT ===\n")
      cat("x:", click$x, "\n")
      cat("y:", click$y, "\n")
      cat("==================\n")
    }

    if (!values$selecting_region && input$constraint_mode == "region") {
      if (!is.null(click)) {
        x <- click$x
        cat("Recherche région pour x =", x, "\n")

        for (region in values$regions) {
          cat("  Vérification région", region$id, ": [", region$xmin, ", ", region$xmax, "]\n")
          if (x >= region$xmin && x <= region$xmax) {
            values$selected_region_id <- region$id
            debug_selection("Clic sur région", region$xmin, region$xmax)
            update_region_fields(region$xmin, region$xmax)
            updateRadioButtons(session, "region_monot", selected = as.character(region$monot))
            updateRadioButtons(session, "region_conv", selected = as.character(region$conv))
            updateRadioButtons(session, "region_der3", selected = as.character(region$der3))
            showNotification(paste("Région", region$id, "sélectionnée"), type = "message")
            break
          }
        }
      }
    }
  })

  # ============ OUTPUTS ============

  output$fit_info <- renderPrint({
    if (is.null(values$fit)) { cat("Aucune régression") } else {
      cat("Statut: Réussi\n")
      cat("Degré:", attr(values$fit, "degree"), "\n")
      cat("τ:", input$tau, "\n")
      cat("Nœuds:", length(values$knots), "\n")
      cat("Coefficients:", length(attr(values$fit, "coefficients")), "\n")
    }
  })

  output$coef_info <- renderPrint({
    if (is.null(values$fit)) { cat("Aucun coefficient") } else {
      coefs <- attr(values$fit, "coefficients")
      cat("Min:", round(min(coefs), 4), "\n")
      cat("Max:", round(max(coefs), 4), "\n")
      cat("Moyenne:", round(mean(coefs), 4), "\n")
    }
  })

  output$curve_count <- renderText({ length(values$curve_lines) })

  output$selection_status <- renderText({
    if (values$selecting_region) "Mode sélection actif" else
      if (!is.null(values$selected_region_id)) paste("Région sélectionnée:", values$selected_region_id) else
        "Aucune sélection"
  })

  output$regions_list_ui <- renderUI({
    if (length(values$regions) == 0) return(p("Aucune région", style = "color: #999;"))
    tags$div(lapply(values$regions, function(r) {
      is_selected <- !is.null(values$selected_region_id) && values$selected_region_id == r$id
      tags$div(
        class = "region-box",
        style = if (is_selected) "border: 3px solid #ff0000; background-color: rgba(255,0,0,0.1);",
        tags$div(style = "display: flex; justify-content: space-between;",
                 tags$span(style = "font-weight: bold;",
                           paste0("Région ", r$id, " [", round(r$xmin, 3), ", ", round(r$xmax, 3), "]")),
                 actionButton(paste0("del_", r$id), "✕", class = "btn-sm btn-danger",
                              style = "padding: 0px 6px;",
                              onclick = paste0("Shiny.setInputValue('delete_region', ", r$id, ")"))
        ),
        tags$div(style = "font-size: 12px; color: #555;",
                 paste0("M:", c("↘","✗","↗",)[r$monot+2], " | C:", c("∩","✗","∪")[r$conv+2],
                        " | D3:", c("-","✗","+")[r$der3+2])
        )
      )
    }))
  })

  output$data_summary <- renderPrint({
    if (is.null(values$xtab)) { cat("Aucune donnée") } else {
      cat("Source:", values$data_name, "\n")
      cat("Points:", length(values$xtab), "\n")
      cat("x: [", min(values$xtab), ",", max(values$xtab), "]\n")
      cat("y: [", min(values$ytab), ",", max(values$ytab), "]")
    }
  })

  output$knots_info <- renderPrint({
    if (is.null(values$knots)) { cat("Aucun nœud") } else {
      cat("Nœuds:", length(values$knots), "\n")
      print(round(values$knots, 4))
    }
  })

  output$data_table <- renderDT({
    if (is.null(values$xtab)) return(NULL)
    datatable(data.frame(x = round(values$xtab, 4), y = round(values$ytab, 4)),
              options = list(pageLength = 10, scrollX = TRUE))
  })

  output$regions_info <- renderPrint({
    if (length(values$regions) == 0) { cat("Aucune région") } else {
      for (r in values$regions) {
        cat(r$id, ": [", round(r$xmin, 3), ", ", round(r$xmax, 3), "]  M=",
            c("✗","↗","↘")[r$monot+2], " C=", c("✗","∪","∩")[r$conv+2],
            " D3=", c("✗","+","-")[r$der3+2], "\n")
      }
    }
  })

  # ============ CODE R ============
    observeEvent(input$copy_code, {
    runjs('
    var text = document.getElementById("r_code_text").innerText;
    navigator.clipboard.writeText(text).then(
      function() { alert("Code copié!"); },
      function(err) { alert("Erreur de copie: " + err); }
    );
  ')
  })

  output$r_code <- renderText({
    if (is.null(values$fit)) return("# Lancez d'abord une régression")
    constraints <- build_constraints()
    if (is.null(constraints)) return("# Erreur: contraintes non définies")
    paste0(
      "library(BsplineQuantReg)\n\n",
      "x <- c(", paste(round(values$xtab, 4), collapse = ", "), ")\n",
      "y <- c(", paste(round(values$ytab, 4), collapse = ", "), ")\n",
      "knots <- c(", paste(round(values$knots, 4), collapse = ", "), ")\n\n",
      "fit <- quantile_spline(x, y, knots,\n",
      "                       tau = ", input$tau, ",\n",
      "                       degree = ", input$degree, ",\n",
      "                       monot = c(", paste(constraints$monot, collapse = ", "), "),\n",
      "                       convcons = c(", paste(constraints$conv, collapse = ", "), "),\n",
      "                       der3cons = c(", paste(constraints$der3, collapse = ", "), "),\n",
      "                       solver = '", input$solver, "',\n",
      "                       callable = TRUE)\n\n",
      "x_eval <- seq(min(x), max(x), length.out = 300)\n",
      "y_eval <- fit(x_eval)\n\n",
      "plot(x, y, pch = 16, cex = 0.5, col = 'gray')\n",
      "lines(x_eval, y_eval, col = '", input$curve_color, "', lwd = 2)"
    )
  })


  # ============ ACTIONS ============

  observeEvent(input$apply_color, {
    showNotification(paste("Couleur:", input$curve_color), type = "message")
  })

  observeEvent(input$clear_curves, {
    values$curve_lines <- list()
    showNotification("Courbes effacées", type = "message")
  })

  observeEvent(input$clear_all, {
    values$xtab <- NULL; values$ytab <- NULL; values$knots <- NULL
    values$fit <- NULL; values$curve_lines <- list()
    values$regions <- list(); values$region_id <- 0
    values$selected_region_id <- NULL
    values$data_name <- "Aucune donnée"
    showNotification("Tout effacé", type = "message")
  })

  # ============ DEMOS ============

  lapply(c("demo_comp", "demo_log", "demo_temp", "demo_conv"), function(id) {
    observeEvent(input[[id]], {
      demo_name <- switch(id,
                          demo_comp = "comprehensive",
                          demo_log = "logistic",
                          demo_temp = "temperature",
                          demo_conv = "convexity")
      showNotification(paste("Démo:", demo_name), type = "message")
      demo(demo_name, package = "BsplineQuantReg")
    })
  })
}
1
# Run app
shinyApp(ui = ui, server = server)
