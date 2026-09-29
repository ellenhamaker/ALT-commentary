# PREAMBLE ----
library(shiny)
library(mvtnorm)
library(ggplot2)

# HELPER FUNCTIONS ----

vector_bracket_display <- function(label_html, vec) {

  div(
    class = "equation-row",
    span(class = "equation-label", HTML(label_html)),
    span(class = "equation-equals", "="),
    div(
      class = "matrix-container",
      div(class = "matrix-bracket left-bracket"),
      div(
        class = "vector-input",
        div(class = "matrix-value", round(vec[1], 3)),
        div(class = "matrix-value", round(vec[2], 3))
      ),
      div(class = "matrix-bracket right-bracket")
    )
  )
}

# WEB PAGE (UI) ----

ui <- fluidPage(

  sidebarLayout(

    sidebarPanel(

      # POPULATION MODEL ----

      h4("Population Model"),

      p(
        "Choose which model you want to use as the population model:\n"
      ),
      radioButtons(
        "model_choice",
        label = NULL,
        choices = c(
          "LCM-SR" = "lcmToAlt",
          "ALT"    = "altToLcm"
        ),
        selected = "lcmToAlt"
      ),

      # GROWTH FACTORS ----

      p(
        "Set the values of the growth factors:"
      ),

      conditionalPanel(
        condition = "input.model_choice == 'lcmToAlt'",

        div(
          class = "equation-pair",
          div(
            class = "equation-row",
            span(class = "equation-label", HTML("&Delta;")),
            span(class = "equation-equals", "="),
            div(
              class = "matrix-container",
              div(class = "matrix-bracket left-bracket"),
              div(
                class = "vector-input",
                numericInput("D1", label = NULL, value = 5),
                numericInput("D2", label = NULL, value = 15)
              ),
              div(class = "matrix-bracket right-bracket")
            )
          ),

          div(
            class = "equation-row",
            span(class = "equation-label", HTML("&Gamma;")),
            span(class = "equation-equals", "="),
            div(
              class = "matrix-container",
              div(class = "matrix-bracket left-bracket"),
              div(
                class = "vector-input",
                numericInput("G1", label = NULL, value = 1),
                numericInput("G2", label = NULL, value = 3)
              ),
              div(class = "matrix-bracket right-bracket")
            )
          )
        )
      ),

      conditionalPanel(
        condition = "input.model_choice == 'altToLcm'",

        div(
          class = "equation-pair",

          div(
            class = "equation-row",
            span(class = "equation-label", "A"),
            span(class = "equation-equals", "="),
            div(
              class = "matrix-container",
              div(class = "matrix-bracket left-bracket"),
              div(
                class = "vector-input",
                numericInput("A1", label = NULL, value = 5),
                numericInput("A2", label = NULL, value = 15)
              ),
              div(class = "matrix-bracket right-bracket")
            )
          ),

          div(
            class = "equation-row",
            span(class = "equation-label", "B"),
            span(class = "equation-equals", "="),
            div(
              class = "matrix-container",
              div(class = "matrix-bracket left-bracket"),
              div(
                class = "vector-input",
                numericInput("B1", label = NULL, value = 1),
                numericInput("B2", label = NULL, value = 3)
              ),
              div(class = "matrix-bracket right-bracket")
            )
          )
        )
      ),

      br(),

      # LAGGED PARAMETERS ----
      p(
        "Set the values of lagged parameters:"
      ),

      div(
        class = "equation-row",
        span(class = "equation-label", "R"),
        span(class = "equation-equals", "="),
        div(
          class = "matrix-container",
          div(class = "matrix-bracket left-bracket"),
          div(
            class = "matrix-input",
            numericInput("R11", label = NULL, value = 0.5),
            numericInput("R12", label = NULL, value = 0.4),
            numericInput("R21", label = NULL, value = 0),
            numericInput("R22", label = NULL, value = 0.5)
          ),
          div(class = "matrix-bracket right-bracket")
        )
      ),

      tags$hr(),

      # MODEL CHARACTERISTICS ----

      h4("Model Characteristics"),

      p(
        "The above population values imply:"
      ),

      uiOutput("eigenvalues"),

      tags$hr(),

      # ADVANCED SETTINGS ----
      tags$details(

        tags$summary(
          span(class = "chevron", HTML("&#9656;")),
          h4("Advanced Settings")
        ),

        br(),

        numericInput(
          "nt",
          "Number of occasions",
          value = 10,
          min = 2,
          max = 50,
          step = 1
        ),

        h5("Residual covariance matrix E:"),

        div(
          class = "equation-row",
          span(class = "equation-label", HTML("&Sigma;<sub>E</sub>")),
          span(class = "equation-equals", "="),
          div(
            class = "matrix-container",
            div(class = "matrix-bracket left-bracket"),
            div(
              class = "matrix-input",
              numericInput("E11", label = NULL, value = 5, min = 0.001),
              numericInput("E12", label = NULL, value = 2, step = 0.5),
              numericInput("E21", label = NULL, value = 2, step = 0.5),
              numericInput("E22", label = NULL, value = 5, min = 0.001)
            ),
            div(class = "matrix-bracket right-bracket")
          )
        ),

        br(),

        numericInput(
          "seed",
          "Random seed",
          value = 123,
          min = 1,
          step = 1
        ),

        br(),

        actionButton(
          "sample_btn",
          "Sample"
        )
      )
    ),

    # MAIN PANEL ----

    mainPanel(

      # PARAMETER COMPARISON ----
      h4("Parameter Comparison"),

      uiOutput("comparison_text"),

      uiOutput("comparison_matrices"),

      tags$hr(),

      # VISUALIZATION ----

      withMathJax(
        p(
          "The trends in \\(X\\) and \\(Y\\) as implied by the population model and the set parameter values are visualized below."
        )
      ),

      plotOutput(
        "timeseries",
        height = "500px"
      )
    )
  ),

  # CSS ----
  tags$head(
    tags$style(
      HTML("

      /* =================================================
         EQUATION ROW (label = matrix)
         ================================================= */

      .equation-row {
        display: flex;
        align-items: center;
        gap: 8px;
        margin-bottom: 10px;
      }

      .equation-row .matrix-container {
        margin: 0;
        justify-content: flex-start;
      }

      .equation-label {
        font-size: 18px;
        font-weight: 600;
        min-width: 1.4em;
        text-align: right;
      }

      .equation-equals {
        font-size: 18px;
      }


      /* =================================================
         EQUATION PAIR (two matrices grouped side by side)
         ================================================= */

      .equation-pair {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-start;
        gap: 24px;
        margin-bottom: 4px;
      }

      .equation-pair .equation-row {
        margin-bottom: 6px;
      }


      /* =================================================
         VECTOR AND MATRIX LAYOUT
         ================================================= */

      .matrix-container {
        display: flex;
        align-items: stretch;
        justify-content: center;
        margin-top: 8px;
        margin-bottom: 12px;
      }


      /* =================================================
         VECTOR INPUT
         ================================================= */

      .vector-input {
        display: flex;
        flex-direction: column;
        gap: 3px;
        margin-left: 6px;
        margin-right: 6px;
      }


      /* =================================================
         MATRIX INPUT
         ================================================= */

      .matrix-input {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 3px 6px;
        margin-left: 6px;
        margin-right: 6px;
      }


      /* =================================================
         REMOVE SHINY LABEL SPACE
         ================================================= */

      .vector-input .form-group,
      .matrix-input .form-group {
        margin-bottom: 0;
      }


      /* =================================================
         INPUT BOXES
         ================================================= */

      .vector-input input,
      .matrix-input input {
        width: 65px;
        height: 32px;
        padding: 4px;
        text-align: center;
        font-size: 15px;
      }


      /* =================================================
         BRACKETS
         ================================================= */

      .matrix-bracket {
        position: relative;
        width: 10px;
        align-self: stretch;
      }


      /* =================================================
         LEFT BRACKET
         ================================================= */

      .left-bracket {
        border-left: 1px solid #444;
      }

      .left-bracket::before,
      .left-bracket::after {
        content: '';
        position: absolute;
        left: 0;
        width: 6px;
        height: 0;
        border-top: 1px solid #444;
      }

      .left-bracket::before { top: 0; }
      .left-bracket::after  { bottom: 0; }


      /* =================================================
         RIGHT BRACKET
         ================================================= */

      .right-bracket {
        border-right: 1px solid #444;
      }

      .right-bracket::before,
      .right-bracket::after {
        content: '';
        position: absolute;
        right: 0;
        width: 6px;
        height: 0;
        border-top: 1px solid #444;
      }

      .right-bracket::before { top: 0; }
      .right-bracket::after  { bottom: 0; }


      /* =================================================
         OUTPUT VECTORS / MATRICES
         ================================================= */

      .output-vector {
        display: flex;
        justify-content: center;
      }

      pre {
        border: none;
        background-color: transparent;
        font-size: 15px;
      }


      /* =================================================
         PRINTED MATRIX VALUES (read-only bracket display)
         ================================================= */

      .matrix-value {
        width: 65px;
        height: 32px;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 15px;
      }


      /* =================================================
         COLLAPSIBLE DETAILS / SUMMARY (Advanced Settings)
         ================================================= */

      details > summary {
        display: flex;
        align-items: center;
        gap: 6px;
        cursor: pointer;
        list-style: none;
        margin-bottom: 8px;
      }

      details > summary::-webkit-details-marker {
        display: none;
      }

      details > summary h4 {
        margin: 0;
      }

      details > summary .chevron {
        display: inline-block;
        font-size: 14px;
        transition: transform 0.15s ease-in-out;
      }

      details[open] > summary .chevron {
        transform: rotate(90deg);
      }

      details[open] > summary {
        margin-bottom: 12px;
      }

    ")
    )
  )
)

# SERVER ----

server <- function(input, output, session) {

  # R MATRIX ----
  R <- reactive({
    matrix(
      c(input$R11, input$R12, input$R21, input$R22),
      nrow = 2,
      ncol = 2,
      byrow = TRUE
    )
  })

  # COVARIANCE MATRIX OF E ----
  SigmaE <- reactive({
    matrix(
      c(input$E11, input$E12, input$E21, input$E22),
      nrow = 2,
      ncol = 2,
      byrow = TRUE
    )
  })


  # INPUT VECTORS ----
  D_in <- reactive(matrix(c(input$D1, input$D2), nrow = 2, ncol = 1))
  G_in <- reactive(matrix(c(input$G1, input$G2), nrow = 2, ncol = 1))
  A_in <- reactive(matrix(c(input$A1, input$A2), nrow = 2, ncol = 1))
  B_in <- reactive(matrix(c(input$B1, input$B2), nrow = 2, ncol = 1))

  # DELTA (D) ----
  # Delta = (I - R)^{-1} A - (I - R)^{-2} R B
  D <- reactive({
    if (input$model_choice == "altToLcm") {
      IR_inv <- solve(diag(2) - R())
      IR_inv %*% A_in() - IR_inv %*% IR_inv %*% R() %*% B_in()
    } else {
      D_in()
    }
  })

  # GAMMA (G) ----
  # Gamma = (I - R)^{-1} B
  G <- reactive({
    if (input$model_choice == "altToLcm") {
      solve(diag(2) - R()) %*% B_in()
    } else {
      G_in()
    }
  })

  # A ----
  # A = (I - R) Delta + R Gamma
  A <- reactive({
    if (input$model_choice == "lcmToAlt") {
      (diag(2) - R()) %*% D_in() + R() %*% G_in()
    } else {
      A_in()
    }
  })

  # B ----
  # B = (I - R) Gamma
  B <- reactive({
    if (input$model_choice == "lcmToAlt") {
      (diag(2) - R()) %*% G_in()
    } else {
      B_in()
    }
  })

  # DISPLAY D / DELTA ----
  output$D_display <- renderPrint({
    round(D(), 3)
  })

  # DISPLAY G / GAMMA ----
  output$G_display <- renderPrint({
    round(G(), 3)
  })

  # DISPLAY A ----
  output$A_display <- renderPrint({
    round(A(), 3)
  })

  # DISPLAY B ----
  output$B_display <- renderPrint({
    round(B(), 3)
  })

  # PARAMETER COMPARISON: TEXT ----

  # Select label for user-suplied and computed model
  output$comparison_text <- renderUI({
    selected_label   <- if (input$model_choice == "lcmToAlt") "LCM-SR" else "ALT"
    unselected_label <- if (input$model_choice == "lcmToAlt") "ALT"    else "LCM-SR"

    p(
      paste0(
        "The ", selected_label,
        " population values imply the following values of the growth components in the ",
        unselected_label, " model: "
      )
    )
  })

  # PARAMETER COMPARISON: MATRICES ----
  output$comparison_matrices <- renderUI({
    if (input$model_choice == "lcmToAlt") {
      div(
        class = "equation-pair",
        vector_bracket_display("A", A()),
        vector_bracket_display("B", B())
      )
    } else {
      div(
        class = "equation-pair",
        vector_bracket_display("&Delta;", D()),
        vector_bracket_display("&Gamma;", G())
      )
    }
  })

  # STATIONARITY ----
  output$eigenvalues <- renderUI({
    eig <- eigen(R())$values
    max_eig <- max(Mod(eig))
    eig_text <- paste(
      format(round(eig, 3)),
      collapse = ", "
    )

    stationarity_text <- if (max_eig < 1) {
      "R is stationary"
    } else {
      "R is NOT stationary"
    }

    tags$ul(
      tags$li(paste0("Eigenvalues: ", eig_text)),
      tags$li(stationarity_text)
    )
  })

  # SAMPLE BUTTON ----

  # Count presses on "Sample" button (works through random seed)
  sample_count <- reactiveVal(0)

  observeEvent(input$sample_btn, {
    sample_count(sample_count() + 1)
  })

  # SIMULATION ----
  simulated_data <- reactive({
    set.seed(input$seed + sample_count())

    # Number of occasions
    nt <- input$nt

    # Sample residual errors
    E <- rmvnorm(nt, c(0, 0), SigmaE())

    # Initialize innovation and observed vectors
    Z <- matrix(0, nrow = nt, ncol = 2)
    Y <- matrix(0, nrow = nt, ncol = 2)

    # Initial occasion
    Z[1, ] <- E[1, ]
    Y[1, ] <- D() + Z[1, ]

    # Loop through occasions
    for (t in 2:nt) {
      Z[t, ] <- R() %*% Z[t - 1, ] + E[t, ]
      Y[t, ] <- D() + G() * (t - 1) + Z[t, ]
    }

    colnames(Y) <- c("X", "Y")

    Y
  })

  # PLOT ----
  output$timeseries <- renderPlot({

    # Simulate data
    Y <- simulated_data()

    # Number of occasions
    nt <- input$nt

    # x-axis values
    occasions <- 0:(nt - 1)

    # Model trend values
    trend_X <- D()[1] + G()[1] * occasions
    trend_Y <- D()[2] + G()[2] * occasions

    plot_data <- data.frame(
      occasion = rep(occasions, 2),
      value    = c(Y[, "X"], Y[, "Y"]),
      variable = factor(
        rep(c("X", "Y"), each = nt),
        levels = c("X", "Y")
      )
    )

    # Create plot data for model trend
    trend_data <- data.frame(
      occasion = rep(occasions, 2),
      value    = c(trend_X, trend_Y),
      variable = factor(
        rep(c("X", "Y"), each = nt),
        levels = c("X", "Y")
      )
    )

    # Colors for X and Y (shared by data points and trend lines)
    series_colors <- c(
      "X" = "blue",
      "Y" = "orange"
    )

    # Line types for the trend lines only (not shown in the legend)
    series_linetypes <- c(
      "X" = "solid",
      "Y" = "dashed"
    )

    ggplot() +
      geom_point(
        data  = plot_data,
        aes(x = occasion, y = value, color = variable),
        alpha = 0.4,
        size  = 2
      ) +
      geom_line(
        data  = trend_data,
        aes(x = occasion, y = value, color = variable, linetype = variable),
        linewidth = 1
      ) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
      scale_color_manual(values = series_colors, name = NULL) +
      scale_linetype_manual(values = series_linetypes, guide = "none") +
      labs(x = "Occasion", y = "Value") +
      theme_minimal(base_size = 18) +
      theme(
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = NA)
      )

  })
}

# RUN APP ----
shinyApp(
  ui = ui,
  server = server
)