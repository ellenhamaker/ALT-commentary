# PREAMBLE ----
library(shiny)

# SETTINGS ----

# Order of the growth components in all 4 x 4 matrices
lcm_labels <- c("&Delta;<sub>X</sub>", "&Delta;<sub>Y</sub>",
                "&Gamma;<sub>X</sub>", "&Gamma;<sub>Y</sub>")
alt_labels <- c("A<sub>X</sub>", "A<sub>Y</sub>",
                "B<sub>X</sub>", "B<sub>Y</sub>")

# Default LCM-SR covariance matrix (numerical example in the supplement)
Sigma_default <- matrix(
  c(1.00, 0.50, 0.25, 0.20,
    0.50, 1.00, 0.10, 0.05,
    0.25, 0.10, 1.00, 0.30,
    0.20, 0.05, 0.30, 1.00),
  nrow = 4, byrow = TRUE
)

# HELPER FUNCTIONS ----

# Wrap a set of 16 cells (row-wise) in brackets, with row and column labels
labeled_matrix <- function(label_html, cells, row_labels, col_labels) {
  div(
    class = "equation-row",
    span(class = "equation-label", HTML(label_html)),
    span(class = "equation-equals", "="),
    div(
      class = "labeled-matrix",
      div(),
      div(class = "col-labels",
          lapply(col_labels, function(l) div(class = "axis-label", HTML(l)))),
      div(class = "row-labels",
          lapply(row_labels, function(l) div(class = "axis-label", HTML(l)))),
      div(
        class = "matrix-container",
        div(class = "matrix-bracket left-bracket"),
        div(class = "matrix-grid-4", cells),
        div(class = "matrix-bracket right-bracket")
      )
    )
  )
}

# Read-only display of a 4 x 4 matrix
matrix4_display <- function(label_html, mat, labels, highlight_negative = FALSE) {
  vals <- round(mat, 3)
  cells <- lapply(seq_len(16), function(k) {
    v <- t(vals)[k] # row-wise
    style <- if (highlight_negative && v < 0) "color: #d32f2f; font-weight: 600;" else NULL
    div(class = "matrix-value", style = style, v)
  })
  labeled_matrix(label_html, cells, labels, labels)
}

# Input grid for a symmetric 4 x 4 matrix: the lower triangle (incl. the
# diagonal) is editable; the upper triangle mirrors it
matrix4_input <- function(label_html, default, labels) {
  cells <- list()
  for (i in 1:4) {
    for (j in 1:4) {
      cells[[length(cells) + 1]] <- if (i >= j) {
        numericInput(
          paste0("S", i, j), label = NULL, value = default[i, j], step = 0.1,
          min = if (i == j) 0.001 else NA
        )
      } else {
        uiOutput(paste0("S", i, j, "_mirror"))
      }
    }
  }
  labeled_matrix(label_html, cells, labels, labels)
}

# 2 x 2 input for R
R_input <- function() {
  div(
    class = "equation-row",
    span(class = "equation-label", "R"),
    span(class = "equation-equals", "="),
    div(
      class = "matrix-container",
      div(class = "matrix-bracket left-bracket"),
      div(
        class = "matrix-input",
        numericInput("R11", label = NULL, value = 0.5, step = 0.1),
        numericInput("R12", label = NULL, value = 0.4, step = 0.1),
        numericInput("R21", label = NULL, value = 0,   step = 0.1),
        numericInput("R22", label = NULL, value = 0.5, step = 0.1)
      ),
      div(class = "matrix-bracket right-bracket")
    )
  )
}

# Transformation matrix T such that (A', B')' = T (Delta', Gamma')', with
# A = (I - R) Delta + R Gamma and B = (I - R) Gamma
compute_T <- function(R) {
  I <- diag(2)
  rbind(
    cbind(I - R, R),
    cbind(matrix(0, 2, 2), I - R)
  )
}

# WEB PAGE (UI) ----

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      width = 5,

      h4("Population Model (LCM-SR)"),
      p(
        "Set the covariance matrix of the intercepts and slopes of X and Y. ",
        "Edit the diagonal and the lower triangle; the upper triangle is ",
        "filled in automatically."
      ),

      matrix4_input("&Sigma;", Sigma_default, lcm_labels),

      uiOutput("pd_check"),

      br(),

      p("Set the values of the lagged parameters:"),

      R_input(),

      uiOutput("stationarity_check")
    ),

    # MAIN PANEL ----
    mainPanel(
      width = 7,

      h4("Implied Covariance Matrix (ALT Model)"),
      p(
        "The LCM-SR population values imply the following covariance matrix ",
        "of the intercepts and slopes in the ALT model. Negative values are ",
        "shown in red."
      ),

      uiOutput("alt_cov"),

      tags$hr(class = "thin-hr"),

      h4("Implied Correlations"),
      p("For easier comparison, the same matrices expressed as correlations:"),

      div(
        class = "equation-pair",
        uiOutput("lcm_cor"),
        uiOutput("alt_cor")
      )
    )
  ),

  # CSS ----
  tags$head(
    tags$style(
      HTML("

      :root {
        --cell-w: 54px;
        --cell-h: 32px;
        --bracket-w: 8px;
      }

      /* =================================================
         EQUATION ROW (label = matrix)
         ================================================= */

      .equation-row {
        display: flex;
        align-items: center;
        gap: 6px;
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

      .equation-pair {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-start;
        gap: 24px;
      }

      /* =================================================
         LABELED 4 x 4 MATRIX
         ================================================= */

      .labeled-matrix {
        display: grid;
        grid-template-columns: auto auto;
        column-gap: 4px;
        row-gap: 2px;
      }

      .col-labels {
        display: grid;
        grid-template-columns: repeat(4, var(--cell-w));
        column-gap: 4px;
        margin-left: calc(var(--bracket-w) + 4px);
      }

      .row-labels {
        display: flex;
        flex-direction: column;
        row-gap: 3px;
      }

      .axis-label {
        height: var(--cell-h);
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 13px;
        color: #555;
      }

      .col-labels .axis-label {
        height: auto;
      }

      .matrix-grid-4 {
        display: grid;
        grid-template-columns: repeat(4, var(--cell-w));
        gap: 3px 4px;
        margin-left: 4px;
        margin-right: 4px;
      }

      .matrix-grid-4 .shiny-input-container,
      .matrix-input .shiny-input-container {
        width: auto;
      }

      /* =================================================
         VECTOR AND MATRIX LAYOUT
         ================================================= */

      .matrix-container {
        display: flex;
        align-items: stretch;
        justify-content: center;
      }

      .matrix-input {
        display: grid;
        grid-template-columns: var(--cell-w) var(--cell-w);
        gap: 3px 4px;
        margin-left: 4px;
        margin-right: 4px;
      }

      .matrix-grid-4 .form-group,
      .matrix-input .form-group {
        margin-bottom: 0;
      }

      /* =================================================
         INPUT BOXES
         ================================================= */

      .matrix-grid-4 input,
      .matrix-input input {
        box-sizing: border-box;
        width: var(--cell-w);
        height: var(--cell-h);
        padding: 2px;
        text-align: center;
        font-size: 14px;
      }

      .matrix-input-mirrored {
        background-color: #f5f5f5 !important;
        color: #666;
        cursor: not-allowed;
      }

      /* =================================================
         BRACKETS
         ================================================= */

      .matrix-bracket {
        position: relative;
        width: var(--bracket-w);
        align-self: stretch;
      }

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
         PRINTED MATRIX VALUES (read-only bracket display)
         ================================================= */

      .matrix-value {
        box-sizing: border-box;
        width: var(--cell-w);
        height: var(--cell-h);
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 14px;
      }

      /* =================================================
         MESSAGES AND DIVIDERS
         ================================================= */

      .check-ok   { color: #2e7d32; }
      .check-warn { color: #d32f2f; font-weight: 600; }

      .thin-hr {
        border: none;
        border-top: 1px solid #ccc;
        margin: 18px 0;
      }
    ")
    )
  )
)

# SERVER ----
server <- function(input, output, session) {

  # LAGGED PARAMETERS ----
  R_mat <- reactive({
    matrix(
      c(input$R11, input$R12, input$R21, input$R22),
      nrow = 2, ncol = 2, byrow = TRUE
    )
  })

  # LCM-SR COVARIANCE MATRIX (symmetric, built from the lower triangle) ----
  Sigma_LCM <- reactive({
    S <- matrix(NA_real_, 4, 4)
    for (i in 1:4) {
      for (j in 1:i) {
        v <- input[[paste0("S", i, j)]]
        S[i, j] <- S[j, i] <- if (is.null(v) || is.na(v)) NA_real_ else v
      }
    }
    validate(need(!anyNA(S), "Please fill in all elements of \u03a3."))
    S
  })

  # Mirrored (read-only) upper-triangle cells
  for (i in 1:3) {
    for (j in (i + 1):4) {
      local({
        ii <- i
        jj <- j
        output[[paste0("S", ii, jj, "_mirror")]] <- renderUI({
          tags$input(
            type = "number",
            class = "form-control matrix-input-mirrored",
            value = input[[paste0("S", jj, ii)]],
            disabled = NA
          )
        })
      })
    }
  }

  # Is the LCM-SR covariance matrix positive definite?
  is_pd <- reactive({
    min(eigen(Sigma_LCM(), symmetric = TRUE, only.values = TRUE)$values) > 0
  })

  output$pd_check <- renderUI({
    if (is_pd()) {
      p(class = "check-ok", "\u03a3 is positive definite.")
    } else {
      p(class = "check-warn",
        "\u03a3 is not positive definite, so it is not a valid covariance matrix.")
    }
  })

  output$stationarity_check <- renderUI({
    req(!anyNA(R_mat()))
    eig <- eigen(R_mat(), only.values = TRUE)$values
    txt <- paste0("Eigenvalues of R: ", paste(format(round(eig, 3)), collapse = ", "))
    if (max(Mod(eig)) < 1) {
      p(class = "check-ok", txt, " (stationary).")
    } else {
      p(class = "check-warn", txt,
        " (NOT stationary: the ALT model and LCM-SR are then not equivalent).")
    }
  })

  # IMPLIED ALT COVARIANCE MATRIX ----
  # Sigma_ALT = T Sigma_LCM T', which contains Sigma_AA, Sigma_AB (upper right),
  # Sigma_BA (lower left) and Sigma_BB as its blocks
  Sigma_ALT <- reactive({
    req(!anyNA(R_mat()))
    Tm <- compute_T(R_mat())
    Tm %*% Sigma_LCM() %*% t(Tm)
  })

  output$alt_cov <- renderUI({
    matrix4_display("&Sigma;<sub>ALT</sub>", Sigma_ALT(), alt_labels,
                    highlight_negative = TRUE)
  })

  # CORRELATIONS ----
  output$lcm_cor <- renderUI({
    validate(need(is_pd(), "Correlations require a positive definite \u03a3."))
    matrix4_display("P<sub>LCM-SR</sub>", cov2cor(Sigma_LCM()), lcm_labels,
                    highlight_negative = TRUE)
  })

  output$alt_cor <- renderUI({
    validate(need(is_pd(), ""))
    matrix4_display("P<sub>ALT</sub>", cov2cor(Sigma_ALT()), alt_labels,
                    highlight_negative = TRUE)
  })
}

# RUN APP ----
shinyApp(
  ui = ui,
  server = server
)