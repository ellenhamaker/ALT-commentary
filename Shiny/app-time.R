# PREAMBLE ----
library(shiny)
library(bslib)
library(mvtnorm)
library(ggplot2)

# HELPER FUNCTIONS ----

vector_bracket_display <- function(label_html, vec, highlight_negative = FALSE) {

  # Round the values to 3 decimal places
  val1 <- round(vec[1], 3)
  val2 <- round(vec[2], 3)

  # Highlight negative values in red
  style1 <- if (highlight_negative && val1 < 0) "color: #d32f2f; font-weight: 600;" else NULL
  style2 <- if (highlight_negative && val2 < 0) "color: #d32f2f; font-weight: 600;" else NULL

  div(
    class = "equation-row",
    span(class = "equation-label", HTML(label_html)),
    span(class = "equation-equals", "="),
    div(
      class = "matrix-container",
      div(class = "matrix-bracket left-bracket"),
      div(
        class = "vector-input",
        div(class = "matrix-value", style = style1, val1),
        div(class = "matrix-value", style = style2, val2)
      ),
      div(class = "matrix-bracket right-bracket")
    )
  )
}

# Function to display a matrix in a bracketed format
matrix_bracket_display <- function(label_html, mat, highlight_negative = FALSE) {

  vals <- round(mat, 3)

  cell_style <- function(v) {
    if (highlight_negative && v < 0) "color: #d32f2f; font-weight: 600;" else NULL
  }

  div(
    class = "equation-row",
    span(class = "equation-label", HTML(label_html)),
    span(class = "equation-equals", "="),
    div(
      class = "matrix-container",
      div(class = "matrix-bracket left-bracket"),
      div(
        class = "matrix-value-grid",
        div(class = "matrix-value", style = cell_style(vals[1, 1]), vals[1, 1]),
        div(class = "matrix-value", style = cell_style(vals[1, 2]), vals[1, 2]),
        div(class = "matrix-value", style = cell_style(vals[2, 1]), vals[2, 1]),
        div(class = "matrix-value", style = cell_style(vals[2, 2]), vals[2, 2])
      ),
      div(class = "matrix-bracket right-bracket")
    )
  )
}

# Function to compute ALT growth factors from LCM-SR
compute_AB <- function(R, D, G) {
  I <- diag(2)
  list(
    A = (I - R) %*% D + R %*% G,
    B = (I - R) %*% G
  )
}

# WEB PAGE (UI) ----

ui <- page_sidebar(

  # SIDEBAR ----
  sidebar = sidebar(
    width = 300,

    h4("Population Model"),
      p(
        "Set population model parameters for the LCM-SR for Scenario 1, which involves taking quarterly measurements."
      ),
    p(
      "Set the values of ", HTML("R"), ", ", HTML("&Delta;"), ", and ", HTML("&Gamma;"),
      " for Scenario 1, involving quarterly measurement occasions. The corresponding values ",
      "for Scenario 2 (biannual) and Scenario 3 (annual) are derived from ",
      "these and shown alongside them."
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
          numericInput("R11", label = NULL, value = 0.5, step = 0.1),
          numericInput("R12", label = NULL, value = 0.4, step = 0.1),
          numericInput("R21", label = NULL, value = 0, step = 0.1),
          numericInput("R22", label = NULL, value = 0.5, step = 0.1)
        ),
        div(class = "matrix-bracket right-bracket")
      )
    ),

    div(
      class = "equation-row",
      span(class = "equation-label", HTML("&Delta;")),
      span(class = "equation-equals", "="),
      div(
        class = "matrix-container",
        div(class = "matrix-bracket left-bracket"),
        div(
          class = "vector-input",
          numericInput("D1", label = NULL, value = 5, step = 0.1),
          numericInput("D2", label = NULL, value = 15, step = 0.1)
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
          numericInput("G1", label = NULL, value = 1, step = 0.1),
          numericInput("G2", label = NULL, value = 3, step = 0.1)
        ),
        div(class = "matrix-bracket right-bracket")
      )
    )
  ),

  # MAIN AREA ----
  
  layout_columns(
    col_widths = c(4, 4, 4),

    # SCENARIO 1 CARD ----
    card(
      card_header("Scenario 1 (Quarterly)"),
      card_body(
        fillable = FALSE,
        uiOutput("scenario1_matrices")
      )
    ),

    # SCENARIO 2 CARD ----
    card(
      card_header("Scenario 2 (Biannual)"),
      card_body(
        fillable = FALSE,
        uiOutput("scenario2_matrices")
      )
    ),

    # SCENARIO 3 CARD ----
    card(
      card_header("Scenario 3 (Annual)"),
      card_body(
        fillable = FALSE,
        uiOutput("scenario3_matrices")
      )
    )
  ),

  # CSS ----
  tags$head(
    tags$style(
      HTML("

      /* =================================================
         CARD BODY (tighten padding so the matrices have
         room to breathe inside each scenario card)
         ================================================= */

      .card-body {
        padding: 0.75rem !important;
        overflow-x: hidden;
      }

      /* =================================================
         EQUATION ROW (label = matrix)
         ================================================= */

      .equation-row {
        display: flex;
        align-items: center;
        gap: 4px;
        margin-bottom: 10px;
      }

      .equation-row .matrix-container {
        margin: 0;
        justify-content: flex-start;
      }

      .equation-label {
        font-size: 15px;
        font-weight: 600;
        min-width: 1.3em;
        text-align: right;
      }

      .equation-equals {
        font-size: 15px;
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
        margin-left: 3px;
        margin-right: 3px;
      }

      /* =================================================
         MATRIX INPUT
         ================================================= */

      .matrix-input {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 2px 2px;
        margin-left: 3px;
        margin-right: 3px;
      }

      /* =================================================
         READ-ONLY MATRIX VALUE GRID (e.g. R display)
         ================================================= */

      .matrix-value-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 2px 2px;
        margin-left: 3px;
        margin-right: 3px;
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
        box-sizing: border-box;
        width: 55px;
        height: 30px;
        padding: 2px;
        text-align: center;
        font-size: 14px;
      }

      /* =================================================
         BRACKETS
         ================================================= */

      .matrix-bracket {
        position: relative;
        width: 7px;
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
        box-sizing: border-box;
        width: 42px;
        height: 28px;
        display: flex;
        align-items: center;
        justify-content: center;
        font-size: 13px;
      }

      /* =================================================
         THIN SECTION DIVIDER
         ================================================= */

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

  # SCENARIO 1: QUARTERLY (user-editable base scenario, set in sidebar) ----
  R1 <- reactive({
    matrix(
      c(input$R11, input$R12, input$R21, input$R22),
      nrow = 2, ncol = 2, byrow = TRUE
    )
  })

  D1 <- reactive(matrix(c(input$D1, input$D2), nrow = 2, ncol = 1))
  G1 <- reactive(matrix(c(input$G1, input$G2), nrow = 2, ncol = 1))

  AB1 <- reactive(compute_AB(R1(), D1(), G1()))

  # SCENARIO 2: BIANNUAL ----
  # Doubling the interval squares R. Note: for an upper-triangular R with
  # equal diagonal elements a and off-diagonal b, the off-diagonal of R^2 is
  # 2ab, which equals b when a = 0.5. With the defaults (a = 0.5, b = 0.4) the
  # cross-lagged effect therefore stays 0.4 by coincidence; it is not a bug.
  R2 <- reactive(R1() %*% R1())
  D2 <- reactive(D1())
  G2 <- reactive(G1() * 2)

  AB2 <- reactive(compute_AB(R2(), D2(), G2()))

  # SCENARIO 3: ANNUAL ----
  R3 <- reactive(R2() %*% R2())
  D3 <- reactive(D1())
  G3 <- reactive(G2() * 2)

  AB3 <- reactive(compute_AB(R3(), D3(), G3()))

  # SCENARIO 1: FULL DISPLAY (read-only, copied from sidebar inputs) ----
  output$scenario1_matrices <- renderUI({
    ab <- AB1()
    div(
      matrix_bracket_display("R<sub>Qu</sub>", R1()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("&Delta;<sub>Qu</sub>", D1()),
      vector_bracket_display("&Gamma;<sub>Qu</sub>", G1()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("A<sub>Qu</sub>", ab$A, highlight_negative = TRUE),
      vector_bracket_display("B<sub>Qu</sub>", ab$B, highlight_negative = TRUE)
    )
  })

  # SCENARIO 2: FULL DISPLAY ----
  output$scenario2_matrices <- renderUI({
    ab <- AB2()
    div(
      matrix_bracket_display("R<sub>Bi</sub>", R2()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("&Delta;<sub>Bi</sub>", D2()),
      vector_bracket_display("&Gamma;<sub>Bi</sub>", G2()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("A<sub>Bi</sub>", ab$A, highlight_negative = TRUE),
      vector_bracket_display("B<sub>Bi</sub>", ab$B, highlight_negative = TRUE)
    )
  })

  # SCENARIO 3: FULL DISPLAY ----
  output$scenario3_matrices <- renderUI({
    ab <- AB3()
    div(
      matrix_bracket_display("R<sub>An</sub>", R3()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("&Delta;<sub>An</sub>", D3()),
      vector_bracket_display("&Gamma;<sub>An</sub>", G3()),
      tags$hr(class = "thin-hr"),
      vector_bracket_display("A<sub>An</sub>", ab$A, highlight_negative = TRUE),
      vector_bracket_display("B<sub>An</sub>", ab$B, highlight_negative = TRUE)
    )
  })
}

# RUN APP ----
shinyApp(
  ui = ui,
  server = server
)