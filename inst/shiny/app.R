# Shiny app for cqtcpower: analytical power, sample size and simulation for
# C-QTc analyses of parallel SAD studies. Launch with cqtcpower::cqtc_app().

library(shiny)
library(bslib)
library(cqtcpower)

# ---- helpers -------------------------------------------------------------------

# Parse "1, 3, 1/144, 0.5" into numbers (only digits, '.', '/', 'e', '-' allowed).
parse_nums <- function(txt, name) {
  parts <- trimws(strsplit(txt, "[,;[:space:]]+")[[1]])
  parts <- parts[nzchar(parts)]
  ok <- grepl("^[0-9.eE+-]+(/[0-9.eE+-]+)?$", parts)
  validate(need(length(parts) > 0 && all(ok),
                sprintf("%s: enter numbers separated by commas (fractions like 1/144 are fine).", name)))
  vapply(parts, function(p) {
    x <- strsplit(p, "/", fixed = TRUE)[[1]]
    v <- as.numeric(x)
    if (length(v) == 2) v[1] / v[2] else v
  }, numeric(1), USE.NAMES = FALSE)
}

model_cols <- c("RI" = "#1b6ca8", "RI + Z" = "#d1495b", "RS0" = "#3c8d2f",
                "RS0 + Z" = "#8e5ea2", "RS" = "#e08e0b", "RS + Z" = "#5c5c5c")

fmt_pct <- function(p) ifelse(is.na(p), "\u2013", sprintf("%.1f%%", 100 * p))

try_or_validate <- function(expr) {
  tryCatch(expr, error = function(e) validate(need(FALSE, conditionMessage(e))))
}

# ---- UI ------------------------------------------------------------------------

design_panel <- accordion_panel(
  "Study design", icon = icon("users"),
  radioButtons("conc_type", "Concentrations",
               c("One-compartment PK" = "pk", "Fixed profile" = "profile"), inline = TRUE),
  conditionalPanel(
    "input.conc_type == 'pk'",
    textInput("dose", "Doses (mg), one per cohort", "1, 3, 10, 30, 100, 200, 400, 800"),
    layout_column_wrap(
      width = 1/3,
      numericInput("Ka", "Ka (1/h)", 0.8, min = 0, step = 0.1),
      numericInput("CL", "CL", 12000, min = 0, step = 1000),
      numericInput("V", "V", 40000, min = 0, step = 1000)
    ),
    numericInput("iiv_sd", "Between-subject SD (log scale)", 0.1, min = 0, step = 0.05),
    textInput("times", "ECG times (h)", "1, 2, 4, 8, 12, 16, 24")
  ),
  conditionalPanel(
    "input.conc_type == 'profile'",
    numericInput("cmax", "Cmax of the top cohort", 0.008, min = 0, step = 0.001),
    textInput("frac", "Cohort Cmax as fraction of Cmax",
              "1/144, 1/48, 1/16, 1/8, 1/4, 1/2, 3/4, 1"),
    textInput("profile", "Time profile (fraction of cohort Cmax)",
              "0.45, 0.9, 0.95, 0.85, 0.75, 0.4, 0.1, 0.01")
  ),
  layout_column_wrap(
    width = 1/2,
    numericInput("n_active", "Active / cohort", 6, min = 1, step = 1),
    numericInput("n_placebo", "Placebo / cohort", 2, min = 1, step = 1)
  )
)

truth_panel <- accordion_panel(
  "True effect & variability", icon = icon("heart-pulse"),
  numericInput("U1", "True \u0394\u0394QTc at C_supra, U1 (ms)", 5, step = 0.5),
  layout_column_wrap(
    width = 1/2,
    numericInput("C_supra", "C_supra", 0.008, min = 0, step = 0.001),
    numericInput("theta_Z", "Treatment effect \u03b8Z (ms)", 0.5, step = 0.5)
  ),
  layout_column_wrap(
    width = 1/2,
    numericInput("sigma_e", "Residual SD \u03c3e", 5, min = 0.01, step = 0.5),
    numericInput("sigma0", "Random intercept SD \u03c30", 8, min = 0, step = 0.5)
  ),
  radioButtons("slope_mode", "Random slope given as",
               c("CV of slope" = "cv", "SD \u03c31" = "sd"), inline = TRUE),
  conditionalPanel("input.slope_mode == 'cv'",
                   numericInput("slope_cv", "CV of random slope", 0.75, min = 0, step = 0.05)),
  conditionalPanel("input.slope_mode == 'sd'",
                   numericInput("sigma1", "\u03c31 (ms per concentration unit)", 600, min = 0, step = 50)),
  numericInput("rho", "Intercept\u2013slope correlation \u03c1 (simulation only)", 0,
               min = -1, max = 1, step = 0.1)
)

rule_panel <- accordion_panel(
  "Decision rule", icon = icon("scale-balanced"),
  layout_column_wrap(
    width = 1/2,
    numericInput("U0", "Threshold U0 (ms)", 10, step = 1),
    numericInput("alpha", "One-sided \u03b1", 0.05, min = 0.001, max = 0.49, step = 0.01)
  ),
  sliderInput("target", "Target power", min = 0.5, max = 0.99, value = 0.8, step = 0.01)
)

ui <- page_sidebar(
  title = "C-QTc power and sample size (parallel SAD)",
  theme = bs_theme(version = 5, primary = "#1b6ca8", base_font = font_collection("system-ui", "sans-serif")),
  fillable = FALSE,
  sidebar = sidebar(
    width = 360,
    checkboxGroupInput("models", "Models (analytical)",
                       c("RI", "RI + Z", "RS0", "RS0 + Z"),
                       selected = c("RI", "RI + Z", "RS0", "RS0 + Z"), inline = TRUE),
    accordion(design_panel, truth_panel, rule_panel, open = "Study design", multiple = TRUE)
  ),
  navset_card_tab(
    id = "tabs",
    nav_panel(
      "Power",
      p(class = "text-muted small mb-1",
        "Analytical power at the design's cohort size. Green: at or above the target power."),
      uiOutput("power_boxes"),
      card(
        card_header("Power by cohort size"),
        sliderInput("n_max", "Show cohort sizes up to (subjects per cohort)",
                    min = 8, max = 160, value = 40, step = 4, width = "100%"),
        plotOutput("curve_plot", height = "380px"),
        downloadButton("dl_curve", "Download power curve (CSV)", class = "btn-sm")
      ),
      card(card_header("Details at the current design"), tableOutput("power_table"))
    ),
    nav_panel(
      "Sample size",
      card(
        card_header(textOutput("ss_header", inline = TRUE)),
        tableOutput("ss_table"),
        downloadButton("dl_ss", "Download (CSV)", class = "btn-sm")
      ),
      card(
        card_header("Subjects per cohort across true \u0394\u0394QTc"),
        sliderInput("u_range", "True \u0394\u0394QTc range (ms)", min = 0, max = 9, value = c(0, 5),
                    step = 1, width = "100%"),
        plotOutput("ss_grid_plot", height = "380px"),
        p(class = "text-muted small",
          "Each bar is the smallest cohort size whose analytical power reaches the target, ",
          "keeping the active:placebo ratio. With 'CV of slope', \u03c31 changes with U1.")
      )
    ),
    nav_panel(
      "Simulation",
      card(
        card_header("Check by simulation (fits mixed models with lmerTest)"),
        layout_column_wrap(
          width = 1/4,
          numericInput("nsim", "Simulated studies", 100, min = 10, max = 20000, step = 50),
          numericInput("sim_n", "Subjects per cohort", 8, min = 1, step = 1),
          numericInput("seed", "Seed", 1, step = 1),
          selectInput("ddf", "Degrees of freedom",
                      c("Wald", "Satterthwaite", "Kenward-Roger"),
                      selected = c("Wald", "Satterthwaite"), multiple = TRUE)
        ),
        checkboxGroupInput("sim_models", "Models",
                           names(model_cols), selected = c("RI + Z", "RS0 + Z", "RS + Z"),
                           inline = TRUE),
        p(class = "text-muted small",
          "Each simulated study fits every selected model; random-slope models and ",
          "Kenward-Roger are slower. Start with 100 studies (about a minute)."),
        actionButton("run_sim", "Run simulation", icon = icon("play"), class = "btn-primary")
      ),
      card(
        card_header("Results"),
        tableOutput("sim_table"),
        uiOutput("sim_note"),
        downloadButton("dl_sim", "Download (CSV)", class = "btn-sm")
      )
    ),
    nav_panel(
      "Concentrations",
      card(
        card_header("Typical concentration\u2013time profiles (no between-subject variability)"),
        plotOutput("conc_plot", height = "420px"),
        uiOutput("conc_note")
      )
    ),
    nav_panel(
      "About",
      card(markdown("
**What is computed.** A parallel SAD study concludes *no QT effect* when the upper bound of the
one-sided (1 \u2212 \u03b1) confidence interval of the predicted \u0394\u0394QTc at the supratherapeutic concentration,
U = \u03b8Z + \u03b8C\u00b7C_supra (or \u03b8C\u00b7C_supra without a treatment term), is below U0.
Analytical power is \u03a6((U0 \u2212 U1)/SE(\u00db) \u2212 z\u2081\u208b\u03b1), with SE(\u00db) from closed-form approximations
of the covariance of (\u03b8\u0302Z, \u03b8\u0302C) in linear mixed-effects models.

**Models.** All include fixed effects for baseline QTc and each time point.

| Model | Fixed effects | Random effects |
|---|---|---|
| RI | concentration | intercept |
| RI + Z | treatment + concentration | intercept |
| RS0 | concentration | intercept + slope, uncorrelated |
| RS0 + Z | treatment + concentration | intercept + slope, uncorrelated |
| RS | concentration | intercept + slope, correlated (simulation only) |
| RS + Z | treatment + concentration | intercept + slope, correlated (simulation only) |

**Notes.** RI models ignore any random slope in the data. The analytical formulas assume \u03c1 = 0;
\u03c1 only affects the simulation. Power does not depend on intercept, time or baseline effects, which
are fixed at the paper's values. Random PK concentrations are averaged over 1000 simulated studies.

Built on the `cqtcpower` R package.
"))
    )
  )
)

# ---- server --------------------------------------------------------------------

server <- function(input, output, session) {

  design_raw <- reactive({
    req(input$n_active, input$n_placebo)
    conc <- if (input$conc_type == "pk") {
      try_or_validate(conc_pk_1cmt(dose = parse_nums(input$dose, "Doses"),
                                   Ka = input$Ka, CL = input$CL, V = input$V,
                                   iiv_sd = input$iiv_sd))
    } else {
      try_or_validate(conc_profile(cmax = input$cmax,
                                   cmax_fraction = parse_nums(input$frac, "Cohort fractions"),
                                   profile = parse_nums(input$profile, "Time profile")))
    }
    times <- if (input$conc_type == "pk") parse_nums(input$times, "ECG times") else NULL
    try_or_validate(cqtc_design(conc, n_active = input$n_active,
                                n_placebo = input$n_placebo, times = times))
  })
  design <- debounce(design_raw, 500)

  params_for <- function(U1) {
    req(input$C_supra, input$sigma_e)
    slope <- if (input$slope_mode == "cv") list(slope_cv = input$slope_cv) else list(sigma1 = input$sigma1)
    try_or_validate(do.call(cqtc_params, c(list(
      U1 = U1, C_supra = input$C_supra, theta_Z = input$theta_Z,
      sigma_e = input$sigma_e, sigma0 = input$sigma0, rho = input$rho,
      U0 = input$U0, alpha = input$alpha), slope)))
  }
  params <- debounce(reactive(params_for(input$U1)), 500)

  models <- reactive({
    validate(need(length(input$models) > 0, "Select at least one model."))
    input$models
  })

  # ---- power ----
  curve <- reactive({
    d <- design()
    req(input$n_max)
    n_grid <- seq(d$block, max(d$block, input$n_max), by = d$block)
    try_or_validate(cqtc_power_curve(d, params(), models = models(),
                                     n_per_cohort = n_grid, target = input$target))
  })
  power_now <- reactive(try_or_validate(cqtc_power(design(), params(), models = models())))

  output$power_boxes <- renderUI({
    pw <- as.data.frame(unclass(power_now()))
    d <- design()
    boxes <- lapply(seq_len(nrow(pw)), function(i) {
      value_box(
        title = pw$model[i],
        value = fmt_pct(pw$power[i]),
        p(sprintf("%d per cohort (%d:%d), N = %d", pw$n_per_cohort[i], pw$n_active[i],
                  pw$n_placebo[i], pw$n_total[i])),
        theme = if (pw$power[i] >= input$target) "success" else "secondary"
      )
    })
    do.call(layout_column_wrap, c(list(width = 1 / min(4, length(boxes))), boxes))
  })

  output$curve_plot <- renderPlot({
    pc <- as.data.frame(unclass(curve()))
    op <- par(mar = c(4.5, 4.5, 1, 1)); on.exit(par(op))
    plot(range(pc$n_per_cohort), c(0, 1), type = "n", las = 1,
         xlab = "Subjects per cohort", ylab = "Power")
    abline(h = input$target, lty = 2, col = "grey50")
    abline(v = design()$n_active + design()$n_placebo, lty = 3, col = "grey70")
    for (m in unique(pc$model)) {
      s <- pc[pc$model == m, ]
      lines(s$n_per_cohort, s$power, col = model_cols[[m]], lwd = 2.5)
      points(s$n_per_cohort, s$power, col = model_cols[[m]], pch = 19, cex = 0.8)
    }
    legend("bottomright", legend = unique(pc$model), col = model_cols[unique(pc$model)],
           lwd = 2.5, bty = "n")
  }, res = 96)

  output$power_table <- renderTable({
    pw <- as.data.frame(unclass(power_now()))
    data.frame(Model = pw$model,
               `Per cohort (A:P)` = sprintf("%d (%d:%d)", pw$n_per_cohort, pw$n_active, pw$n_placebo),
               `Total N` = pw$n_total, `SE(U), ms` = sprintf("%.3f", pw$se_U),
               Power = fmt_pct(pw$power), check.names = FALSE)
  }, striped = TRUE, align = "l")

  output$dl_curve <- downloadHandler(
    filename = function() "cqtc_power_curve.csv",
    content = function(f) utils::write.csv(as.data.frame(unclass(curve())), f, row.names = FALSE))

  # ---- sample size ----
  ss <- reactive({
    validate(need(input$U1 < input$U0, "U1 must be below U0 for the power to reach a target."))
    try_or_validate(cqtc_sample_size(design(), params(), power = input$target, models = models()))
  })
  output$ss_header <- renderText(sprintf("Smallest cohort size for %.0f%% power", 100 * input$target))
  output$ss_table <- renderTable({
    s <- as.data.frame(unclass(ss()))
    data.frame(Model = s$model,
               `Per cohort (A:P)` = sprintf("%d (%d:%d)", s$n_per_cohort, s$n_active, s$n_placebo),
               `Total N` = s$n_total, `Power achieved` = fmt_pct(s$power),
               `SE(U), ms` = sprintf("%.3f", s$se_U), check.names = FALSE)
  }, striped = TRUE, align = "l")
  output$dl_ss <- downloadHandler(
    filename = function() "cqtc_sample_size.csv",
    content = function(f) utils::write.csv(as.data.frame(unclass(ss())), f, row.names = FALSE))

  ss_grid <- reactive({
    d <- design()
    req(input$u_range)
    u_seq <- seq(input$u_range[1], input$u_range[2], by = 1)
    u_seq <- u_seq[u_seq < input$U0]
    validate(need(length(u_seq) > 0, "Choose true ΔΔQTc values below U0."))
    cap <- 500
    rows <- lapply(u_seq, function(u) {
      s <- tryCatch(suppressWarnings(cqtc_sample_size(d, params_for(u), power = input$target,
                                                      models = models(),
                                                      max_n_per_cohort = cap)),
                    error = function(e) NULL)
      if (is.null(s)) return(NULL)
      data.frame(U1 = u, model = s$model, n = s$n_per_cohort,
                 capped = s$power < input$target)
    })
    do.call(rbind, rows)
  })
  observeEvent(input$U0, {
    mx <- max(1, ceiling(input$U0) - 1)
    updateSliderInput(session, "u_range", max = mx,
                      value = pmin(input$u_range, mx))
  })
  output$ss_grid_plot <- renderPlot({
    g <- ss_grid()
    validate(need(!is.null(g) && nrow(g) > 0, "No sample sizes to show."))
    ms <- unique(g$model)
    mat <- t(sapply(ms, function(m) g$n[g$model == m][match(sort(unique(g$U1)), g$U1[g$model == m])]))
    mat <- matrix(mat, nrow = length(ms), dimnames = list(ms, sort(unique(g$U1))))
    op <- par(mar = c(4.5, 4.5, 1, 1)); on.exit(par(op))
    bp <- barplot(mat, beside = TRUE, col = model_cols[ms], border = NA, las = 1,
                  xlab = "True \u0394\u0394QTc at C_supra, U1 (ms)", ylab = "Subjects per cohort",
                  ylim = c(0, max(mat, na.rm = TRUE) * 1.15))
    capm <- t(sapply(ms, function(m) g$capped[g$model == m][match(sort(unique(g$U1)), g$U1[g$model == m])]))
    capm <- matrix(capm, nrow = length(ms))
    lab <- ifelse(!is.na(capm) & capm, paste0("≥", mat), mat)
    text(bp, mat, labels = lab, pos = 3, cex = 0.75, col = "grey30")
    legend("topleft", legend = ms, fill = model_cols[ms], border = NA, bty = "n")
  }, res = 96)

  # ---- simulation ----
  sim_res <- eventReactive(input$run_sim, {
    d <- design()
    p <- params()
    ms <- input$sim_models
    validate(need(length(ms) > 0, "Select at least one model."),
             need(length(input$ddf) > 0, "Select at least one degrees-of-freedom method."))
    n <- input$sim_n
    validate(need(n %% d$block == 0,
                  sprintf("Subjects per cohort must be a multiple of %d (active:placebo %d:%d).",
                          d$block, d$n_active, d$n_placebo)))
    nsim <- input$nsim
    chunk <- max(1, ceiling(nsim / 25))
    starts <- seq(1, nsim, by = chunk)
    reps <- list()
    withProgress(message = "Simulating studies", value = 0, {
      for (k in seq_along(starts)) {
        m <- min(chunk, nsim - starts[k] + 1)
        s <- cqtc_simulate(d, p, nsim = m, models = ms, ddf = input$ddf, n_per_cohort = n,
                           seed = input$seed + starts[k] - 1, n_conc_sim = 10)
        if (!is.null(s$replicates)) {
          s$replicates$rep <- s$replicates$rep + starts[k] - 1
          reps[[k]] <- s$replicates
        }
        incProgress(m / nsim, detail = sprintf("%d of %d", starts[k] + m - 1, nsim))
      }
    })
    reps <- do.call(rbind, reps)
    ddf <- cqtcpower:::.match_ddf(input$ddf)
    summ <- cqtcpower:::.summarise_sim(reps, ms, ddf, p, nsim)
    an <- intersect(ms, c("RI", "RI + Z", "RS0", "RS0 + Z"))
    summ$analytical <- NA_real_
    if (length(an)) {
      ap <- cqtc_power(d, p, models = an, n_per_cohort = n)
      summ$analytical <- ap$power[match(summ$model, ap$model)]
    }
    list(summary = summ, nsim = nsim, n = n)
  })

  output$sim_table <- renderTable({
    s <- sim_res()$summary
    data.frame(Model = s$model, df = s$ddf,
               `Empirical power` = fmt_pct(s$power),
               `MC SE` = sprintf("%.1f%%", 100 * s$power_mcse),
               Analytical = fmt_pct(s$analytical),
               Coverage = fmt_pct(s$coverage),
               `Rel. bias` = ifelse(is.na(s$rel_bias), "\u2013", sprintf("%.1f%%", 100 * s$rel_bias)),
               `Not converged` = s$n_nonconv, Singular = s$n_singular, Failed = s$n_failed,
               check.names = FALSE)
  }, striped = TRUE, align = "l")
  output$sim_note <- renderUI({
    r <- sim_res()
    p(class = "text-muted small",
      sprintf("%d simulated studies with %d subjects per cohort. Coverage is for the two-sided %.0f%% interval. ",
              r$nsim, r$n, 100 * (1 - 2 * input$alpha)),
      "Singular = boundary fits (a variance component estimated as 0); these are valid estimates.")
  })
  output$dl_sim <- downloadHandler(
    filename = function() "cqtc_simulation.csv",
    content = function(f) utils::write.csv(sim_res()$summary, f, row.names = FALSE))

  observeEvent(design(), {
    d <- design()
    updateNumericInput(session, "sim_n", value = d$n_active + d$n_placebo,
                       step = d$block, min = d$block)
    updateSliderInput(session, "n_max", min = 2 * d$block, max = 40 * d$block, step = d$block,
                      value = min(max(input$n_max, 2 * d$block), 40 * d$block))
  })

  # ---- concentrations ----
  output$conc_plot <- renderPlot({
    d <- design()
    ct <- conc_typical(d)
    M <- ct$conc
    pos <- M[M > 0]
    validate(need(length(pos) > 0, "No positive concentrations."))
    cols <- grDevices::hcl.colors(ncol(M), "viridis")
    op <- par(mar = c(4.5, 5.5, 1, 1)); on.exit(par(op))
    matplot(d$times, pmax(M, min(pos) / 10), type = "b", pch = 19, lty = 1, col = cols, log = "y",
            las = 1, xlab = if (d$conc$type == "pk_1cmt") "Time (h)" else "Time point",
            ylab = "", ylim = range(c(pos, input$C_supra)))
    mtext("Concentration (log scale)", side = 2, line = 4)
    abline(h = input$C_supra, lty = 2, col = "#d1495b", lwd = 2)
    lab <- if (d$conc$type == "pk_1cmt") paste(d$conc$dose, "mg") else paste("cohort", seq_len(ncol(M)))
    legend("bottomleft", legend = c(lab, "C_supra"), col = c(cols, "#d1495b"),
           lty = c(rep(1, ncol(M)), 2), pch = c(rep(19, ncol(M)), NA),
           ncol = 3, bg = "white", box.col = "grey85", cex = 0.8)
  }, res = 96)
  output$conc_note <- renderUI({
    d <- design()
    ct <- conc_typical(d)
    txt <- sprintf("Highest sampled mean concentration: %.4g.", max(ct$cmax))
    if (!is.null(ct$true_cmax)) {
      txt <- paste(txt, sprintf("True Cmax of the top dose: %.4g at %.2f h.", ct$true_cmax, ct$true_tmax))
    }
    p(class = "text-muted small", txt, " The dashed line is C_supra.")
  })
}

shinyApp(ui, server)
