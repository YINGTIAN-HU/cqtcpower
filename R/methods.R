# ---------------------------------------------------------------------------
# print / plot methods
# ---------------------------------------------------------------------------

.header <- function(params, design) {
  cat(sprintf("Design: %d cohorts, active:placebo = %d:%d, %d ECG times\n",
              design$n_cohort, design$n_active, design$n_placebo, design$n_time))
  cat(sprintf("Truth:  U1 = %g ms at C_supra = %g (theta_Z = %g); sigma_e = %g, sigma0 = %g, sigma1 = %.4g\n",
              params$U1, params$C_supra, params$theta_Z, params$sigma_e,
              params$sigma0, params$sigma1))
  cat(sprintf("Rule:   one-sided %.0f%% upper bound < U0 = %g ms\n\n",
              100 * (1 - params$alpha), params$U0))
}

#' @export
print.cqtc_power <- function(x, digits = 3, ...) {
  cat("Analytical C-QTc power\n")
  .header(attr(x, "params"), attr(x, "design"))
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df$allocation <- sprintf("%d (%d:%d)", df$n_per_cohort, df$n_active, df$n_placebo)
  out <- data.frame(model = df$model, `n/cohort (A:P)` = df$allocation,
                    `total N` = df$n_total, `SE(U)` = round(df$se_U, digits),
                    power = round(df$power, digits),
                    check.names = FALSE, stringsAsFactors = FALSE)
  print(out, row.names = FALSE)
  invisible(x)
}

#' @export
print.cqtc_sample_size <- function(x, digits = 3, ...) {
  method <- attr(x, "method")
  cat(sprintf("C-QTc sample size (%s), target power %.0f%%\n",
              method, 100 * attr(x, "target")))
  .header(attr(x, "params"), attr(x, "design"))
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  out <- data.frame(model = df$model,
                    `n/cohort (A:P)` = sprintf("%d (%d:%d)", df$n_per_cohort, df$n_active, df$n_placebo),
                    `total N` = df$n_total, power = round(df$power, digits),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (method == "analytical") {
    out$`SE(U)` <- round(df$se_U, digits)
    if (any(df$n_formula != df$n_per_cohort)) out$`n/cohort (closed form)` <- df$n_formula
  } else {
    out$`MC SE` <- round(df$power_mcse, digits)
    cat(sprintf("Empirical power from %d simulated studies per step, %s df.\n\n",
                attr(x, "nsim"), df$ddf[1]))
  }
  print(out, row.names = FALSE)
  invisible(x)
}

#' @export
print.cqtc_sim <- function(x, digits = 3, ...) {
  cat(sprintf("Simulated C-QTc operating characteristics (%d studies, %d subjects per cohort)\n",
              x$nsim, x$n_per_cohort))
  .header(x$params, x$design)
  s <- x$summary
  out <- data.frame(model = s$model, ddf = s$ddf,
                    power = round(s$power, digits),
                    `MC SE` = round(s$power_mcse, digits),
                    analytical = round(s$analytical, digits),
                    coverage = round(s$coverage, digits),
                    `rel. bias` = round(s$rel_bias, digits),
                    nonconv = s$n_nonconv, singular = s$n_singular,
                    failed = s$n_failed,
                    check.names = FALSE, stringsAsFactors = FALSE)
  print(out, row.names = FALSE)
  cat("\nnonconv: fits with a convergence warning after trying all optimisers;",
      "singular: boundary fits (a variance component estimated as 0).\n")
  invisible(x)
}

#' Plot power against cohort size
#'
#' @param x A `cqtc_power` object, typically from [cqtc_power_curve()].
#' @param target Target power drawn as a dashed line (`NULL` for none).
#' @param xvar Plot against `"n_per_cohort"` or `"n_total"`.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.cqtc_power <- function(x, target = attr(x, "target"),
                            xvar = c("n_per_cohort", "n_total"), ...) {
  xvar <- match.arg(xvar)
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  models <- unique(df$model)
  cols <- c("#1b6ca8", "#d1495b", "#3c8d2f", "#8e5ea2")[seq_along(models)]
  graphics::plot(range(df[[xvar]]), c(0, 1), type = "n",
                 xlab = if (xvar == "n_per_cohort") "Subjects per cohort" else "Total subjects",
                 ylab = "Power", las = 1, ...)
  if (!is.null(target)) graphics::abline(h = target, lty = 2, col = "grey50")
  for (i in seq_along(models)) {
    d <- df[df$model == models[i], ]
    d <- d[order(d[[xvar]]), ]
    graphics::lines(d[[xvar]], d$power, col = cols[i], lwd = 2)
    graphics::points(d[[xvar]], d$power, col = cols[i], pch = 19 + (i - 1) %% 6)
  }
  graphics::legend("bottomright", legend = models, col = cols, lwd = 2,
                   pch = 19 + (seq_along(models) - 1) %% 6, bty = "n")
  invisible(x)
}

#' @export
as.data.frame.cqtc_sim <- function(x, ...) x$summary
