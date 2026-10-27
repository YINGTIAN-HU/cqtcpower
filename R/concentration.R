# ---------------------------------------------------------------------------
# Concentration generation
# ---------------------------------------------------------------------------

# Concentrations for the active subjects of one cohort: a T x n matrix.
.sim_conc_cohort <- function(conc, cohort, times, n) {
  if (n == 0) return(matrix(0, length(times), 0))
  if (conc$type == "pk_1cmt") {
    dose <- conc$dose[cohort]
    Ka <- conc$Ka * exp(rnorm(n, 0, conc$iiv_sd))
    CL <- conc$CL * exp(rnorm(n, 0, conc$iiv_sd))
    V  <- conc$V  * exp(rnorm(n, 0, conc$iiv_sd))
    Ke <- CL / V
    C <- vapply(seq_len(n), function(i) {
      dose * Ka[i] / (V[i] * (Ka[i] - Ke[i])) *
        (exp(-Ke[i] * times) - exp(-Ka[i] * times))
    }, numeric(length(times)))
    C <- matrix(C, nrow = length(times))
    if (conc$resid_cv > 0) {
      C <- C * exp(matrix(rnorm(length(C), 0, conc$resid_cv), nrow = nrow(C)))
    }
  } else {
    mu <- conc$cmax * conc$cmax_fraction[cohort] * conc$profile
    C <- matrix(mu, nrow = length(mu), ncol = n)
    if (conc$cv > 0) {
      C <- C + matrix(rnorm(length(C), 0, rep(mu * conc$cv, n)), nrow = nrow(C))
    }
  }
  C
}

# Concentrations for the whole study: a T x N matrix (subjects in columns,
# cohort by cohort, active subjects first) plus subject-level indicators.
.study_conc <- function(design, n_active, n_placebo) {
  T <- design$n_time
  blocks <- lapply(seq_len(design$n_cohort), function(g) {
    cbind(.sim_conc_cohort(design$conc, g, design$times, n_active),
          matrix(0, T, n_placebo))
  })
  list(C = do.call(cbind, blocks),
       treatment = rep(rep(c(1, 0), c(n_active, n_placebo)), design$n_cohort),
       cohort = rep(seq_len(design$n_cohort), each = n_active + n_placebo))
}

#' Theoretical concentration-time profile
#'
#' Mean concentration of each cohort at the design's ECG times (no
#' between-subject variability), and the maximum mean concentration. Useful for
#' choosing `C_supra`.
#'
#' @param design A [cqtc_design()].
#' @return A list with `conc` (time x cohort matrix) and `cmax` (maximum per
#'   cohort over the sampled times; for PK designs also the true `Cmax` and
#'   `Tmax` of the highest dose).
#' @examples
#' conc_typical(cqtc_design())$cmax
#' @export
conc_typical <- function(design) {
  stopifnot(inherits(design, "cqtc_design"))
  conc <- design$conc
  if (conc$type == "pk_1cmt") {
    conc$iiv_sd <- 0
    conc$resid_cv <- 0
  } else {
    conc$cv <- 0
  }
  M <- vapply(seq_len(design$n_cohort), function(g) {
    .sim_conc_cohort(conc, g, design$times, 1)[, 1]
  }, numeric(design$n_time))
  M <- matrix(M, nrow = design$n_time,
              dimnames = list(paste0("t", design$times), paste0("cohort", seq_len(design$n_cohort))))
  out <- list(conc = M, cmax = apply(M, 2, max))
  if (conc$type == "pk_1cmt") {
    Ke <- conc$CL / conc$V
    tmax <- (log(conc$Ka) - log(Ke)) / (conc$Ka - Ke)
    dose <- max(conc$dose)
    out$true_cmax <- dose * conc$Ka / (conc$V * (conc$Ka - Ke)) *
      (exp(-Ke * tmax) - exp(-conc$Ka * tmax))
    out$true_tmax <- tmax
  }
  out
}
