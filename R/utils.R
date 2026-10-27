# Evaluate `code` with a fixed seed and restore the caller's RNG state.
.with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  })
  set.seed(seed)
  code
}
