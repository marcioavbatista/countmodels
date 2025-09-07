#'@export
confint.cms <- function(object, parm = NULL, level = 0.95, ...) {
  V <- vcov(object)
  par.hat <- object$fit$par[-1]
  alpha <- 1 - level
  d <- stats::qnorm(1 - alpha / 2) * sqrt(diag(V))
  lower <- par.hat - d
  upper <- par.hat + d
  CI <- cbind(lower, upper)
  labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
  colnames(CI) <- paste0(labels, "%")
  if (is.null(parm)) {
    return(CI)
  } else {
    CI <- CI[parm, , drop = FALSE]
    return(CI)
  }
}
#'@export
confint.cmai <- function(object, parm = NULL, level = 0.95, ...) {
  p <- object$p
  q <- object$q
  V <- vcov(object)
  par.hat <- object$fit$par[-1]
  alpha <- 1 - level
  d <- stats::qnorm(1 - alpha / 2) * sqrt(diag(V))
  lower <- par.hat - d
  upper <- par.hat + d
  ci <- cbind(lower, upper)
  labels <- round(100 * (c(alpha / 2, 1 - alpha / 2)), 1)
  colnames(ci) <- paste0(labels, "%")
  # if(!is.null(parm)){
  #   ci <- ci[parm, ,drop = FALSE]
  # }
  CI <- list(
    "Degenerated dist." = ci[1:p, ],
    "Bell dist." = ci[(q + 1):(q + p), ]
  )
  return(CI)
}
