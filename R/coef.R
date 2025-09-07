#'@export
coef.cms <- function(object, ...) {
  coeffs <- object$fit$par[-1]
  names(coeffs) <- object$labels
  return(coeffs)
}
#'@export
coef.cmai <- function(object, ...) {
  coefs <- object$fit$par[-1]
  p <- object$p
  q <- object$q
  coeffs1 <- coefs[1:q]
  coeffs2 <- coefs[(q + 1):(q + p)]
  names(coeffs1) <- object$labels1
  names(coeffs2) <- object$labels2
  coeffs1
  coeffs2
  coeffs <- list("Degenerated dist." = coeffs1, "Bell dist." = coeffs2)
  return(coeffs)
}
