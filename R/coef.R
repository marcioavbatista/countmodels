#'@export
coef.poisreg <- function(object, ...) {
  coeffs <- object$fit$par
  names(coeffs) <- object$labels
  return(coeffs)
}
#'@export
coef.zapoisreg <- function(object, ...) {
  coeffs <- object$fit$par
  names(coeffs) <- c("psi", object$labels)
  return(coeffs)
}
#'@export
coef.zipoisreg <- function(object, ...) {
  coefs <- object$fit$par
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
