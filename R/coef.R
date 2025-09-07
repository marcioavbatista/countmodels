#---------------------------------------------
#' S3 class to Estimated regression coefficients for countmodeg model
#'
#' @aliases coef.cms
#' @export
#' @param object an object of the class cms
#' @param ... further arguments passed to or from other methods
#' @return  a vector with the estimated regression coefficients.
#'
#' @examples
#'
#'
coef.cms <- function(object, ...) {
  coeffs <- object$fit$par[-1]
  names(coeffs) <- object$labels
  return(coeffs)
}
#---------------------------------------------
#' S3 class to Estimated regression coefficients for cmai models
#'
#' @aliases coef.cmai
#' @export
#' @param object an object of the class cmai
#' @param ... further arguments passed to or from other methods
#' @return  a list containing the the estimated regression coefficients associated with the degenerated and count distributions, respectively.
#'
#' @examples
#'
#'
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
