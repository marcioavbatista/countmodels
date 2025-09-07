#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.cms
#' @export
#' @param object an object of the class cms
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required
#' @param ... further arguments passed to or from other methods
#' @return  A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#'
#' @examples
#----------------------------------------------

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
#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.cai
#' @export
#' @param object an object of the class cmai
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required
#' @param ... further arguments passed to or from other methods
#' @return  A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#'
#' @examples
#----------------------------------------------
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
