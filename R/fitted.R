#---------------------------------------------
#' Predictions for mean for each observation
#'
#' @aliases fitted.cms
#' @export
#' @param object an object of the class cms
#' @param ... further arguments passed to or from other methods
#' @return  vector of the prediction for each observation
#'
#' @examples
#----------------------------------------------
fitted.cms <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
}
#---------------------------------------------
#' Predictions for mean for each observation
#'
#' @aliases fitted.cmai
#' @export
#' @param object an object of the class cms
#' @param ... further arguments passed to or from other methods
#' @return  vector of the prediction for each observation
#'
#' @examples
#----------------------------------------------
fitted.cmai <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  return(mu)
}
#---------------------------------------------
#' Residuals for mean for each observation
#'
#' @aliases residuals.cms
#' @export
#' @param object an object of the class cms
#' @param ... further arguments passed to or from other methods
#' @return  vector of the prediction for each observation
#'
#' @examples
#----------------------------------------------
residuals.cms <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  y <- model.response(mf)

  return(y - mu)
}
#---------------------------------------------
#' Residuals for mean for each observation
#'
#' @aliases residuals.cmai
#' @export
#' @param object an object of the class cms
#' @param ... further arguments passed to or from other methods
#' @return  vector of the prediction for each observation
#'
#' @examples
#----------------------------------------------
residuals.cmai <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  y <- model.response(mf)

  return(y - mu)
}
