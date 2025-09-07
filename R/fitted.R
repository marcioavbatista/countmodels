#'@export
fitted.cms <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
}
#'@export
fitted.cmai <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  return(mu)
}
#'@export
residuals.cms <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  y <- model.response(mf)

  return(y - mu)
}
#'@export
residuals.cmai <- function(object, ...) {
  mf <- eval(object$call$data)
  mu <- predict(object, mf)[, 1]
  y <- model.response(mf)

  return(y - mu)
}
