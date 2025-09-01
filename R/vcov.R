#'@export
vcov.poisreg <- function(object, ...) {
  Delta <- object$Delta
  V <- MASS::ginv(object$fit$hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- object$labels
  rownames(V) <- object$labels
  return(V)
}
#'@export
vcov.zapoisreg <- function(object, ...) {
  Delta <- object$Delta
  V <- MASS::ginv(object$fit$hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- c("psi", object$labels)
  rownames(V) <- c("psi", object$labels)
  return(V)
}
#'@export
vcov.zipoisreg <- function(object, ...) {
  Delta <- object$Delta
  V <- MASS::ginv(object$fit$hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- with(object, c(labels1, labels2))
  rownames(V) <- with(object, c(labels1, labels2))
  return(V)
}
