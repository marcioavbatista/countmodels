#'@export
vcov.cms <- function(object, ...) {
  Delta <- object$Delta
  hessian <- object$fit$hessian[-(object$p + 1), -(object$p + 1)]
  V <- MASS::ginv(hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- object$labels
  rownames(V) <- object$labels
  return(V)
}
#'@export
vcov.cmai <- function(object, ...) {
  Delta <- object$Delta
  hessian <- object$fit$hessian[
    -(object$p + object$q + 1),
    -(object$p + object$q + 1)
  ]
  V <- MASS::ginv(hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- with(object, c(labels1, labels2))
  rownames(V) <- with(object, c(labels1, labels2))
  return(V)
}
