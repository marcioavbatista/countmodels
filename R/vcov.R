#---------------------------------------------
#' Variance-covariance matrix for a standard count model
#'
#' @aliases vcov.cms
#' @description This function extracts and returns the variance-covariance matrix associated with the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class cms.
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
vcov.cms <- function(object, ...) {
  Delta <- object$Delta
  hessian <- object$fit$hessian[-(object$p + 1), -(object$p + 1)]
  V <- MASS::ginv(hessian)
  V <- Delta %*% V %*% t(Delta)
  colnames(V) <- object$labels
  rownames(V) <- object$labels
  return(V)
}
#---------------------------------------------
#' Variance-covariance matrix for a inflated or hurdle count model
#'
#' @aliases vcov.cmai
#' @description This function extracts and returns the variance-covariance matrix associated with the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class cmai.
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
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
