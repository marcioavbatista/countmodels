#'@export
predict.cms <- function(
  object,
  newdata,
  level = 0.95,
  type = c("link", "response", "terms"),
  ...
) {
  link <- object$link
  labels <- object$labels
  betas <- object$fit$par[-1]
  type <- match.arg(type)
  formula <- object$formula
  f_rhs <- update(formula, NULL ~ .)
  vcov <- vcov(object)
  alpha <- 1 - level
  X <- model.matrix(f_rhs, newdata)

  if (!all(colnames(X) == labels)) {
    stop("newdata needs to have the same columnames of the model matrix")
  }

  if (type %in% c("response", "terms")) {
    warning("still in development, using link interval")
  }
  # forces link use
  type <- "link"

  pred_ <- X %*% betas

  # defines inverce functions

  inv_link <- function(pred, link) {
    if (link == "log") {
      return(exp(pred))
    } else if (link == "sqrt") {
      return(pred^2)
    } else if (link == "identity") {
      return(pred)
    } else {
      stop("error in link")
    }
  }

  vars <- diag(X %*% vcov %*% t(X))

  d <- stats::qnorm(1 - alpha / 2) * sqrt(vars)

  pred <- inv_link(pred_, link)
  lower <- inv_link(pred_ - d, link)
  upper <- inv_link(pred_ + d, link)

  CI <- cbind(pred = pred, lower = lower, upper = upper)

  return(CI)
}


#'@export
predict.cmai <- function(
  object,
  newdata,
  level = 0.95,
  type = c("link", "response", "terms"),
  ...
) {
  link1 <- object$link1
  link2 <- object$link2
  p <- object$p
  q <- object$q
  case <- object$case
  labels1 <- object$labels1
  labels2 <- object$labels2
  beta <- object$fit$par[2:(p + 1)]
  psi <- object$fit$par[(p + 2):(p + q + 1)]
  type <- match.arg(type)
  formula <- update(object$formula, NULL ~ .)
  Z <- stats::model.matrix(formula, data = newdata, rhs = 1)
  X <- stats::model.matrix(formula, data = newdata, rhs = 2)
  vcov <- vcov(object)
  vcovX <- vcov[1:p, 1:p]
  vcovZ <- vcov[(p + 1):(p + q), (p + 1):(p + q)]
  alpha <- 1 - level

  if (!all(colnames(X) == labels1) && !all(colnames(Z) == labels2)) {
    stop("newdata needs to have the same columnames of the model matrix")
  }

  if (type %in% c("response", "terms")) {
    warning("still in development, forcing use of link interval")
  }
  # forces link use
  type <- "link"

  predX <- X %*% beta
  predZ <- Z %*% psi

  # finding variance for each linear predictor

  covpredX <- diag(X %*% vcovX %*% t(X))
  covpredZ <- diag(Z %*% vcovZ %*% t(Z))

  d1 <- stats::qnorm(1 - alpha / 2) * sqrt(covpredX)
  d2 <- stats::qnorm(1 - alpha / 2) * sqrt(covpredZ)

  #finding the interval limits

  lowerX <- predX - d1
  upperX <- predX + d1
  lowerZ <- predZ - d2
  upperZ <- predZ + d2

  # defines inverce functions --

  inv_link1 <- function(pred, link) {
    if (link == "log") {
      return(exp(pred))
    } else if (link == "sqrt") {
      return(pred^2)
    } else if (link == "identity") {
      return(pred)
    } else {
      stop("error in link")
    }
  }

  inv_link2 <- function(pred, link) {
    if (link == "logit") {
      return(1 / (1 + exp(pred)))
    } else if (link == "probit") {
      return(pnorm(pred))
    } else if (link == "cloglog") {
      return(1 - exp(-exp(pred)))
    } else if (link == "cauchy") {
      return(pcauchy(pred))
    } else {
      stop("error in link")
    }
  }

  inv_mean <- function(mu1, mu2, case) {
    if (case == "inflated") {
      mu <- (1 - mu2) * mu1
    } else {
      mu <- ((1 - mu2) * mu1) / (1 - exp(-mu1))
    }
    return(mu)
  }

  #--

  mu1 <- inv_link1(predX, link = link1)
  mu2 <- inv_link2(predZ, link = link2)

  lowermu1 <- inv_link1(lowerX, link = link1)
  lowermu2 <- inv_link2(lowerZ, link = link2)

  uppermu1 <- inv_link1(upperX, link = link1)
  uppermu2 <- inv_link2(upperZ, link = link2)

  mu <- inv_mean(mu1, mu2, case = case)
  lowermu <- inv_mean(lowermu1, lowermu2, case = case)
  uppermu <- inv_mean(uppermu1, uppermu2, case = case)

  CI <- cbind(mu2, pred = mu, lower = lowermu, upper = uppermu)

  return(CI)
}
