#--------------------------------------------------
#' Regression for count models
#' @export
#' @aliases countmodel
#' @description
#' Fits an regression model for poisson and  Negative Binomial an their inflated
#'  and hurdle versions
#'
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted.
#' @param family family distribution for the regression model (poisson, negbinom)
#' ; default is poisson
#' @param case case of the regression (standard, inflated, hurdle); default is
#'  standard
#' @param data data an optional data.frame, list or enviroment containing the
#' variables in the model
#' @param approach approach to be used to fit the model (mle: maximum likelihood;
#'  bayes : Bayesian approach); default is mle
#' @param hessian logical ; if TRUE (default), the hessian matrix is returned
#' when approach = "mle"
#' @param link2 assumed link function for degenerate distribution (logit, probit,
#'  cloglog, cauchy); default is logit.
#' @param link1 assumed link function for count distribution (log, sqrt or
#' identiy); default is log.
#' @param hyperparsa list containing the hyperparameters associated with the
#' prior distribution of the regression coefficients;
#'  if not specified then default choice is an minor informative prior
#' @param ... further arguments
#'
#' @returns an object of class (cms,cmai,cmps,cmnb,cmzips,cmzaps,cmzinb,cmzanb)
#' depending on the family and case used
#'
#'
#' @examples
#'
#'
#'
#--------------------------------------------------
countmodel <- function(
  formula,
  family = c("poisson", "negbinom"),
  case = c("standard", "inflated", "hurdle"),
  data = NULL,
  approach = c("mle", "bayes"),
  hessian = TRUE,
  link1 = c("log", "sqrt", "identity"),
  link2 = c("logit", "probit", "cloglog", "cauchy"),
  hyperpars = list(
    mu_psi = NULL,
    sigma_psi = NULL,
    mu_beta = NULL,
    sigma_beta = NULL,
    a_theta = NULL,
    b_theta = NULL
  ),
  ...
) {
  hyperpars_aux <- function(hyperpars, p = NULL, q = NULL, simp = FALSE) {
    hp <- hyperpars # cópia local
    names_hp <- names(hp)

    # defaults comuns
    default_a_b <- 0.01
    default_sigma_diag_p <- function() diag(10, p, p)
    default_sigma_diag_q <- function() diag(10, q, q)
    default_mu_beta <- function() rep(0, p)
    default_mu_psi <- function() rep(0, q)

    if ("mu_beta" %in% names_hp && is.null(hp$mu_beta)) {
      hp$mu_beta <- default_mu_beta()
    }
    if ("sigma_beta" %in% names_hp && is.null(hp$sigma_beta)) {
      hp$sigma_beta <- default_sigma_diag_p()
    }
    if ("a_theta" %in% names_hp && is.null(hp$a_theta)) {
      hp$a_theta <- default_a_b
    }
    if ("b_theta" %in% names_hp && is.null(hp$b_theta)) {
      hp$b_theta <- default_a_b
    }

    # quando existe psi (modelos zereg)
    if (!simp) {
      if ("mu_psi" %in% names_hp && is.null(hp$mu_psi)) {
        hp$mu_psi <- default_mu_psi()
      }
      if ("sigma_psi" %in% names_hp && is.null(hp$sigma_psi)) {
        hp$sigma_psi <- default_sigma_diag_q()
      }
    }

    return(hp)
  }
  case <- match.arg(case)
  dist <- switch(family, "poisson" = 1, "negbinom" = 2)
  case_id <- switch(case, "inflated" = 1, "hurdle" = 2)
  approach <- match.arg(approach)
  family <- match.arg(family)

  if (case == "standard") {
    classr <- switch(family, "poisson" = "cmps", "negbinom" = "cmnb")

    link <- match.arg(link1)
    mf <- stats::model.frame(formula = formula, data = data)
    Terms <- stats::terms(mf)
    X <- as.matrix(stats::model.matrix(attr(mf, "terms"), data = mf))
    labels <- colnames(X)
    y <- stats::model.response(mf)
    n <- nrow(X)
    p <- ncol(X)
    hyperpars <- hyperpars_aux(hyperpars, p = p, simp = T)

    if (match("(Intercept)", labels) == 1) {
      X_std <- scale(X[, -1])
      x_mean <- array(c(0, attr(X_std, "scaled:center")))
      x_sd <- array(c(1, attr(X_std, "scaled:scale")))
      X_std <- cbind(1, X_std)
      Delta <- diag(1 / x_sd)
      Delta[1, ] <- Delta[1, ] - x_mean / x_sd
    } else {
      X_std <- scale(X)
      x_mean <- array(attr(X_std, "scaled:center"))
      x_sd <- array(attr(X_std, "scaled:scale"))
      Delta <- diag(1 / x_sd)
    }

    Link <- switch(link, "log" = 1, "sqrt" = 2, "identity" = 3)

    stan_data <- list(
      y = y,
      X = X_std,
      n = n,
      p = p,
      x_mean = x_mean,
      x_sd = x_sd,
      mu_beta = hyperpars$mu_beta,
      sigma_beta = hyperpars$sigma_beta,
      a_theta = hyperpars$a_theta,
      b_theta = hyperpars$b_theta,
      approach = 0,
      link = Link,
      dist = dist
    )

    if (approach == "mle") {
      fit <- rstan::optimizing(
        stanmodels$reg,
        hessian = hessian,
        data = stan_data,
        verbose = FALSE,
        ...
      )
      if (hessian == TRUE) {
        fit$hessian <- -fit$hessian
      }
      fit$par <- fit$theta_tilde[((p + 1):(2 * p + 1))]
      AIC <- -2 * fit$value + 2 * p
      fit <- list(fit = fit, logLik = fit$value, AIC = AIC, Delta = Delta)
    } else {
      stan_data$approach <- 1
      fit <- rstan::sampling(
        stanmodels$reg,
        data = stan_data,
        verbose = FALSE,
        ...
      )
      fit <- list(fit = fit)
    }

    fit$mf <- mf
    fit$n <- n
    fit$p <- p
    # fit$x_mean <- x_mean
    # fit$x_sd <- x_sd

    fit$call <- match.call()
    fit$formula <- stats::formula(Terms)
    fit$terms <- stats::terms.formula(formula)
    fit$labels <- labels
    fit$approach <- approach
    fit$link <- link
    class(fit) <- c("cms", classr)
    return(fit)
  } else {
    if (family == "poisson") {
      if (case == "inflated") {
        classr <- "cmzips"
      } else {
        classr <- "cmzaps"
      }
    } else {
      if (case == "inflated") {
        classr <- "cmzinb"
      } else {
        classr <- "cmzanb"
      }
    }

    link1 <- match.arg(link1)
    link2 <- match.arg(link2)
    formula <- Formula::Formula(formula)
    mf <- stats::model.frame(formula = formula, data = data)
    Terms <- stats::terms(mf)
    Z <- stats::model.matrix(formula, data = mf, rhs = 2)
    X <- stats::model.matrix(formula, data = mf, rhs = 1)
    Xlabels <- colnames(X)
    Zlabels <- colnames(Z)
    y <- stats::model.response(mf)
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Z)
    hyperpars <- hyperpars_aux(hyperpars, p = p, q = q)
    if (p > 1) {
      if (match("(Intercept)", Xlabels) == 1) {
        X_std <- scale(X[, -1])
        x_mean <- array(c(0, attr(X_std, "scaled:center")))
        x_sd <- array(c(1, attr(X_std, "scaled:scale")))
        X_std <- cbind(1, X_std)
        Delta_x <- diag(1 / x_sd)
        Delta_x[1, ] <- Delta_x[1, ] - x_mean / x_sd
      } else {
        X_std <- scale(X)
        x_mean <- array(attr(X_std, "scaled:center"))
        x_sd <- array(attr(X_std, "scaled:scale"))
        Delta_x <- diag(1 / x_sd)
      }
    } else {
      X_std <- X
      x_mean <- array(0)
      x_sd <- array(1)
      Delta_x <- diag(n)
    }

    if (q > 1) {
      if (match("(Intercept)", Zlabels) == 1) {
        Z_std <- scale(Z[, -1])
        z_mean <- array(c(0, attr(Z_std, "scaled:center")))
        z_sd <- array(c(1, attr(Z_std, "scaled:scale")))
        Z_std <- cbind(1, Z_std)
        Delta_z <- diag(1 / z_sd)
        Delta_z[1, ] <- Delta_z[1, ] - z_mean / z_sd
      } else {
        Z_std <- scale(Z)
        z_mean <- array(attr(Z_std, "scaled:center"))
        z_sd <- array(attr(Z_std, "scaled:scale"))
        Delta_z <- diag(1 / z_sd)
      }
    } else {
      Z_std <- Z
      z_mean <- array(0)
      z_sd <- array(1)
      Delta_z <- diag(n)
    }

    Link2 <- switch(
      link2,
      "logit" = 1,
      "probit" = 2,
      "cloglog" = 3,
      "cauchy" = 4
    )

    Link1 <- switch(link1, "log" = 1, "sqrt" = 2, "identity" = 3)

    stan_data <- list(
      y = y,
      X = X_std,
      Z = Z_std,
      n = n,
      p = p,
      q = q,
      x_mean = x_mean,
      x_sd = x_sd,
      z_mean = z_mean,
      z_sd = z_sd,
      mu_beta = hyperpars$mu_beta,
      sigma_beta = hyperpars$sigma_beta,
      mu_psi = hyperpars$mu_psi,
      sigma_psi = hyperpars$sigma_psi,
      a_theta = hyperpars$a_theta,
      b_theta = hyperpars$b_theta,
      approach = 0,
      link1 = Link2,
      link2 = Link1,
      dist = dist,
      case_id = case_id
    )

    if (approach == "mle") {
      fit <- rstan::optimizing(
        stanmodels$zereg,
        hessian = hessian,
        data = stan_data,
        verbose = FALSE,
        ...
      )
      if (hessian == TRUE) {
        fit$hessian <- -fit$hessian
      }
      fit$par <- fit$theta_tilde[(p + q + 1):(2 * (p + q) + 1)]
      AIC <- -2 * fit$value + 2 * (p + q)
      fit <- list(
        fit = fit,
        logLik = fit$value,
        AIC = AIC,
        Delta = magic::adiag(Delta_z, Delta_x)
      )
    } else {
      stan_data$approach <- 1
      fit <- rstan::sampling(
        stanmodels$zereg,
        data = stan_data,
        verbose = FALSE,
        ...
      )
      fit <- list(fit = fit)
    }

    fit$n <- n
    fit$p <- p
    fit$q <- q
    # fit$x_mean <- x_mean
    # fit$x_sd <- x_sd
    # fit$z_mean <- z_mean
    # fit$z_sd <- z_sd
    # fit$v_sd <- c(z_sd, x_sd)

    fit$call <- match.call()
    fit$formula <- stats::formula(Terms)
    fit$terms <- stats::terms.formula(formula)
    fit$labels1 <- Zlabels
    fit$labels2 <- Xlabels
    fit$approach <- approach
    fit$link1 <- link1
    fit$link2 <- link2
    fit$case <- case
    class(fit) <- c("cmai", classr)
    return(fit)
  }
}
