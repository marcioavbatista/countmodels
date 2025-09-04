#'@export
countmodel <- function(
  formula,
  family = c("poisson", "negbinom"),
  case = c("standard", "inflated", "hurdle"),
  data = NULL,
  approach = c("mle", "bayes"),
  hessian = TRUE,
  link1 = c("log", "sqrt", "identity"),
  link2 = c("log", "sqrt", "identity"),
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
  hyperpars_aux <- function(hyperpars, simp = F) {
    if (simp) {
      for (nome in names(hyperpars)) {
        if (startsWith(nome, "mu_beta") && !is.null(hyperpars[[nome]])) {
          hyperpars[[nome]] <- rep(0, p) # vetor de zeros
        } else if (
          startsWith(nome, "sigma_beta") && !is.null(hyperpars[[nome]])
        ) {
          hyperpars[[nome]] <- diag(10, p, p) # matriz pxp com diagonal 10
        } else if (
          (startsWith(nome, "a") || startsWith(nome, "b")) &&
            !is.null(hyperpars[[nome]])
        ) {
          hyperpars[[nome]] <- 0.01 # valor numérico
        }
      }
      return(hyperpars)
    } else {
      for (nome in names(hyperpars)) {
        if (startsWith(nome, "mu") && !is.null(hyperpars[[nome]])) {
          hyperpars[[nome]] <- rep(0, p) # vetor de zeros
        } else if (
          startsWith(nome, "sigma_psi") && !is.null(hyperpars[[nome]])
        ) {
          hyperpars[[nome]] <- diag(10, q, q) # matriz pxp com diagonal 10
        } else if (
          startsWith(nome, "sigma_beta") && !is.null(hyperpars[[nome]])
        ) {
          hyperpars[[nome]] <- diag(10, p, p) # matriz pxp com diagonal 10
        } else if (
          (startsWith(nome, "a") || startsWith(nome, "b")) &&
            !is.null(hyperpars[[nome]])
        ) {
          hyperpars[[nome]] <- 0.01 # valor numérico
        }
      }

      return(hyperpars)
    }
  }

  approach <- match.arg(approach)
  family <- match.arg(family)
  case <- match.arg(case)

  if (case == "standard") {
    link1 <- match.arg(link1)
    mf <- stats::model.frame(formula = formula, data = data)
    Terms <- stats::terms(mf)
    X <- as.matrix(stats::model.matrix(attr(mf, "terms"), data = mf))
    labels <- colnames(X)
    y <- stats::model.response(mf)
    n <- nrow(X)
    p <- ncol(X)
    hyperpars <- hyperpars_aux(hyperpars, simp = T)

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

    Link <- switch(link1, "log" = 1, "sqrt" = 2, "identity" = 3)

    if (family == "poisson") {
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
        dist = 1
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
        fit$par <- fit$theta_tilde[((p + 1):(2 * p))]
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
      fit$link <- link1
      class(fit) <- "cmpois"
      return(fit)
    } else {
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
        dist = 2
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
        fit$par <- fit$theta_tilde[((p + 1):(2 * p))]
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
      fit$link <- link1
      class(fit) <- "cmnb"
      return(fit)
    }
  } else {
    link1 <- match.arg(link1)
    link2 <- match.arg(link2)
    formula <- Formula::Formula(formula)
    mf <- stats::model.frame(formula = formula, data = data)
    Terms <- stats::terms(mf)
    Z <- stats::model.matrix(formula, data = mf, rhs = 1)
    X <- stats::model.matrix(formula, data = mf, rhs = 2)
    Xlabels <- colnames(X)
    Zlabels <- colnames(Z)
    y <- stats::model.response(mf)
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Z)
    hyperpars <- hyperpars_aux(hyperpars)
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

    Link1 <- switch(
      link1,
      "logit" = 1,
      "probit" = 2,
      "cloglog" = 3,
      "cauchy" = 4
    )

    Link2 <- switch(link2, "log" = 1, "sqrt" = 2, "identity" = 3)

    if (case == "inflated") {
      if (family == "poisson") {
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
          link1 = Link1,
          link2 = Link2,
          dist = 1,
          case = 1
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
          fit$par <- fit$theta_tilde[(p + q + 1):(2 * (p + q))]
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
        class(fit) <- "cmzipois"
        return(fit)
      } else {
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
          link1 = Link1,
          link2 = Link2,
          dist = 2,
          case = 1
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
          fit$par <- fit$theta_tilde[(p + q + 1):(2 * (p + q))]
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
        class(fit) <- "cmzinb"
        return(fit)
      }
    } else {
      if (family == "poisson") {
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
          link1 = Link1,
          link2 = Link2,
          dist = 1,
          case = 2
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
          fit$par <- fit$theta_tilde[(p + q + 1):(2 * (p + q))]
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
        class(fit) <- "cmzapois"
        return(fit)
      } else {
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
          link1 = Link1,
          link2 = Link2,
          dist = 2,
          case = 2
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
          fit$par <- fit$theta_tilde[(p + q + 1):(2 * (p + q))]
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
        class(fit) <- "cmzanb"
        return(fit)
      }
    }
  }
}
