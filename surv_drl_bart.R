#' @title DR-learner with NFT-BART


surv_drl_bart <- function(X, Y, W, D, t0, W.hat = NULL, cen.fit = "nftbart",
                          surf1, surf0, C.hat, k.folds = 10, new.args.grf.nuisance = list()) {

  # set parameters for survival forest
  args.grf.nuisance <- list(failure.times = NULL,
                            num.trees = max(50, 2000 / 4),
                            sample.weights = NULL,
                            clusters = NULL,
                            equalize.cluster.weights = FALSE,
                            sample.fraction = 0.5,
                            mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                            min.node.size = 15,
                            honesty = TRUE,
                            honesty.fraction = 0.5,
                            honesty.prune.leaves = TRUE,
                            alpha = 0.05,
                            prediction.type = "Nelson-Aalen",
                            compute.oob.predictions = TRUE,
                            num.threads = NULL,
                            seed = runif(1, 0, .Machine$integer.max))
  args.grf.nuisance[names(new.args.grf.nuisance)] <- new.args.grf.nuisance

  Tgrf1 <- 1-surf1
  Tgrf0 <- 1-surf0

  # IPCW weights
  U <- pmin(Y, t0)                         # truncated follow-up time by t0
  if (cen.fit == "Kaplan-Meier") {
    fold.id <- sample(rep(seq(k.folds), length = nrow(X)))
    C.hat <- rep(NA, length(fold.id))
    for (z in 1:k.folds) {
      c.fit <- survival::survfit(survival::Surv(Y[!fold.id == z], 1 - D[!fold.id == z]) ~ 1)
      kmc <- summary(c.fit, times = U[fold.id == z])
      C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
    }
  } else if (cen.fit == "survival.forest") {
    args.grf.nuisance$alpha <- 0.05
    c.fit <- do.call(grf::survival_forest, c(list(X = cbind(W, X), Y = Y, D = 1 - D), args.grf.nuisance))
    C.hat <- predict(c.fit, failure.times = U, prediction.times = "time")$predictions
  } else if (cen.fit == "nftbart") {
    C.hat <- C.hat
  }
  if (any(C.hat == 0)) {
    stop("Some or all uncensored probabilities are exactly zeros. Check input variables or consider adjust the time of interest t0.")
  }

  # Propensity score
  if (is.null(W.hat)) {
    stop("propensity score needs to be supplied")
  } else if (length(W.hat) == 1) {
    W.hat <- rep(W.hat, length(W))
  } else if (length(W.hat) != nrow(X)) {
    stop("W.hat has incorrect length.")
  }

  # CATE function
  D.t0 <- D
  D.t0[D == 1 & Y > t0] <- 0
  D.t0 <- D.t0[D == 1 | Y > t0]
  W.t0 <- W[D == 1 | Y > t0]
  X.t0 <- X[D == 1 | Y > t0,, drop = FALSE]
  Tgrf0.t0 <- Tgrf0[D == 1 | Y > t0]
  Tgrf1.t0 <- Tgrf1[D == 1 | Y > t0]
  sample.weights.t0 <- 1 / C.hat[D == 1 | Y > t0]
  W.hat.t0 <- W.hat[D == 1 | Y > t0]

  Z <- W.t0 * D.t0 / W.hat.t0 - (1 - W.t0) * D.t0 / (1 - W.hat.t0)+
    Tgrf1.t0-W.t0 *Tgrf1.t0/W.hat.t0 - Tgrf0.t0 +
    ((1 - W.hat.t0) * Tgrf0.t0) / (1 - W.hat.t0)

  tau.fit <- grf::regression_forest(X.t0, Z, sample.weights = sample.weights.t0)
  tau.hat <- -predict(tau.fit, X)$predictions


  ret <- list(tau.fit = tau.fit,
              tau.hat = tau.hat)
  class(ret) <- "surv_drl_bart"
  ret
}

predict.surv_drl_bart <- function(object,
                                newdata = NULL,
                                ...) {
  if (!is.null(newdata)) {
    tau.hat <- -predict(object$tau.fit, newdata)$predictions
  } else {
    tau.hat <- object$tau.hat
  }
  return(tau.hat)
}

