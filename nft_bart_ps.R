

nft_bart_ps <- function(X, Y, W, D, t0, W.hat = NULL, x.test,
                         new.args.grf.nuisance = list()) {

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
  
  W.hat = ifelse(W == 0, 1-W.hat, W.hat)
  ps = ifelse(x.test$W == 0, 1-x.test$W.hat, x.test$W.hat)
  
  # fit model on W == 1
  X1<-cbind(X[W == 1,, drop = FALSE],W.hat[W == 1])
  X0<-cbind(X[W == 0,, drop = FALSE],W.hat[W == 0])
  grffit1 <- nftbart::nft2(X1,
                           X1,
                           Y[W == 1],
                           D[W == 1], K=0)
  surf1 <- predict(grffit1,cbind(x.test$X,ps),cbind(x.test$X,ps))
  
  # fit model on W == 0
  grffit0 <- nftbart::nft2(X0,
                           X0,
                           Y[W == 0],
                           D[W == 0], K=0)
  surf0 <- predict(grffit0,cbind(x.test$X,ps),cbind(x.test$X,ps))
  
  S1.hat <- apply(exp(surf1[["f.test"]]), 2, function(x) mean(x > t0))
  S0.hat <- apply(exp(surf0[["f.test"]]), 2, function(x) mean(x > t0))
  
  Tgrf1 <- S1.hat
  Tgrf0 <- S0.hat

  tau.hat <- Tgrf1 - Tgrf0

  ret <- list(tau.hat = tau.hat)

  class(ret) <- "nft_bart_ps"
  ret
}


