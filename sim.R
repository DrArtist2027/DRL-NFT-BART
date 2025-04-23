
# *** Comparison methods ***
estimators <- list(cate_ml = cate_ml,
                   cate_xl_bart = cate_xl_bart,
                   cate_drl_bart = cate_drl_bart,
                   cate_drl_sf = cate_drl_sf,
                   cate_psdrl_bart = cate_psdrl_bart,
                   cate_psdrl_sf = cate_psdrl_sf,
                   cate_rl_bart = cate_rl_bart,
                   cate_nftbart_ps = cate_nftbart_ps,
                   cate_bart_ps = cate_bart_ps,
                   cate_csf_probs = cate_csf_probs)

  
myfun <- function(i, n = n, t0 = t0, e = e, overlap = overlap,
                  f.b=f.b, f.i = f.i, b.b = b.b, b.i = b.i,censor = censor,
                  dgps = dgps){

  set.seed(i)
  data.test = data_nft(n = 2000, t0 = t0, e = e, overlap = overlap,
                       f.b=f.b, f.i = f.i, b.b = b.b, b.i = b.i,censor = censor,
                       dgps = dgps, train = FALSE)
  # Generate data
  data = data_nft(n = n, t0 = t0, e = e, overlap = overlap,
                  f.b=f.b, f.i = f.i, b.b = b.b, b.i = b.i,censor = censor,
                  dgps = dgps, n.mc=1)
  
  true.catesp <- data.test$catesp
  true.catesp.sign <- data.test$catesp.sign

  # Store predictions
  predictions <- matrix(NA, nrow = 2000, ncol = length(estimators))
  estimator.output <- list()

  # Loop through estimators
  for (j in 1:length(estimators)) {
    estimator.name <- names(estimators)[j]
    message("Running estimator: ", estimator.name)

    # Get predictions for CATE
    predictions[, j] <- estimators[[estimator.name]](data, data.test, t0 = t0)

    # Compute performance metrics
    correct.classification <- sign(predictions[,j]) == true.catesp.sign

    # calibration slope
    calib.fit <- lm(predictions[, j] ~ true.catesp)

    # Save results for the estimator

    dfj <- data.frame(estimator.name = estimator.name,
                      true.catesp.var = var(true.catesp),
                      mse = mean((predictions[,j] - true.catesp)^2),
                      bias = mean(abs(predictions[,j] - true.catesp)),
                      rcorr = cor(predictions[,j], true.catesp),
                      taucorr = cor(predictions[,j], true.catesp, method = "kendall"),  # use Kendall's tau for concordance
                      calib.coef = calib.fit$coefficients[2],
                      classif.rate = mean(correct.classification, na.rm = TRUE) # NA: to ignore X1 < 0.3 in DGP 4.
    )
    dfj$rcorr[is.na(dfj$rcorr)==TRUE] <- 0  # assign correlation to 0 when CATE = ATE
    estimator.output[[j]] <- dfj
  }
  df <- do.call(rbind, estimator.output)
  df$n <- n
  df$censor <- censor
  df$horizon <- t0
  df$overlap <- overlap
  df$e <- e
  df$f.b <- f.b
  df$f.i <- f.i
  df$b.b <- b.b
  df$b.i <- b.i
  df$dgps <- dgps

  ret <- df
  return(ret)
}


