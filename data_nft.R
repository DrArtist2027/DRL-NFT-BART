library(survival)
library(rms)
library(mvtnorm)
library(SuperLearner)

# data generation function
data_nft = function(n=200, t0 = 0.5, e = 0.5, overlap = "strong",
                    f.b="L", f.i = "L", b.b = "Y", b.i = "N",censor = "25%",
                    dgps = "nft",n.mc=10000, train = TRUE) {

  # generate continuous covariates X1~X8
  rho=0.5
  R <- (1-rho) * diag(8) + rho * matrix(1,8,8)
  Xn <- rmvnorm(n, rep(0,8), R)
  Xc <- Xn[,1:3]


  # generate binary covariates X9~X15
  Xi = cbind(rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5),
             rbinom(n=n,size=1,p=0.5))
  Xb = Xi[,1:3]

  X=cbind(Xn,Xi)

  # regression coefficients
  betac = c(0.2, 0.3, 0.3)
  betab = c(-0.2,-0.3,-0.2)
  
  if (e==0.5) {
    # propensity score
    if (overlap=="strong") {
      e = c(plogis(0.4 + Xc %*% c(1*betac) + Xb %*% c(1*betab)))
    }
    
    if (overlap=="medium") {
      e = c(plogis(1.1 + Xc %*% c(3.0*betac) + Xb %*% c(3.0*betab)))
    }
    
    if (overlap=="weak") {
      e = c(plogis(1.8 + Xc %*% c(5*betac) + Xb %*% c(5*betab)))
    }
    # generate intervention status
    z <- rbinom(n, 1, e)
  }else if(e==0.1) {
    e = c(plogis(-0.8 + Xc %*% c(5*betac) + Xb %*% c(5*betab)))
    z <- rbinom(n, 1, e)
  }
  
  

  # generate T
  sigma <- exp(-2+0.8*X[,12]+1.6*X[,5]-1.8*X[,5]*X[,12])
  LPc = 0.5*X[,5] + 0.3*X[,13]
  scalec = 1
  if (f.b == "L" & f.i == "L") {
    if (b.b == "Y" & b.i == "N") {
      if (censor == "25%" & dgps == "nft"){
        ft <- exp(-0.3 + c(Xc %*% c(0.3,0.2,0.1)) + c(Xb %*% c(-0.3,-0.3,-0.2))+
                    (-0.5+X[,4]) * z + rnorm(n, mean = 0, sd = sigma))
        C = 0.35+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
      }else if (censor == "25%" & dgps == "aft"){
        ft <- exp(-0.6 + c(Xc %*% c(0.3,0.2,0.1)) + c(Xb %*% c(-0.3,-0.3,-0.2))+
                    (-0.5+X[,4]) * z + rnorm(n, mean = 0, sd = 1))
        C = 0.3+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
      }else if (censor == "50%"){

        ft <- exp(0.3+c(Xc %*% c(0.3,0.2,0.1)) + c(Xb %*% c(-0.3,-0.3,-0.2))+
                    (-0.5+X[,4]) * z + rnorm(n, mean = 0, sd = sigma))
        C = 0.18+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
      }
      
    } else if (b.b == "N" & b.i == "Y") {

      ft <- exp(-0.4 + 0.2*X[,4]+(-1 + c(Xc %*% c(0.3,0.2,0.1)) 
                                  + c(Xb %*% c(-0.3,-0.3,-0.2))) * z + 
                  rnorm(n, mean = 0, sd = sigma))
      C = 0.21+
        (-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
    } else if (b.b == "Y" & b.i == "Y") {

      ft <- exp(0.1 + 0.2*X[,4]+c(Xb %*% c(-0.3,-0.3,-0.2))+
                  (-1 + c(Xc %*% c(0.3,0.2,0.1))) * z +  
                  rnorm(n, mean = 0, sd = sigma))
      C = 0.34+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
    }
  } else if (f.b == "L" & f.i == "N") {
    ft <- exp(-0.3 + c(Xc %*% c(0.3,0.2,0.1)) + c(Xb %*% c(-0.3,-0.3,-0.2))+
                (-0.3 + sin(X[,4])) * z + rnorm(n, mean = 0, sd = sigma))
    C = 0.35+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
  } else if (f.b == "N" & f.i == "N") {
    ft <- exp(0.2 + c(sin(Xc) %*% c(0.3,0.2,0.1)) + c(cos(Xc) %*% c(-0.3,-0.3,-0.2))+
                (-1 + sin(X[,4])) * z + rnorm(n, mean = 0, sd = sigma))
    C = 0.53+(-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
  }
  
  # generate Time
  Time = ft

  Event = as.numeric(Time<C)
  Time = ifelse(Time<C,Time,C)
  Time <- pmax(rep(0.001, length(Time)), Time)
  
  # cate true
  cate <- rep(NA, n)
  catesp <- rep(NA, n)
  eps <- rnorm(n.mc)
  for (i in 1:n) {
    if (f.b == "L" & f.i == "L") {
      if (b.b == "Y" & b.i == "N") {
        if (censor == "25%" & dgps == "nft"){
          
          ft0 <- exp(-0.3 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                       exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
          ft1 <- exp(-0.3 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                      (-0.5+X[i,4]) + exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
          
        }else if (censor == "25%" & dgps == "aft"){
          
          ft0 <- exp(-0.6 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                       eps)
          ft1 <- exp(-0.6 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                       (-0.5+X[i,4]) + eps)
         
        }else if (censor == "50%"){
          
          ft0 <- exp(c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                       exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
          ft1 <- exp(c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                       (-0.5+X[i,4]) + exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
          
        }
        
      } else if (b.b == "N" & b.i == "Y") {
        
        ft0 <- exp(-0.4 + 0.2*X[i,4] + 
                    exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
        ft1 <- exp(-0.4 + 0.2*X[i,4]+(-1 + c(Xc[i,] %*% c(0.3,0.2,0.1)) 
                                    + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))) + 
                    exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
        
      } else if (b.b == "Y" & b.i == "Y") {
        
        ft0 <- exp(0.1 + 0.2*X[i,4]+c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                    exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
        ft1 <- exp(0.1 + 0.2*X[i,4]+c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                    (-1 + c(Xc[i,] %*% c(0.3,0.2,0.1))) +  
                    exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
        
      }
    } else if (f.b == "L" & f.i == "N") {
      
      ft0 <- exp(-0.3 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                 exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
      ft1 <- exp(-0.3 + c(Xc[i,] %*% c(0.3,0.2,0.1)) + c(Xb[i,] %*% c(-0.3,-0.3,-0.2))+
                  (-0.3 + sin(X[i,4])) + exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
      
      
    } else if (f.b == "N" & f.i == "N") {
      
      ft0 <- exp(0.2 + c(sin(Xc[i,]) %*% c(0.3,0.2,0.1)) + c(cos(Xc[i,]) %*% c(-0.3,-0.3,-0.2))+
                 exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
      ft1 <- exp(0.2 + c(sin(Xc[i,]) %*% c(0.3,0.2,0.1)) + c(cos(Xc[i,]) %*% c(-0.3,-0.3,-0.2))+
                  (-1 + sin(X[i,4])) + exp(-2+0.3*X[i,12]+0.2*X[i,5])*eps)
      
      
    }
    
    cate[i] <- mean(ft1 - ft0)
    catesp[i] <- mean(ft1 > t0) - mean(ft0 > t0)
    
  }

  catesp.sign <- sign(catesp)
  X = as.matrix(X)
  Y = Time
  W = z
  D = Event

  # w.hat
  SL.library = c("SL.glm", "SL.ranger", "SL.glmnet", "SL.xgboost", "SL.gam")
  X = as.data.frame(X)
  W.hat <- predict(SuperLearner(Y = W, X = X, SL.library = SL.library,
                                verbose = TRUE, method = "method.NNLS",
                                family = binomial()))$pred

  if (train == TRUE) {

    # fit model on W == 1
    grffit1 <- nftbart::nft2(X[W == 1,, drop = FALSE],
                             X[W == 1,, drop = FALSE],
                             Y[W == 1],
                             D[W == 1], K=0)
    surf1 <- predict(grffit1,X,X)

    # fit model on W == 0
    grffit0 <- nftbart::nft2(X[W == 0,, drop = FALSE],
                             X[W == 0,, drop = FALSE],
                             Y[W == 0],
                             D[W == 0], K=0)
    surf0 <- predict(grffit0,X,X)

    S1.hat <- apply(exp(surf1[["f.test"]]), 2, function(x) mean(x > t0))
    S0.hat <- apply(exp(surf0[["f.test"]]), 2, function(x) mean(x > t0))


    # cen.fit
    cen.fit <- nftbart::nft2(cbind(W, X),
                             cbind(W, X),
                             Y,
                             1 - D, K=0)
    C.hat <- predict(cen.fit,cbind(W, X),cbind(W, X))
    Chat <- apply(exp(C.hat[["f.test"]]), 2, function(x) mean(x > t0))
    
  }else{
    S1.hat = NULL
    S0.hat = NULL
    grffit1 = NULL
    grffit0 = NULL
    Chat = NULL
  }


  return(
    list(
      X = X,
      Y = Y,
      W = W,
      D = D,
      catesp = catesp,
      catesp.sign = catesp.sign,
      W.hat = W.hat,
      surf1 = S1.hat,
      surf0 = S0.hat,
      fit1 = grffit1,
      fit0 = grffit0,
      C.hat = Chat
    )
  )

}

