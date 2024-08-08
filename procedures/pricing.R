proba.def.Q.quarters.fction.of.surplus <- function(values.fs,
                                                   alpha,
                                                   nb.Q){
  # values.fs are the values of the (log) fiscal space
  E.exp.mAlpha.max <- exp(-alpha*pmax(values.fs,0))
  P <- 1 - E.exp.mAlpha.max^nb.Q
  return(P)
}

PDef <- function(model.sol,X,max.h,indic.P=1){
  if(indic.P==1){
    PD <- compute.physical.PD(model.sol,X,max.h)
  }else{# Under Q
    aux <- compute.CDS(model.sol,X,
                       max.h = max.h,
                       indic.cpp = 1)
    PD <- 1 - aux$E.n/aux$P.rf
  }
  return(PD)
}

compute.debt.limits <- function(model.sol,X0,Proba,h,indic.P=1,
                                MAX.d=20,MIN.d=0,
                                max.tol=.001){
  # This function computes the debt value that would be consistent
  # with an horizon-h proba of default equal to Proba.
  n_x <- dim(X0)[2]
  
  d <- X0[,n_x-1]
  
  min.d <- d - 1
  X <- X0
  X[,c(n_x-1,n_x)] <- min.d
  min.PD <- PDef(model.sol,X,h,indic.P)[,h]
  if(sum(is.na(min.PD))==length(d)){
    return(NaN*d)
  }else{
    while((sum(min.PD>Proba)>0)&(sum(min.d<MIN.d)!=length(d))){
      min.d <- min.d - 1
      X <- X0
      X[,c(n_x-1,n_x)] <- min.d
      min.PD <- PDef(model.sol,X,h,indic.P)[,h]
    }
    
    max.d <- d + 8
    X <- X0
    X[,c(n_x-1,n_x)] <- max.d
    max.PD <- PDef(model.sol,X,h,indic.P)[,h]
    while((sum(max.PD<Proba)>0)&(sum(min.d>MAX.d)!=length(d))){
      max.d <- max.d + 1
      X <- X0
      X[,c(n_x-1,n_x)] <- max.d
      max.PD <- PDef(model.sol,X,h,indic.P)[,h]
    }
    
    for(i in 1:20){
      mid.d <- .5*(min.d + max.d)
      X <- X0
      X[,c(n_x-1,n_x)] <- mid.d
      middle.PD <- PDef(model.sol,X,h,indic.P)[,h]
      
      min.d <- mid.d * (middle.PD < Proba) + min.d * (middle.PD >= Proba)
      max.d <- max.d * (middle.PD < Proba) + mid.d * (middle.PD >= Proba)
    }
    d <- .5*(min.d + max.d)
    X <- X0
    X[,c(n_x-1,n_x)] <- d
    PD <- PDef(model.sol,X,h,indic.P)[,h]
    
    erreur <- abs(PD-Proba)
    print(max(abs(PD-Proba)))
    
    d[(erreur>max.tol)&(d>MIN.d)&(d<MAX.d)] <- NaN
    d[d<MIN.d] <- MIN.d
    d[d>MAX.d] <- MAX.d
  }
  return(d)
}

compute.debt.limits_OLD <- function(model.sol,X0,Proba,h,indic.P=1,max.tol=.001){
  # This function computes the debt value that would be consistent
  # with an horizon-h proba of default equal to Proba.
  n_x <- dim(X0)[2]
  
  epsilon <- .01
  
  d <- X0[,n_x-1] + 5
  
  for(i in 1:200){
    # Compute f(d):
    X <- X0
    X[,c(n_x-1,n_x)] <- d
    PD <- PDef(model.sol,X,h,indic.P)
    PD <- PD[,h]
    f <- PD - Proba
    
    # Compute f(d_perturb)
    d_perturb <- X[,n_x-1] + epsilon
    X_perturb <- X0
    X_perturb[,c(n_x-1,n_x)] <- d_perturb
    PD_perturb <- PDef(model.sol,X_perturb,h,indic.P)
    PD_perturb <- PD_perturb[,h]
    f_perturb <- PD_perturb - Proba
    df <- (f_perturb - f)/epsilon
    factor <- min(max(1/(10*abs(f/df)),.02),.2)
    #factor <- min(max(1/(40*abs(f/df)),.01),.1)
    #print(cbind(max(factor),max(abs(f)),max(abs(df))))
    d <- pmax(d - factor*f/df,0)
  }
  print(max(abs(f)))
  d[abs(f)>max.tol] <- 0
  return(d)
}

A <- function(model,u.matrix){
  # u.matrix is of dimension n * k, where k is
  #     the number of vectors at which we want to evaluate A
  A <- t(model$Phi) %*% u.matrix
  return(A)
}

B <- function(model,u.matrix){
  # u.matrix is of dimension n * k, where k is
  #     the number of vectors at which we want to evaluate B
  aux <-  t(u.matrix)
  B <- 1/2 * apply(aux * t(u.matrix),1,sum)
  return(B)
}


make.Gamma <- function(Psi){
  # Psi is of dimension (H+1)x(H+1) (see paper)
  # The sum of the entries (j+1,1),(j+2,2),...,(n,n-j) of Psi is Gamma_{n,j}
  # This function computes the matrix Mat.Gamma whose component (i,j) is Gamma_{j,i-1}
  # That is, Gamma_{i+j-1,i-1} is the component (i,j) of Mat.Gamma,
  #   or Gamma_{n,k} is the component  (k+1,n-k) of Mat.Gamma
  
  n <- dim(Psi)[1]
  indic.entries.Psi.to.keep <- which(!upper.tri(Psi))
  
  AUX <- (!lower.tri(Psi))[,n:1]
  indic.entries.AUX.to.fill <- which(AUX)
  
  AUX <- matrix(0,n,n) #NaN
  AUX[indic.entries.AUX.to.fill] <- Psi[indic.entries.Psi.to.keep]
  
  Mat.Gamma <- t(apply(AUX,1,cumsum))
  return(Mat.Gamma)
}

make.indices <- function(M){
  n.m <- dim(M)[1]
  mat.indices <- matrix(NaN, (n.m-1), (n.m-1))
  for (i in 1:(n.m-1)){
    mat.indices[i:(n.m-1),i] <- 2:(n.m-i+1) + (i-1)*n.m
  }
  return(mat.indices)
}



compute.AB.nom.yields <- function(model.sol){
  
  Lambda <- model.sol$gamma * model.sol$Lambda_c + model.sol$Lambda_pi
  Tau     <- log(model.sol$delta) -
    model.sol$mu_c - model.sol$mu_pi +
    .5 * sum((model.sol$gamma * model.sol$sigma_c + model.sol$sigma_pi)^2)
  
  n_w       <- ncol(model.sol$Phi)
  
  all.b.h <- NULL
  all.a.h <- NULL
  all.B.h <- NULL
  all.A.h <- NULL
  a.h <- matrix(0,n_w,1)
  b.h <- 0
  
  for (h in 1:model.sol$max.h){
    # This loop also computes the (a,b) loadings to compute E(exp( - Lambda (w_{t+1}+...+w_{t+h})))
    a.h <- model.sol$Phi %*% ( a.h - Lambda )
    b.h <- b.h + 1/2 * sum(( a.h - Lambda )^2)
    all.a.h <- cbind(all.a.h,- 1/h * a.h)
    all.b.h <- cbind(all.b.h, - 1/h * b.h - Tau)
    all.A.h <- cbind(all.A.h, a.h)
    all.B.h <- cbind(all.B.h, b.h + h * Tau)
  }
  
  return(list(
    all.a.h = all.a.h,
    all.b.h = all.b.h,
    all.A.h = all.A.h,
    all.B.h = all.B.h))
}


compute.max.Sharpe.Ratio <- function(model.sol,X){
  # this function computes the maximum Sharpe Ratio for a given X
  
  T <- dim(X)[1]
  n_x   <- dim(X)[2]
  n.eta <- length(model.sol$sigma_c)
  n_w   <- nrow(model.sol$Phi)
  
  Mu_x    <- model.sol$Mu_x
  Phi_x   <- model.sol$Phi_x
  Sigma_x <- model.sol$Sigma_x
  
  Lambda_l <- model.sol$a.l
  mu_l     <- c(model.sol$b.l)
  
  Lambda_c <- model.sol$Lambda_c
  sigma_c  <- model.sol$sigma_c
  
  b_c   <- model.sol$b_c
  alpha <- model.sol$alpha
  gamma <- model.sol$gamma
  
  Lambda.l.nx <- matrix(c(Lambda_l,rep(0,n_x-n_w)),ncol=1)
  
  # Extend the state vector with w_{t-1} (dimension n_w), d_{t-1} (dimension 1) and xi_t (dimension 1):
  Mu_x  <- matrix(c(Mu_x,rep(0,n_w + 2)),ncol=1)
  Phi_x <- matrix(0,n_x+n_w+2,n_x+n_w+2)
  Phi_x[1:n_x,1:n_x] <- model.sol$Phi_x
  Phi_x[(n_x+1):(n_x+n_w),1:n_w] <- diag(n_w) # w_{t-1}
  Phi_x[n_x+n_w+1,n_x - 2] <- 1 # d_{t-1}
  Sigma_x <- rbind(cbind(Sigma_x,0),
                   matrix(0,n_w+2,dim(Sigma_x)[2]+1))
  
  Sigma_x[n_x+n_w+2,dim(Sigma_x)[2]] <- model.sol$sigma.nu
  
  Lambda.l.nx <- matrix(c(rep(0,n_x),Lambda_l,0,1),ncol=1)
  Lambda.d.nx <- matrix(0, n_x+n_w+2,1)
  Lambda.d.nx[n_x+n_w+1] <- 1
  a <- alpha*(Lambda.d.nx - Lambda.l.nx)
  b <- - alpha*mu_l
  
  X <- cbind(X,
             matrix(0,dim(X)[1],n_w+2))
  if(T > 1){
    X[1,(n_x+1):(n_x+n_w+1)]   <- X[1,c(1:n_w,n_x-2)]
    X[2:T,(n_x+1):(n_x+n_w+1)] <- X[1:(T-1),c(1:n_w,n_x-2)]
  }else{
    X[1,(n_x+1):(n_x+n_w+1)]   <- X[1,c(1:n_w,n_x-2)]
  }
  
  # Prepare u and v:
  u <- matrix(0,n_x + n_w + 2,1)
  u[1:3] <- - gamma * Lambda_c
  u[5:7] <- - gamma * sigma_c
  v <- gamma * b_c
  
  res.uv <- compute.cond.E.Var(c(u),v,
                               a,b,
                               Phi_x,Mu_x,Sigma_x,
                               X)
  res.2uv <- compute.cond.E.Var(c(2 * u),2 * v,
                                a,b,
                                Phi_x,Mu_x,Sigma_x,
                                X)
  
  E.exp.2m <- exp(res.2uv$E.ux.vD + .5 * res.2uv$Var.ux.vD)
  E.exp.m  <- exp(res.uv$E.ux.vD  + .5 * res.uv$Var.ux.vD)
  
  Var.M <- E.exp.2m - E.exp.m^2
  E.M   <- E.exp.m
  
  max.SR <- sqrt(Var.M) / E.M
  
  return(max.SR)
}

compute.cond.E.Var <- function(u,v,
                               a,b,
                               Phi_x,Mu_x,Sigma_x,
                               X){
  # compute E_t(u'x_{t+1} + v D_{t+1})
  # compute Var_t(u'x_{t+1} + v D_{t+1})
  
  T <- dim(X)[1]
  vec.1 <- matrix(1,T,1)
  
  # Conditional expectation of x_{t+1} (E_t(x_{t+1})):
  E.xp1 <- vec.1 %*% matrix(Mu_x,nrow=1) + X %*% t(Phi_x)
  
  mu.lambda.t <- c(b + E.xp1 %*% matrix(a,ncol=1))
  sigma.lambda <- c(sqrt(matrix(a,nrow=1) %*% Sigma_x %*% t(Sigma_x) %*% matrix(a,ncol=1)))
  
  E.ux.vD <- E.xp1 %*% matrix(u,ncol=1) +
    v * pnorm( mu.lambda.t/sigma.lambda) * mu.lambda.t +
    v * dnorm(-mu.lambda.t/sigma.lambda) * sigma.lambda
  
  Var.lambda <- pnorm( mu.lambda.t/sigma.lambda) * sigma.lambda^2
  Cov.ux.lambda <- pnorm( mu.lambda.t/sigma.lambda) *
    c(matrix(u,nrow=1) %*% Sigma_x %*% t(Sigma_x) %*% matrix(a,ncol=1))
  E.lambda <- 
    pnorm( mu.lambda.t/sigma.lambda) * mu.lambda.t +
    dnorm(-mu.lambda.t/sigma.lambda) * sigma.lambda
  
  # Law of total variance:
  Var.ux.vD <-
    c(matrix(u,nrow=1) %*% Sigma_x %*% t(Sigma_x) %*% matrix(u,ncol=1)) + 
    v^2 * Var.lambda +
    2 * v * Cov.ux.lambda +
    v^2 * E.lambda
  
  return(list(
    E.ux.vD   = E.ux.vD,
    Var.ux.vD = Var.ux.vD
  ))
}


compute.CDS <- function(model.sol,X,
                        max.h,
                        indic.cpp = 1){
  
  T <- dim(X)[1]
  n_x   <- dim(X)[2]
  n_w   <- nrow(model.sol$Phi)
  
  Mu_x    <- model.sol$Mu_x
  Phi_x   <- model.sol$Phi_x
  Sigma_x <- model.sol$Sigma_x
  
  sigma  <- model.sol$gamma * model.sol$sigma_c + model.sol$sigma_pi
  Tau     <- log(model.sol$delta) - model.sol$mu_c - model.sol$mu_pi
  b_c  <- model.sol$b_c
  b_pi <- model.sol$b_pi
  
  a.dot <- matrix(c(sigma,rep(0,n_x-n_w)), ncol=1) 
  b.dot <- 0
  
  alpha  <- model.sol$alpha
  beta   <- model.sol$beta
  delta  <- model.sol$delta
  gamma  <- model.sol$gamma
  d_star <- model.sol$d_star
  
  a           <- matrix(0,n_x,1)
  a[1:n_w]    <- model.sol$alpha * model.sol$sigma_s -
    model.sol$alpha * model.sol$sigma_star
  # a[1:n_w]    <-  - model.sol$alpha * model.sol$sigma_star
  a[n_w + 3]  <- model.sol$alpha * model.sol$beta
  b           <- model.sol$alpha * (- model.sol$beta * model.sol$d_star - model.sol$mu_star)
  
  # Compute conditional expectations:
  if(indic.cpp == 0){
    res.cond.Expect <- compute.condit.Exp(a,b,
                                          a.dot,b.dot,
                                          Phi_x,Mu_x,Sigma_x,
                                          X,max.h)
  }else{
    a     <- matrix(a,ncol=1)
    a.dot <- matrix(a.dot,ncol=1)
    res.cond.Expect <- compute_condit_Exp_cpp(a,b,
                                              a.dot,b.dot,
                                              Phi_x,Mu_x,Sigma_x,
                                              X,max.h)
    res.cond.Expect.bis <- compute_condit_Exp_cpp(0.00001*a,0.00001*b,
                                                  a.dot,b.dot,
                                                  Phi_x,Mu_x,Sigma_x,
                                                  X,max.h)
  }
  
  E.n      <- res.cond.Expect$E.n
  E.n.star <- res.cond.Expect$E.n.star
  
  matrix.Tau          <- t(matrix(exp(Tau)^(1:max.h),max.h,dim(X)[1]))
  
  E.n      <- matrix.Tau * E.n
  
  # all.a.h  <- res.cond.Expect$all.a.h
  # all.b.h  <- res.cond.Expect$all.b.h
  # E.n.rf   <- matrix.Tau * t(matrix(exp(all.b.h),max.h,dim(X)[1])) *
  #  exp(X %*% all.a.h)
  
  E.n.bis  <- res.cond.Expect.bis$E.n # E.n but with default intensity = 0
  E.n.rf   <- matrix.Tau * E.n.bis
  
  P.n      <- E.n # E(Mx(1-D))
  P.rf <- exp(gamma * b_c + b_pi) * E.n.rf +
    (1 - exp(gamma * b_c + b_pi)) * E.n  # E(M)
  P.n.star <- P.rf - P.n  # E(MxD)
  yds.rf  = - log(P.rf) * t(matrix(1/(1:max.h),max.h,dim(X)[1]))
  yds.gov = - log(E.n.rf) * t(matrix(1/(1:max.h),max.h,dim(X)[1]))
  
  cumsum.P.n      <- t(apply(P.n,1,cumsum))
  cumsum.P.n.star <- t(apply(P.n.star,1,cumsum))
  
  CDS.spds <- (1-model.sol$RR) * P.n.star/cumsum.P.n
  
  return(list(
    CDS.spds = CDS.spds,
    cumsum.P.n.star = cumsum.P.n.star,
    P.n = P.n,
    P.n.star = P.n.star,
    E.n = E.n,
    E.n.star = E.n.star,
    E.n.rf = E.n.rf,
    cumsum.P.n = cumsum.P.n,
    P.rf = P.rf,
    yds.rf = yds.rf,
    yds.gov = yds.gov
  ))
}







compute.LT_ESSAI <- function(model.sol,
                       X,      # state variables
                       u,      # argument
                       max.h   # maximum horizon
){
  T <- dim(X)[1]
  n_x   <- dim(X)[2]
  n_w   <- nrow(model.sol$Phi)
  
  Mu_x    <- model.sol$Mu_x
  Phi_x   <- model.sol$Phi_x
  Sigma_x <- model.sol$Sigma_x
  
  b_c    <- model.sol$b_c
  b_pi   <- model.sol$b_pi
  
  sigma  <- model.sol$gamma * model.sol$sigma_c + model.sol$sigma_pi
  a.dot  <- matrix(c(sigma,rep(0,n_x-n_w)), ncol=1) 
  b.dot  <- 0
  
  alpha  <- model.sol$alpha
  beta   <- model.sol$beta
  delta  <- model.sol$delta
  gamma  <- model.sol$gamma
  d_star <- model.sol$d_star
  
  a           <- matrix(0,n_x,1)
  a[1:n_w]    <- - model.sol$alpha * model.sol$sigma_star
  a[n_w + 3]  <- model.sol$alpha * model.sol$beta
  b           <- model.sol$alpha * (- model.sol$beta * model.sol$d_star - model.sol$mu_star)
  
  # Compute conditional expectations:
  a     <- matrix(a,ncol=1)
  a.dot <- matrix(a.dot,ncol=1)
  
  res.cond.Expect.Q.no.def <- compute.condit.Exp(0.00001*a,0.00001*b,
                                                 a.dot,b.dot,
                                                 Phi_x,Mu_x,Sigma_x,
                                                 X,max.h,u)
  res.cond.Expect.Q <- compute.condit.Exp(a,b,
                                          a.dot,b.dot,
                                          Phi_x,Mu_x,Sigma_x,
                                          X,max.h,u)
  res.cond.Expect.P <- compute.condit.Exp(0.00001*a,0.00001*b,
                                          0*a.dot,0*b.dot,
                                          Phi_x,Mu_x,Sigma_x,
                                          X,max.h,u)
  
  E.n.Q.no.u          <- res.cond.Expect.Q$E.n # E(Mx(1-D))
  E.n.Q.with.u        <- res.cond.Expect.Q$E.n.with.u # E(Mx(1-D)exp(u))
  E.n.Q.no.def.no.u   <- res.cond.Expect.Q.no.def$E.n # E(exp(phi1))
  E.n.Q.no.def.with.u <- res.cond.Expect.Q.no.def$E.n.with.u # E(exp(phi1+u))
  
  LT.P                <- res.cond.Expect.P$E.n.with.u # E(exp(u))
  
  ExpectM <- exp(gamma * b_c + b_pi) * E.n.Q.no.def.no.u +
    (1 - exp(gamma * b_c + b_pi)) * E.n.Q.no.u  # E(M)
  LT.Q <- exp(gamma * b_c + b_pi) * E.n.Q.no.def.with.u +
    (1 - exp(gamma * b_c + b_pi)) * E.n.Q.with.u  # E(M)
  LT.Q <- LT.Q/ExpectM
  
  return(list(
    LT.P = LT.P,
    LT.Q = LT.Q
  ))
}





compute.LT <- function(model.sol,
                           X,      # state variables
                           u,      # argument
                           max.h,   # maximum horizon
                           indic.cpp = 1
){
  T <- dim(X)[1]
  n_x   <- dim(X)[2]
  n_w   <- nrow(model.sol$Phi)
  
  Mu_x    <- model.sol$Mu_x
  Phi_x   <- model.sol$Phi_x
  Sigma_x <- model.sol$Sigma_x
  
  b_c    <- model.sol$b_c
  b_pi   <- model.sol$b_pi
  
  sigma  <- model.sol$gamma * model.sol$sigma_c + model.sol$sigma_pi
  a.dot  <- matrix(c(sigma,rep(0,n_x-n_w)), ncol=1) 
  b.dot  <- 0
  
  alpha  <- model.sol$alpha
  beta   <- model.sol$beta
  delta  <- model.sol$delta
  gamma  <- model.sol$gamma
  d_star <- model.sol$d_star
  
  a           <- matrix(0,n_x,1)
  a[1:n_w]    <- - model.sol$alpha * model.sol$sigma_star
  a[n_w + 3]  <- model.sol$alpha * model.sol$beta
  b           <- model.sol$alpha * (- model.sol$beta * model.sol$d_star - model.sol$mu_star)
  
  # Compute conditional expectations:
  a     <- matrix(a,ncol=1)
  a.dot <- matrix(a.dot,ncol=1)
  if(indic.cpp == 0){
    res.cond.Expect.Q.no.u <- compute.condit.Exp(a,b,
                                                 a.dot,b.dot,
                                                 Phi_x,Mu_x,Sigma_x,
                                                 X,max.h)
    res.cond.Expect.Q.no.def.no.u  <- compute.condit.Exp(0.00001*a,0.00001*b,
                                                         a.dot,b.dot,
                                                         Phi_x,Mu_x,Sigma_x,
                                                         X,max.h)
    res.cond.Expect.Q.no.def.with.u <- compute.condit.Exp(0.00001*a,0.00001*b,
                                                          a.dot+u,b.dot,
                                                          Phi_x,Mu_x,Sigma_x,
                                                          X,max.h)
    res.cond.Expect.Q.with.u <- compute.condit.Exp(a,b,
                                                   a.dot+u,b.dot,
                                                   Phi_x,Mu_x,Sigma_x,
                                                   X,max.h)
    res.cond.Expect.P <- compute.condit.Exp(0.00001*a,0.00001*b,
                                            u,b.dot,
                                            Phi_x,Mu_x,Sigma_x,
                                            X,max.h)
  }else{
    res.cond.Expect.Q.no.u <- compute_condit_Exp_cpp(a,b,
                                                     a.dot,b.dot,
                                                     Phi_x,Mu_x,Sigma_x,
                                                     X,max.h)
    res.cond.Expect.Q.no.def.no.u  <- compute_condit_Exp_cpp(0.00001*a,0.00001*b,
                                                             a.dot,b.dot,
                                                             Phi_x,Mu_x,Sigma_x,
                                                             X,max.h)
    res.cond.Expect.Q.no.def.with.u <- compute_condit_Exp_cpp(0.00001*a,0.00001*b,
                                                              a.dot+u,b.dot,
                                                              Phi_x,Mu_x,Sigma_x,
                                                              X,max.h)
    res.cond.Expect.Q.with.u <- compute_condit_Exp_cpp(a,b,
                                                       a.dot+u,b.dot,
                                                       Phi_x,Mu_x,Sigma_x,
                                                       X,max.h)
    res.cond.Expect.P <- compute_condit_Exp_cpp(0.00001*a,0.00001*b,
                                                u,b.dot,
                                                Phi_x,Mu_x,Sigma_x,
                                                X,max.h)
  }
  
  E.n.Q.no.u          <- res.cond.Expect.Q.no.u$E.n # E(Mx(1-D))
  E.n.Q.with.u        <- res.cond.Expect.Q.with.u$E.n # E(Mx(1-D)exp(u))
  E.n.Q.no.def.no.u   <- res.cond.Expect.Q.no.def.no.u$E.n # E(exp(phi1))
  E.n.Q.no.def.with.u <- res.cond.Expect.Q.no.def.with.u$E.n # E(exp(phi1+u))
  
  LT.P                <- res.cond.Expect.P$E.n # E(exp(u))
  
  ExpectM <- exp(gamma * b_c + b_pi) * E.n.Q.no.def.no.u +
    (1 - exp(gamma * b_c + b_pi)) * E.n.Q.no.u  # E(M)
  LT.Q <- exp(gamma * b_c + b_pi) * E.n.Q.no.def.with.u +
    (1 - exp(gamma * b_c + b_pi)) * E.n.Q.with.u  # E(M)
  LT.Q <- LT.Q/ExpectM
  
  return(list(
    LT.P = LT.P,
    LT.Q = LT.Q
  ))
}


compute.physical.PD <- function(model.sol,X,
                                max.h){
  
  T <- dim(X)[1]
  n_x   <- dim(X)[2]
  n_w   <- nrow(model.sol$Phi)
  
  Mu_x    <- model.sol$Mu_x
  Phi_x   <- model.sol$Phi_x
  Sigma_x <- model.sol$Sigma_x
  
  a.dot <- matrix(0,n_x,1)
  b.dot <- 0
  
  alpha  <- model.sol$alpha
  beta   <- model.sol$beta
  delta  <- model.sol$delta
  gamma  <- model.sol$gamma
  d_star <- model.sol$d_star
  
  a           <- matrix(0,n_x,1)
  a[1:n_w]    <- model.sol$alpha * model.sol$sigma_s - model.sol$alpha * model.sol$sigma_star
  a[n_w + 3]  <- model.sol$alpha * model.sol$beta
  b           <- model.sol$alpha * (- model.sol$beta * model.sol$d_star - model.sol$mu_star)
  
  # Compute conditional expectations:
  res.cond.Expect <- compute_condit_Exp_cpp(a,b,
                                            a.dot,b.dot,
                                            Phi_x,Mu_x,Sigma_x,
                                            X,max.h)
  
  E.n      <- res.cond.Expect$E.n
  
  return(1-E.n)
}










compute.condit.Exp <- function(a,b,
                               a.dot,b.dot,
                               Phi_x,Mu_x,Sigma_x,
                               X,max.h,
                               u = 0){
  # This function computes two types of conditional expectations,
  # (Those appearing in Eqs. I.5 and I.6)
  
  n_x <- dim(Phi_x)[1]
  T <- dim(X)[1]
  
  phi1 <- - matrix(a.dot, ncol=1)
  
  Phi.x.h   <- diag(n_x)
  MU        <- array(NaN, c(n_x, T, max.h))
  vec.1.T   <- matrix(1, T, 1)
  
  Omega <- Sigma_x%*%t(Sigma_x)
  
  GAMMA <- array(NaN, c(n_x,n_x, max.h)) # The i^th layer is Gamma_{i,0}. In particular, the first layer if Gamma_{1,0} = Omega
  
  gamma.h.0 <- matrix(0, n_x, n_x)
  
  KSI.a        <- matrix(NaN, n_x, max.h + 1) # the i^th column of KSI.a will be xi_{i-1}^a = t(Phi_x^i) x a
  KSI.a[,1]    <- a
  KSI.a.dot    <- matrix(NaN, n_x, max.h + 1)  # the i^th column of KSI.a.dot will be xi_{i-1}^a.dot = t(Phi_x^i) x a.dot
  KSI.a.dot[,1]<- a.dot
  
  SIGMA.lambda <- matrix(NaN,max.h,1)
  
  MU.lambda <- matrix(NaN, T, max.h)
  MU.i        <-  matrix(NaN, T, max.h)
  MU.i.with.u <-  matrix(NaN, T, max.h)
  
  all.b.h <- NULL
  all.a.h <- NULL
  a.h <- matrix(0,n_x,1)
  b.h <- 0
  sum.Phi.x.h_1 <- 0
  
  for (h in 1:max.h){
    # This loop also computes the (a,b) loadings to compute E(exp( - Lambda (w_{t+1}+...+w_{t+h})))
    # It also computes matrices MU.lambda and MU.i
    
    b.h <- b.h + 1/2 * c(t( a.h + phi1 ) %*% Omega %*% ( a.h + phi1 )) +
      c(t( a.h + phi1 ) %*% Mu_x)
    a.h <- t(Phi_x) %*% ( a.h + phi1 )
    all.a.h <- cbind(all.a.h,a.h)
    all.b.h <- cbind(all.b.h,b.h)
    
    sum.Phi.x.h_1 <- sum.Phi.x.h_1 + Phi.x.h
    Phi.x.h <- Phi.x.h%*%Phi_x # Phi.x.h is Phi_x^h
    
    MU[,,h] <- (sum.Phi.x.h_1%*%Mu_x)%*%t(vec.1.T) +
      Phi.x.h%*%t(X)
    
    gamma.h.0 <- Omega + Phi_x%*%gamma.h.0%*%t(Phi_x)
    
    GAMMA[,,h] <- gamma.h.0
    
    SIGMA.lambda[h] <- sqrt(t(a)%*%gamma.h.0%*%a)
    
    MU.lambda[, h] <- b +  t(MU[,,h])%*%a
    MU.i[,h]        <- b.dot + t(MU[,,h])%*%a.dot 
    MU.i.with.u[,h] <- b.dot + t(MU[,,h])%*%(a.dot + u)
    
    KSI.a[,h+1]    <- Phi.x.h%*%a
    KSI.a.dot[,h+1]<- Phi.x.h%*%a.dot
  }
  
  SIGMA.lambda <- vec.1.T%*%t(SIGMA.lambda) # This is matrix of dimension T x (max.h - 1)
  
  AUX.lambda <-  MU.lambda/SIGMA.lambda
  
  P <- pnorm(AUX.lambda)
  PSI.a         <- t(KSI.a)%*%Omega%*%KSI.a
  PSI.a.dot     <- t(KSI.a.dot)%*%Omega%*%KSI.a.dot
  PSI.a.dot.a     <- t(KSI.a.dot)%*%Omega%*%KSI.a
  PSI.a.a.dot     <- t(KSI.a)%*%Omega%*%KSI.a.dot
  PSI.dot.dot     <- t(KSI.a + KSI.a.dot)%*%Omega%*%(KSI.a + KSI.a.dot)
  mat.GAMMA           <- make.Gamma(PSI.a)
  mat.GAMMA.dot       <- make.Gamma(PSI.a.dot)
  mat.GAMMA.left.dot  <- make.Gamma(PSI.a.dot.a)
  mat.GAMMA.right.dot <- make.Gamma(PSI.a.a.dot)
  mat.GAMMA.dot.dot   <- make.Gamma(PSI.dot.dot)
  
  MU.0 <- t(X)
  MU.i.lagged <- cbind(b.dot + t(MU.0)%*%a.dot,MU.i[,1:(max.h-1)])
  
  MU.i.lagged        <- MU.i
  MU.i.lagged.with.u <- MU.i.with.u
  
  F.n_1.n      <- matrix(NaN,T,max.h) # The first column is for (n-1,n) = (0,1)
  F.n_1.n.star <- matrix(NaN,T,max.h) # The first column is for (n-1,n) = (0,1)
  
  
  # Treat F.n_1.n: ------------------ 
  F.n_1.n        <- MU.i.lagged        + P*MU.lambda + dnorm(-AUX.lambda)*SIGMA.lambda
  F.n_1.n.with.u <- MU.i.lagged.with.u + P*MU.lambda + dnorm(-AUX.lambda)*SIGMA.lambda
  
  F.n_1.n[,2:max.h] <- F.n_1.n[,2:max.h] -
    1/2*(P[,2:(dim(P)[2])]*(vec.1.T%*%matrix(mat.GAMMA.dot.dot[1,2:(dim(P)[2])],nrow=1))+
           (1-P[,2:(dim(P)[2])])*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1))) 
  F.n_1.n.with.u[,2:max.h] <- F.n_1.n.with.u[,2:max.h] -
    1/2*(P[,2:(dim(P)[2])]*(vec.1.T%*%matrix(mat.GAMMA.dot.dot[1,2:(dim(P)[2])],nrow=1))+
           (1-P[,2:(dim(P)[2])])*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1))) 
  
  
  # Treat F.n_1.n.star: ------------------
  F.n_1.n.star <- MU.i.lagged
  
  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] + P[,1:(dim(P)[2]-1)]*MU.lambda[,1:(dim(P)[2]-1)] +
    dnorm(-AUX.lambda[,1:(dim(P)[2]-1)])*SIGMA.lambda[,1:(dim(P)[2]-1)]
  
  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    1/2*(P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1))+
           (1-P[,2:dim(P)[2]])*(vec.1.T%*%matrix(mat.GAMMA.dot[1,2:(dim(P)[2])],nrow=1))) 
  
  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA.right.dot[2,2:(dim(P)[2])],nrow=1))
  
  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] - 
    1/2*(P[,2:dim(P)[2]]*(vec.1.T%*%matrix(mat.GAMMA[1,2:dim(P)[2]],nrow=1)))
  
  F.n_1.n.star[,2:max.h] <- F.n_1.n.star[,2:max.h] -
    vec.1.T%*%matrix(mat.GAMMA.dot[(dim(mat.GAMMA.dot)[1]),2:(dim(P)[2])],nrow=1) +
    P[,1]*( vec.1.T%*%matrix(mat.GAMMA.right.dot[(dim(mat.GAMMA.right.dot)[1]-1),1:(dim(P)[2]-1)],nrow=1)
    )
  
  
  SUM.GAMMA.P.n_1.n      <- matrix(0, T, max.h)
  SUM.GAMMA.P.n_1.n.star <- matrix(0, T, max.h)
  
  mat.indices <- make.indices(PSI.a.dot)
  
  for (n in 3:max.h){
    
    AUX.SUM.GAMMA.P.n_1.n <- P[,1:(n-2)]*vec.1.T%*%matrix(mat.GAMMA.dot.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1) +
      (1-P[,1:(n-2)])*(vec.1.T%*%matrix(mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1))
    
    AUX.SUM.GAMMA.P.n_1.n.star <- P[,1:(n-2)]*
      (
        vec.1.T%*%matrix(  mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1)+
          vec.1.T%*%matrix(mat.GAMMA.right.dot[mat.indices[(n-2),1:(n-2)]],nrow=1) +
          vec.1.T%*%matrix(mat.GAMMA.left.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1)+
          vec.1.T%*%matrix(mat.GAMMA[mat.indices[(n-2),1:(n-2)]],nrow=1)
      ) +
      (1-P[,1:(n-2)])*(vec.1.T%*%matrix(mat.GAMMA.dot[mat.indices[(n-2),1:(n-2)]+1],nrow=1))
    
    SUM.GAMMA.P.n_1.n[,n]      <- AUX.SUM.GAMMA.P.n_1.n%*%matrix(1,n-2,1)
    SUM.GAMMA.P.n_1.n.star[,n] <- AUX.SUM.GAMMA.P.n_1.n.star%*%matrix(1,n-2,1)
  }
  
  F.n_1.n        <- F.n_1.n        - SUM.GAMMA.P.n_1.n
  F.n_1.n.with.u <- F.n_1.n.with.u - SUM.GAMMA.P.n_1.n
  F.n_1.n.star   <- F.n_1.n.star   - SUM.GAMMA.P.n_1.n.star
  
  cumsum.f.n_1.n        <- t(apply(F.n_1.n, 1, cumsum))
  cumsum.f.n_1.n.with.u <- cumsum.f.n_1.n - F.n_1.n + F.n_1.n.with.u
  cumsum.f.n_1.n.star <- t(apply(F.n_1.n.star, 1, cumsum))
  E.n        <- exp(-cumsum.f.n_1.n)
  E.n.with.u <- exp(-cumsum.f.n_1.n.with.u)
  E.n.star   <- exp(-cumsum.f.n_1.n.star)
  
  return(list(E.n = E.n,
              E.n.with.u = E.n.with.u,
              E.n.star = E.n.star,
              all.a.h = all.a.h,
              all.b.h = all.b.h))
}



