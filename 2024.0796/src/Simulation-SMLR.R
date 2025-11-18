# Parameter setting -------------------------------------------------------


n <- c(50, 100, 200)
d <- c(kronecker(n, c(2, 5, 10)))
k <- c(3, 5, 10)
m <- c(3, 5, 6, 10, 20)
para <- cbind(rep(rep(n, each = 3), 15), rep(d, 15), expand.grid(k, m)[c(rep(1:15, each = 9)), ])
para <- as.matrix(para)
id <- which((para[, 3]-para[,4]) > 0)
para <- para[-id, ]

id <- which(para[, 2] < 1000)
para <- para[id, ]
row.names(para) <- c(1:66)
colnames(para) <- c("n", "d", "k", "m")


# Simulation -------------------------------------------------------------------
sim_func <- function(pa){
  # setwd("~/package")
  # .libPaths(c('~/package',.libPaths()))
  library(glmnet)
  library(msgl)
  library(mvtnorm)
  library(nnet)
  library(e1071)
  
  n <- 50
  d <- pa
  k <- 3
  m <- 3
  
  # Data generation ---------------------------------------------------------
  sim_gen <- function(n, d, k, n_te, m){
    theta <- matrix(0, nr = d, nc = k)
    for (i in 1:m) {
      theta[c(1:5) + 5*(i-1), ] <- -1/(k-1)
      theta[c(1:5) + 5*(i-1), (i%%k+1)] <- 1
      
    }
    #data
    Sigma <- diag(1, d)
    # Sigma[c(1:5), c(1:5)] <- 0.2
    # Sigma[c(6:10), c(6:10)] <- 0
    # Sigma[c(11:15), c(11:15)] <- -0.2
    # diag(Sigma) <- 1
    mean <- rep(0, d)
    
    x_tr <- rmvnorm(n, mean, Sigma)
    y_tr <- NULL
    prob <- exp(x_tr %*% theta)
    
    for (i in 1:n) {
      y_tr[i] <- which.max(prob[i, ])
    }
    
    x_tu <- rmvnorm(n, mean, Sigma)
    y_tu <- NULL
    prob <- exp(x_tu %*% theta)
    for (i in 1:n) {
      y_tu[i] <- which.max(prob[i, ])
    }
    
    x_te <- rmvnorm(n_te, mean, Sigma)
    y_te <- NULL
    prob <- exp(x_te %*% theta)
    for (i in 1:n_te) {
      y_te[i] <- which.max(prob[i, ])
    }
    
    list(x_tr = x_tr, y_tr = y_tr, x_tu = x_tu, y_tu = y_tu,x_te = x_te, y_te = y_te)
  }
  # Basis function ----------------------------------------------------------
  lammax = function(X,Y,k,group,dg){ 
    
    n = nrow(X)
    ng = length(dg)
    ny = colSums(Y)
    
    lam = rep(0,ng)
    nc = as.matrix(colSums(Y)/n)
    bxi = matrix(1,n) %x% t(nc) - Y
    
    M = t(X) %*% bxi/n
    for (t in 1:ng){
      id = which(group==t)
      subm = M[id,]
      subm = as.matrix(subm)
      lam[t] = norm(subm,type="F")/sqrt(dg[t])
    }
    
    lammax = max(lam)
    return(lammax)
  }
  hess_eig = function(X,group,dg){ 
    n = nrow(X)
    ng = length(dg)
    hess = rep(0,ng)
    
    M = t(X) %*% X/n
    for (t in 1:ng){
      id = which(group==t)
      subm = M[id,id]#Ht
      tem = eigen(subm)$values
      hess[t] = max(tem)*(1+1e-6)
    }
    return(hess)
  }
  pen<- function(B,group,dg){
    ng = length(dg)
    pen = 0
    for (t in 1:ng){
      id = which(group==t)
      sB = B[id,]
      sB = as.matrix(sB)
      pen = pen + sqrt(dg[t])*norm(sB,type="F")
    }
    return(pen)
  }
  XI.gen = function(k){
    XI = matrix(0,k-1,k)
    
    XI[,1]=rep((k-1)^(-1/2),k-1)
    for (ii in 2:k)
    {
      XI[,ii]=rep( -(1+sqrt(k))/((k-1)^(1.5)), k-1)
      XI[ii-1,ii]=XI[ii-1,ii]+sqrt(k/(k-1))
    }
    
    XI = XI*sqrt(k-1)/sqrt(k)
    return(XI)
  }
  obj<- function(X,Y,B,b){
    k = ncol(Y)
    n = nrow(Y)
    V = XI.gen(k)
    t0 = X %*% B %*%V+ matrix(1,n)%*%(t(b)%*%V)
    t1 = rowSums(exp(t0))
    obj = sum(log(t1))-sum(t0*Y)
    
    obj = obj/n
    return(obj)
  }
  pb = function(X,B,b){
    k = length(b)+1
    V = XI.gen(k)
    n = nrow(X)
    score = X %*% B %*%V+ matrix(1,n)%*%(t(b)%*%V)
    P = exp(score)
    P = P/rowSums(P)
    return(P)
  }
  dual <- function(P){
    n = nrow(P)
    k = ncol(P)
    dual = 0
    for (i in 1:n){
      for (j in 1:k){
        dual = dual + P[i,j]*log(P[i,j])
      }
    }
    dual = dual/n
    return(dual)
  }
  shrink <- function(x,lam){
    n = length(x)
    t0 = rep(0,n)
    
    for(i in 1:n){
      tm = x[i]
      t0[i] = sign(tm)*max(abs(tm)-lam,0)
    }
    return(t0)
  }
  lammaxlasso = function(X,Y){ 
    n = nrow(X)
    d = ncol(X)
    k = ncol(Y)
    V = XI.gen(k)
    lam = rep(0,d)
    nc = as.matrix(colSums(Y)/n)
    bxi = matrix(1,n) %*% t(nc) - Y
    M = t(X) %*% bxi %*% t(V)/n
    lammax = max(abs(M))
    return(lammax)
  }
  pred <- function(prob){
    val = min(which.max(prob))
    return(val)
  }
  
  # SMLR based methods ------------------------------------------------------
  pen_gmlr <- function(B,b,group,dg){
    ng = length(dg)
    pen = norm(b,type="F")
    for (t in 1:ng){
      id = which(group==t)
      sB = B[id,]
      sB = as.matrix(sB)
      pen = pen + sqrt(dg[t])*norm(sB,type="F")
    }
    return(pen)
  }
  pb_gmlr = function(X,B,b){
    k = length(b)
    n = nrow(X)
    score = X %*% B+ matrix(1,n)%x%(t(b))
    P = exp(score)
    P = P/rowSums(P)
    return(P)
  }
  gmlr_path<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    obj_gmlr <- function(X,Y,B,b){
      k = ncol(Y)
      n = nrow(Y)
      t0 = X %*% B+ matrix(1,n)%x%(t(b))
      t1 = rowSums(exp(t0))
      obj = sum(log(t1))-sum(t0*Y)
      
      obj = obj/n
      return(obj)
    }
    k = ncol(Y)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    B0 = matrix(0,d,k)
    b0 = as.matrix(log(colSums(Y)/n))
    
    # a1 = log(colSums(Y)/n)
    # a1 = a1- mean(a1)
    # a1 = as.matrix(a1,k,1)
    # b0 = a1
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      tobj = rep(0,maxiters)
      tobj[1] = obj_gmlr(X,Y,B0,b0) + pen_gmlr(B0,b0,group,dg)*lambda
      
      err = Inf
      iter = 0
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:ng){
          P = pb_gmlr(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb_gmlr(X,B0,b0)
        bxi = P-Y
        u = t(bxi) %*% matrix(1,n)/n
        b0 = (b0 - 2*u) / 2
        tem = obj_gmlr(X,Y,B0,b0)+ pen_gmlr(B0,b0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0 - Bbold,type = "F")
        if (display)
          print(paste("Rand: Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb_gmlr(X,B0,b0)
      
      if (display)
        print(paste("Rand: lambda ", ii, ", err: ", err, sep =""))
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  rsmlr_bcd<- function(X,Y,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    d = ncol(X)
    n = nrow(X)
    hess = t(X)%*%X/n
    
    lambdas = lambda
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V%*%a1
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+lambda*norm(B0,type="F")^2
      
      err = Inf
      iter = 0
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:d){
          P = pb(X,B0,b0)
          bxi = P-Y
          Xt = X[,t]
          Ut = t(Xt)%*%bxi%*%t(V)/n 
          
          ht = hess[t,t]/2
          B0[t,] = (ht*B0[t,]-Ut)/(ht+2*lambda)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+lambda*norm(B0,type="F")^2
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      
      if (display)
        print(paste("Iteration ", ii, ", err: ", err, sep =""))
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  lsmlr_bcd <- function(X,Y,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    d = ncol(X)
    n = nrow(X)
    
    hess = t(X)%*%X/n
    lambdas = lambda
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V %*% a1
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+lambda*sum(abs(B0))
      
      err = Inf
      iter = 0
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:d){
          P = pb(X,B0,b0)
          bxi = P-Y
          Xt = X[,t]
          Ut = t(Xt)%*%bxi%*% t(V)/n
          
          ht = hess[t,t]/2
          B0[t,] = shrink(ht*B0[t,]-Ut,lambda)/ht
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+lambda*sum(abs(B0))
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      if (display)
        print(paste("lasso: Iteration ", ii, ", err: ", err, sep =""))
      
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  # random + GCD
  gsmlr_rand<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      B0 = matrix(0,d,k-1)  
      a1 = log(colSums(Y)/n)
      a1 = a1- mean(a1)
      a1 = as.matrix(a1,k,1)
      b0 = V%*%a1
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+pen(B0,group,dg)*lambda
      
      err = Inf
      iter = 0
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:ng){
          P = pb(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi%*%t(V)/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+ pen(B0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Rand: Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      
      if (display)
        print(paste("Rand: lambda ", ii, ", err: ", err, sep =""))
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  # path + GCD
  gsmlr_path<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda 
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V%*%a1
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+pen(B0,group,dg)*lambda
      
      err = Inf
      iter = 0
      
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:ng){
          P = pb(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi%*%t(V)/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+ pen(B0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      if(display){
        print(paste("Path: lambda ", ii, ", err: ", err, sep ="")) 
      }
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  # path + gap + GCD
  gsmlr_gap<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda 
    
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    rr = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V%*%a1
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      if (ii>1){
        P0 = pb(X,B0,b0)
        rt = lambdas[ii]/lambdas[ii-1]
        bxi0 = P0-Y
        A = log(P0)/n
        
        # difference here from refined gap
        r1 = dual(rt*bxi0+Y)-dual(P0)
        r0 = sqrt(2*n*r1)
        
        # screening step
        index = rep(1,ng)
        U0 = t(X)%*%bxi0/n
        
        for (t in 1:ng){
          id = which(group==t)
          Ut = U0[id,]
          Ut = as.matrix(Ut)
          value = lambda*sqrt(dg[t])-sqrt(hess[t]/n)*r0-norm(Ut,type = "F")
          if(value>0){
            B0[id,] = 0
            index[t] = 0
          }
          
        }
        ng1 = which(index>0)
        # cat("lambda",ii,ng-length(ng1),"groups are screened out,\n")
      }else{
        ng1 = c(1:ng)
      }
      
      
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+pen(B0,group,dg)*lambda
      
      err = Inf
      iter = 0
      
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        
        for (t in ng1){
          P = pb(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi%*%t(V)/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+ pen(B0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      
      if (ii>1){
        rr[ii] = length(ng1)
      }
      if (display)
        print(paste("Gap: lambda ", ii, ", err: ", err, sep =""))
      
    }
    
    list(B=rB, b=rb, obj=robj, fitted=fitp,lambda=lambdas,time = ct, active = rr)
  }
  # path + refined gap + GCD
  gsmlr_gapref<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda 
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    rr = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V%*%a1
    
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      if (ii>1){
        P0 = pb(X,B0,b0)
        rt = lambdas[ii]/lambdas[ii-1]
        bxi0 = P0-Y
        A = log(P0)/n
        r1 = dual(rt*bxi0+Y)-dual(P0)+(1-rt)*sum(diag(A%*%t(bxi0)))
        r0 = sqrt(2*n*r1)
        
        # screening step
        index = rep(1,ng)
        U0 = t(X)%*%bxi0/n
        
        for (t in 1:ng){
          id = which(group==t)
          Ut = U0[id,]
          Ut = as.matrix(Ut)
          value = lambda*sqrt(dg[t])-sqrt(hess[t]/n)*r0-norm(Ut,type = "F")
          if(value>0){
            B0[id,] = 0
            index[t] = 0
          }
        }
        ng1 = which(index>0)
        # cat("lambda",ii,ng-length(ng1),"groups are screened out,\n")
      }else{
        ng1 = c(1:ng)
      }
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+pen(B0,group,dg)*lambda
      
      err = Inf
      iter = 0
      
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        
        for (t in ng1){
          P = pb(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi%*%t(V)/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+ pen(B0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      
      if (ii>1){
        rr[ii] = length(ng1)
      }
      
      # gB = gnorm(B0,group,dg)
      # ng2 = which(gB>0)
      # rr[ii] = length(ng1)/length(ng2)
      if (display)
        print(paste("Ref: lambda ", ii, ", err: ", err, sep =""))
    }
    
    list(B=rB, b=rb, obj=robj, fitted=fitp,lambda=lambdas,time = ct, active = rr)
  }
  #source('C:/Users/sli16/Dropbox/Simulation-Fast Multinomial Logistic Regression with Group Sparsity/Code/Simulation-func.R')
  
  # main function -----------------------------------------------------------
  
  n_te <- 10000
  group <- rep(1:(d/5), each = 5)
  dg <- rep(5, (d/5))
  epsilon <- 1e-4
  maxiters <- 1000
  nlambda <- 100
  display = "F"
  
  #res_100 <- list()
  #for (rep in 1:100) {
  data <- sim_gen(n, d, k, n_te, m)
  x_tr <- data$x_tr
  y_tr <- data$y_tr
  x_tu <- data$x_tu
  y_tu <- data$y_tu
  x_te <- data$x_te
  y_te <- data$y_te
  Y <- matrix(0, nr = nrow(x_tr), nc = k)
  for (i in 1:nrow(x_tr)) {
    Y[i, y_tr[i]] <- 1
  }
  X <- x_tr
  
  lam0 <- max(lammaxlasso(X,Y),
              lammax(X,Y,k,group,dg),
              lambda(x_tr, y_tr, alpha = 0, d = 2, lambda.min = 0.05)[1])
  lambda.min.ratio <- ifelse(n < d, 0.01, 1e-04)
  lambda <- lam0*seq(1, lambda.min.ratio, length.out=nlambda)
  
  
  time <- list()
  error <- list()
  SF <- list()
  
  #MLR
  data_fra <- data.frame(class = y_tr, x_tr)
  start_time <- Sys.time()
  fit_m <- multinom(class ~ ., data = data_fra)
  end_time <- Sys.time()
  pre_res <- predict(fit_m, newdata = data.frame(x_te))
  time[['MLR']] <- end_time - start_time
  error[['MLR']] <- length(which(pre_res != y_te)) / length(y_te)
  SF[['MLR']] <- d
  
  #SVM
  start_time <- Sys.time()
  fit_svm <- svm(x_tr, as.factor(y_tr))
  end_time <- Sys.time()
  pre_res <- predict(object = fit_svm, newdata = x_te)
  time[['SVM']] <- end_time - start_time
  error[['SVM']] <- length(which(pre_res != y_te)) / length(y_te)
  SF[['SVM']] <- d
  
  #ANN
  start_time <- Sys.time()
  fit_ann <- nnet(x_tr, Y, size = 1)
  end_time <- Sys.time()
  pre_res <- max.col(predict(fit_ann, x_te))
  time[['ANN']] <- end_time - start_time
  error[['ANN']] <- length(which(pre_res != y_te)) / length(y_te)
  SF[['ANN']] <- d
  
  #MLR+Lasso
  start_time <- Sys.time()
  fit_ml <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 1, type.measure = "class", lambda = lambda, intercept = T)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    err_tu[i] <- length(which(as.numeric(predict(fit_ml, newx = x_tu, s = fit_ml$lambda[i], type = "class")) != y_tu))
  }
  index <- which.min(err_tu)
  
  pre_res <- as.numeric(predict(fit_ml, newx = x_te, s = fit_ml$lambda[index], type = "class"))
  time[['LMLR']] <- end_time - start_time
  error[['LMLR']] <- length(which(pre_res != y_te)) / length(y_te)
  beta_lmlr <- 0
  for (i in 1:k) {
    beta_lmlr <- beta_lmlr + abs(c(fit_ml$beta[[i]][, index])) 
  }
  SF[['LMLR']] <- as.numeric(which(beta_lmlr > 0))
  
  
  #MLR+Ridge
  start_time <- Sys.time()
  fit_mr <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 0, type.measure = "class", lambda = lambda, intercept = T)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    err_tu[i] <- length(which(as.numeric(predict(fit_mr, newx = x_tu, s = fit_mr$lambda[i], type = "class")) != y_tu))
  }
  index <- which.min(err_tu)
  pre_res <- as.numeric(predict(fit_mr, newx = x_te, s = fit_mr$lambda[index], type = "class"))
  time[['RMLR']] <- end_time - start_time
  error[['RMLR']] <- length(which(pre_res != y_te)) / length(y_te)
  beta_rmlr <- 0
  for (i in 1:k) {
    beta_rmlr <- beta_rmlr + abs(c(fit_mr$beta[[i]][, index])) 
  }
  SF[['RMLR']] <- as.numeric(which(beta_rmlr > 0))
  
  #MLR+GL
  
  start_time <- Sys.time()
  fit_mg <- fit(x_tr, y_tr, alpha = 0, lambda = lambda, grouping = group, groupWeights = rep(sqrt(5), (d/5)), intercept = F)
  end_time <- Sys.time()
  res_tu <- predict(fit_mg, x_tu)
  index <- which.min(Err(res_tu, classes = y_tu))
  res_pre <- predict(fit_mg, x_te)
  time[['GMLR']] <- end_time - start_time
  error[['GMLR']] <- length(which(apply(res_pre$response[[index]], 2, pred) != y_te)) / length(y_te)
  SF[['GMLR']] <- features(fit_mg)[[index]]
  
  #SMLR+Lasso
  start_time <- Sys.time()
  fit_sl <- lsmlr_bcd(X, Y, epsilon, maxiters, nlambda, display, lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_sl$B[[i]], fit_sl$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['LSMLR']] <- end_time - start_time
  error[['LSMLR']] <- length(which(apply(pb(x_te, fit_sl$B[[index]], fit_sl$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['LSMLR']] <- which(rowSums(abs(fit_sl$B[[index]])) > 0)
  
  #SMLR+Ridge
  start_time <- Sys.time()
  fit_sr <- rsmlr_bcd(X, Y, epsilon, maxiters, nlambda, display, lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_sr$B[[i]], fit_sr$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['RSMLR']] <- end_time - start_time
  error[['RSMLR']] <- length(which(apply(pb(x_te, fit_sr$B[[index]], fit_sr$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['RSMLR']] <- which(rowSums(abs(fit_sr$B[[index]])) > 0)
  
  #GSMLR+random GCD
  start_time <- Sys.time()
  fit_grg <- gsmlr_rand(X,Y,group,dg,epsilon,maxiters,nlambda,display, lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_grg$B[[i]], fit_grg$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['GSMLR_rand']] <- end_time - start_time
  error[['GSMLR_rand']] <- length(which(apply(pb(x_te, fit_grg$B[[index]], fit_grg$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['GSMLR_rand']] <- which(rowSums(abs(fit_grg$B[[index]])) >0)
  
  #GMLR+path GCD
  start_time <- Sys.time()
  fit_gmlr <- gmlr_path(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb_gmlr(x_tu, fit_gmlr$B[[i]], fit_gmlr$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['GMLR_path']] <- end_time - start_time
  error[['GMLR_path']] <- length(which(apply(pb_gmlr(x_te, fit_gmlr$B[[index]], fit_gmlr$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['GMLR_path']] <- which(rowSums(abs(fit_gmlr$B[[index]])) > 0)
  
  #GSMLR+path GCD
  start_time <- Sys.time()
  fit_gpg <- gsmlr_path(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_gpg$B[[i]], fit_gpg$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['GSMLR_path']] <- end_time - start_time
  error[['GSMLR_path']] <- length(which(apply(pb(x_te, fit_gpg$B[[index]], fit_gpg$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['GSMLR_path']] <- which(rowSums(abs(fit_gpg$B[[index]])) > 0)
  
  #GSMLR+path GAP GCD
  start_time <- Sys.time()
  fit_gpgg <- gsmlr_gap(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_gpgg$B[[i]], fit_gpgg$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['GSMLR_gap']] <- end_time - start_time
  error[['GSMLR_gap']] <- length(which(apply(pb(x_te, fit_gpgg$B[[index]], fit_gpgg$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['GSMLR_gap']] <- which(rowSums(abs(fit_gpgg$B[[index]])) > 0)
  
  #GSMLR+path refined GAP GCD
  start_time <- Sys.time()
  fit_gprgg <- gsmlr_gapref(X, Y, group, dg, epsilon, maxiters, nlambda, display, lambda)
  end_time <- Sys.time()
  err_tu <- NULL
  for (i in 1:nlambda) {
    prob_tu <- pb(x_tu, fit_gprgg$B[[i]], fit_gprgg$b[[i]])
    err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
  }
  index <- which.min(err_tu)
  time[['GSMLR_gapref']] <- end_time - start_time
  error[['GSMLR_gapref']] <- length(which(apply(pb(x_te, fit_gprgg$B[[index]], fit_gprgg$b[[index]]), 1, pred) != y_te)) / length(y_te)
  SF[['GSMLR_gapref']] <- which(rowSums(abs(fit_gprgg$B[[index]])) > 0)
  
  res <- list(time = time, error = error, SF = SF)
  return(res)
}
library(parallel)
cl <- detectCores()
clus <- makeCluster(cl - 1)
res_all_sim_1 <- list()
for (i in 1:4) {
  # l_1 <- 6 * (i-1) + 1
  # l_2 <- 6 * (i-1) + 6
  parm <- matrix(c(2:5) * 50, nc = 100, nr = 4, byrow = F)
  #res <- parApply(clus, X = parm, MARGIN = 1, FUN = sim_func_2)
  res <- parLapply(clus, parm[i, ], fun = sim_func)
  res_all_sim_1[[i]] <- res
  save(res_all_sim_1, file = "sim_GSMLR_1.Rdata")
  #save(res_all, file = "sim_SMLR_independent.Rdata")
  print(i)
}
stopCluster(clus)
# sim_1 ---------------------------------------------------------------
res_time <- list()
for (i in 1:4) {
  res_time[[i]] <- list()
  for (j in 1:13) {
    res_time[[i]][[j]] <- 0
  }
  names(res_time[[i]]) <- names(res_all_sim_1[[1]][[1]][["time"]])
}
id <- NULL
for (i in 1:4) {
  if(length(res_all_sim_1[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    for (k in 1:13) {
      res_time[[i]][[k]] <- res_time[[i]][[k]] + res_all_sim_1[[i]][[j]][["time"]][[k]] / 100
    }
  }
  id <- c(id, i)
}
sim_time <- NULL
for (i in id) {
  sim_time <- rbind(sim_time, c(as.numeric(res_time[[i]])))
}
colnames(sim_time) <- c(names(res_all_sim_1[[1]][[1]][["time"]]))
#, "n", "d", "k", "m")
write.csv(sim_time, file = "time_sim_1.csv")

res_error <- list()
for (i in 1:4) {
  res_error[[i]] <- list()
  for (j in 1:13) {
    res_error[[i]][[j]] <- 0
  }
  names(res_error[[i]]) <- names(res_all_sim_1[[1]][[1]][["error"]])
}
id <- NULL
for (i in 1:4) {
  if(length(res_all_sim_1[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    for (k in 1:13) {
      res_error[[i]][[k]] <- res_error[[i]][[k]] + res_all_sim_1[[i]][[j]][["error"]][[k]] / 100
    }
  }
  id <- c(id, i)
}
sim_error <- NULL
for (i in id) {
  sim_error <- rbind(sim_error, c(as.numeric(res_error[[i]])))
  #, para[i, ]))
}
colnames(sim_error) <- c(names(res_all_sim_1[[1]][[1]][["error"]]))
#, "n", "d", "k", "m")
write.csv(sim_error, file = "error_sim_1.csv")

res_SF <- list()
for (i in 1:4) {
  res_SF[[i]] <- list()
  for (j in 1:13) {
    res_SF[[i]][[j]] <- 0
  }
  names(res_SF[[i]]) <- names(res_all_sim_1[[1]][[1]][["SF"]])
}
id <- NULL
for (i in 1:4) {
  if(length(res_all_sim_1[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    for (k in 1:13) {
      res_SF[[i]][[k]] <- res_SF[[i]][[k]] + length(res_all_sim_1[[i]][[j]][["SF"]][[k]]) / 100
    }
  }
  id <- c(id, i)
}
sim_SF <- NULL
for (i in id) {
  sim_SF <- rbind(sim_SF, c(as.numeric(res_SF[[i]])))
}
colnames(sim_SF) <- c(names(res_all_sim_1[[1]][[1]][["SF"]]))
#, "n", "d", "k", "m")
write.csv(sim_SF, file = "SF_sim_1.csv")

#write.csv(sim_SF_gr, file = "SF_gr_app_recog.csv")

# sim_2 -------------------------------------------------------------------
n <- c(50, 100, 200)
d <- c(kronecker(n, c(2, 5, 10)))
k <- c(3, 5, 10)
m <- c(3, 5, 6, 10, 20)
para <- cbind(rep(rep(n, each = 3), 15), rep(d, 15), expand.grid(k, m)[c(rep(1:15, each = 9)), ])
para <- as.matrix(para)
id <- which((para[, 3]-para[,4]) > 0)
para <- para[-id, ]
id <- which(para[, 2] < 1000)
para <- para[id, ]
row.names(para) <- c(1:66)
colnames(para) <- c("n", "d", "k", "m")
para <- para[c(1:25), ]

res_error <- list()
for (i in 1:25) {
  res_error[[i]] <- list()
  for (j in 1:13) {
    res_error[[i]][[j]] <- 0
  }
  names(res_error[[i]]) <- names(res_all_sim_2[[1]][[1]][["error"]])
}
id <- NULL
for (i in 1:25) {
  if(length(res_all_sim_2[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    if(names(res_all_sim_2[[i]][[j]][["error"]])[[1]] == "MLR"){
      for (k in 1:13) {
        res_error[[i]][[k]] <- res_error[[i]][[k]] + res_all_sim_2[[i]][[j]][["error"]][[k]] / 100
      }
    }else{
      for (k in 1:12) {
        res_error[[i]][[k+1]] <- res_error[[i]][[k+1]] + res_all_sim_2[[i]][[j]][["error"]][[k]] / 100
      }
    }
    
  }
  id <- c(id, i)
}
sim_error <- NULL
for (i in id) {
  sim_error <- rbind(sim_error, c(as.numeric(res_error[[i]])))
  #, para[i, ]))
}
colnames(sim_error) <- c(names(res_all_sim_2[[1]][[1]][["error"]]))
write.csv(sim_error, file = "error_sim_2.csv")


res_SF <- list()
for (i in 1:25) {
  res_SF[[i]] <- list()
  for (j in 1:13) {
    res_SF[[i]][[j]] <- 0
  }
  names(res_SF[[i]]) <- names(res_all_sim_2[[1]][[1]][["SF"]])
}
id <- NULL
for (i in 1:25) {
  if(length(res_all_sim_2[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    if(names(res_all_sim_2[[i]][[j]][["SF"]])[[1]] == "MLR"){
      for (k in 1:13) {
        res_SF[[i]][[k]] <- res_SF[[i]][[k]] + length(res_all_sim_2[[i]][[j]][["SF"]][[k]]) / 100
      }
    }else{
      for (k in 1:12) {
        res_SF[[i]][[k+1]] <- res_SF[[i]][[k+1]] + length(res_all_sim_2[[i]][[j]][["SF"]][[k]]) / 100
      }
    }
    
  }
  id <- c(id, i)
}
sim_SF <- NULL
for (i in id) {
  sim_SF <- rbind(sim_SF, c(as.numeric(res_SF[[i]])))
  #, para[i, ]))
}
colnames(sim_SF) <- c(names(res_all_sim_2[[1]][[1]][["SF"]]))
write.csv(sim_SF, file = "SF_sim_2.csv")


res_sel <- list()
for (i in 1:25) {
  res_sel[[i]] <- list()
  for (j in 1:13) {
    res_sel[[i]][[j]] <- 0
  }
  names(res_sel[[i]]) <- names(res_all_sim_2[[1]][[1]][["SF"]])
}
id <- NULL
for (i in 1:25) {
  if(length(res_all_sim_2[[i]]) == 1){
    next
  }
  for (j in 1:100) {
    if(names(res_all_sim_2[[i]][[j]][["SF"]])[[1]] == "MLR"){
      for (k in 1:13) {
        tpr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c(1:(para[i, 4] * 5)))) / (para[i, 4] * 5)
        fpr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((para[i, 4] * 5+1):para[i, 2]))) / (para[i, 2]- para[i, 4] * 5)
        fdr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((para[i, 4] * 5+1):para[i, 2]))) / length(res_all_sim_2[[i]][[j]][["SF"]][[k]])
        
        SF_gr <- NULL
        for (l in 1:(para[i, 2]/5)) {
          if(length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((5*l-4):(5*l)))) != 0){
            SF_gr <- c(SF_gr, l)
          }
        }
        
        tpr_gr <- length(intersect(SF_gr, c(1:para[i, 4]))) / para[i, 4]
        fpr_gr <- length(intersect(SF_gr, c((para[i, 4]+1):(para[i, 2]/5)))) / (para[i, 2]/5- para[i, 4])
        fdr_gr <- length(intersect(SF_gr, c((para[i, 4]+1):(para[i, 2]/5)))) / length(SF_gr)
        
        res_sel[[i]][[k]] <- res_sel[[i]][[k]] + c(tpr, fpr, fdr, tpr_gr, fpr_gr, fdr_gr) / 100
      }
    }else{
      for (k in 1:12) {
        tpr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c(1:(para[i, 4] * 5)))) / (para[i, 4] * 5)
        fpr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((para[i, 4] * 5+1):para[i, 2]))) / (para[i, 2]- para[i, 4] * 5)
        fdr <- length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((para[i, 4] * 5+1):para[i, 2]))) / length(res_all_sim_2[[i]][[j]][["SF"]][[k]])
        
        SF_gr <- NULL
        for (l in 1:(para[i, 2]/5)) {
          if(length(intersect(res_all_sim_2[[i]][[j]][["SF"]][[k]], c((5*l-4):(5*l)))) != 0){
            SF_gr <- c(SF_gr, l)
          }
        }
        
        tpr_gr <- length(intersect(SF_gr, c(1:para[i, 4]))) / para[i, 4]
        fpr_gr <- length(intersect(SF_gr, c((para[i, 4]+1):(para[i, 2]/5)))) / (para[i, 2]/5- para[i, 4])
        fdr_gr <- length(intersect(SF_gr, c((para[i, 4]+1):(para[i, 2]/5)))) / length(SF_gr)
        
        res_sel[[i]][[k+1]] <- res_sel[[i]][[k+1]] + c(tpr, fpr, fdr, tpr_gr, fpr_gr, fdr_gr) / 100
      }
    }
  }
  id <- c(id, i)
}

sim_sel <- NULL
for (i in id) {
  for (j in 1:13) {
    sim_sel <- rbind(sim_sel, c(res_sel[[i]][[j]]))
  }
}
colnames(sim_sel) <- c("tpr", "fpr", "fdr", "tpr_gr", "fpr_gr", "fdr_gr")
rownames(sim_sel) <- rep(names(res_all_sim_2[[1]][[1]][["SF"]]), 25)
para <- cbind(para, matrix(rep(0,25 * 9), nr = 25, nc = 9))

sim_sel <- cbind(sim_sel, c(t(para)))
sim_sel <- cbind(c(t(sim_error)), sim_sel)
write.csv(sim_sel, file = "sim_GSMLR_2.csv")



# Application -------------------------------------------------------------
library("parallel")
app_func <- function(para){
  #setwd("~/package")
  #.libPaths(c('~/package',.libPaths()))
  load("isolet.Rdata")
  #load("USPS.Rdata")
  #load("madelon.Rdata")
  #load("recog_533.Rdata")
  # testdata$X <- testdata$X[, -c(100, 114, 254, 268)]
  library(glmnet)
  library(msgl)
  library(dplyr)
  # Basis function ----------------------------------------------------------
  lammax = function(X,Y,k,group,dg){ 
    
    n = nrow(X)
    ng = length(dg)
    ny = colSums(Y)
    
    lam = rep(0,ng)
    nc = as.matrix(colSums(Y)/n)
    bxi = matrix(1,n) %x% t(nc) - Y
    
    M = t(X) %*% bxi/n
    for (t in 1:ng){
      id = which(group==t)
      subm = M[id,]
      subm = as.matrix(subm)
      lam[t] = norm(subm,type="F")/sqrt(dg[t])
    }
    
    lammax = max(lam)
    return(lammax)
  }
  hess_eig = function(X,group,dg){ 
    n = nrow(X)
    ng = length(dg)
    hess = rep(0,ng)
    
    M = t(X) %*% X/n
    for (t in 1:ng){
      id = which(group==t)
      subm = M[id,id]#Ht
      tem = eigen(subm)$values
      hess[t] = max(tem)*(1+1e-6)
    }
    return(hess)
  }
  pen<- function(B,group,dg){
    ng = length(dg)
    pen = 0
    for (t in 1:ng){
      id = which(group==t)
      sB = B[id,]
      sB = as.matrix(sB)
      pen = pen + sqrt(dg[t])*norm(sB,type="F")
    }
    return(pen)
  }
  XI.gen = function(k){
    XI = matrix(0,k-1,k)
    
    XI[,1]=rep((k-1)^(-1/2),k-1)
    for (ii in 2:k)
    {
      XI[,ii]=rep( -(1+sqrt(k))/((k-1)^(1.5)), k-1)
      XI[ii-1,ii]=XI[ii-1,ii]+sqrt(k/(k-1))
    }
    
    XI = XI*sqrt(k-1)/sqrt(k)
    return(XI)
  }
  obj<- function(X,Y,B,b){
    k = ncol(Y)
    n = nrow(Y)
    V = XI.gen(k)
    t0 = X %*% B %*%V+ matrix(1,n)%*%(t(b)%*%V)
    t1 = rowSums(exp(t0))
    obj = sum(log(t1))-sum(t0*Y)
    
    obj = obj/n
    return(obj)
  }
  pb = function(X,B,b){
    k = length(b)+1
    V = XI.gen(k)
    n = nrow(X)
    score = X %*% B %*%V+ matrix(1,n)%*%(t(b)%*%V)
    P = exp(score)
    P = P/rowSums(P)
    return(P)
  }
  dual <- function(P){
    n = nrow(P)
    k = ncol(P)
    dual = 0
    for (i in 1:n){
      for (j in 1:k){
        dual = dual + P[i,j]*log(P[i,j])
      }
    }
    dual = dual/n
    return(dual)
  }
  shrink <- function(x,lam){
    n = length(x)
    t0 = rep(0,n)
    
    for(i in 1:n){
      tm = x[i]
      t0[i] = sign(tm)*max(abs(tm)-lam,0)
    }
    return(t0)
  }
  lammaxlasso = function(X,Y){ 
    n = nrow(X)
    d = ncol(X)
    k = ncol(Y)
    V = XI.gen(k)
    lam = rep(0,d)
    nc = as.matrix(colSums(Y)/n)
    bxi = matrix(1,n) %*% t(nc) - Y
    M = t(X) %*% bxi %*% t(V)/n
    lammax = max(abs(M))
    return(lammax)
  }
  pred <- function(prob){
    val = min(which.max(prob))
    return(val)
  }
  # SMLR based methods ------------------------------------------------------
  pen_gmlr <- function(B,b,group,dg){
    ng = length(dg)
    pen = norm(b,type="F")
    for (t in 1:ng){
      id = which(group==t)
      sB = B[id,]
      sB = as.matrix(sB)
      pen = pen + sqrt(dg[t])*norm(sB,type="F")
    }
    return(pen)
  }
  pb_gmlr = function(X,B,b){
    k = length(b)
    n = nrow(X)
    score = X %*% B+ matrix(1,n)%x%(t(b))
    P = exp(score)
    P = P/rowSums(P)
    return(P)
  }
  gmlr_path<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    obj_gmlr <- function(X,Y,B,b){
      k = ncol(Y)
      n = nrow(Y)
      t0 = X %*% B+ matrix(1,n)%x%(t(b))
      t1 = rowSums(exp(t0))
      obj = sum(log(t1))-sum(t0*Y)
      
      obj = obj/n
      return(obj)
    }
    k = ncol(Y)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    
    B0 = matrix(0,d,k)
    b0 = as.matrix(log(colSums(Y)/n))
    
    # a1 = log(colSums(Y)/n)
    # a1 = a1- mean(a1)
    # a1 = as.matrix(a1,k,1)
    # b0 = a1
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      tobj = rep(0,maxiters)
      tobj[1] = obj_gmlr(X,Y,B0,b0) + pen_gmlr(B0,b0,group,dg)*lambda
      
      err = Inf
      iter = 0
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        for (t in 1:ng){
          P = pb_gmlr(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb_gmlr(X,B0,b0)
        bxi = P-Y
        u = t(bxi) %*% matrix(1,n)/n
        b0 = (b0 - 2*u) / 2
        tem = obj_gmlr(X,Y,B0,b0)+ pen_gmlr(B0,b0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0 - Bbold,type = "F")
        if (display)
          print(paste("Rand: Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb_gmlr(X,B0,b0)
      
      if (display)
        print(paste("Rand: lambda ", ii, ", err: ", err, sep =""))
    }
    
    list(B = rB, b = rb,obj = robj,fitted = fitp,lambda=lambdas,time = ct)
  }
  # path + refined gap + GCD
  gsmlr_gapref<- function(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda){
    k = ncol(Y)
    V = XI.gen(k)
    ng = length(dg)
    d = ncol(X)
    n = nrow(X)
    
    hess = hess_eig(X,group,dg)
    lambdas = lambda 
    
    rB = list()
    rb = list()
    robj = list()
    fitp = list()
    ct = rep(0,nlambda)
    rr = rep(0,nlambda)
    
    B0 = matrix(0,d,k-1)  
    a1 = log(colSums(Y)/n)
    a1 = a1- mean(a1)
    a1 = as.matrix(a1,k,1)
    b0 = V%*%a1
    
    
    
    for (ii in 1:nlambda){
      lambda = lambdas[ii]
      
      if (ii>1){
        P0 = pb(X,B0,b0)
        rt = lambdas[ii]/lambdas[ii-1]
        bxi0 = P0-Y
        A = log(P0)/n
        r1 = dual(rt*bxi0+Y)-dual(P0)+(1-rt)*sum(diag(A%*%t(bxi0)))
        r0 = sqrt(2*n*r1)
        
        # screening step
        index = rep(1,ng)
        U0 = t(X)%*%bxi0/n
        
        for (t in 1:ng){
          id = which(group==t)
          Ut = U0[id,]
          Ut = as.matrix(Ut)
          value = lambda*sqrt(dg[t])-sqrt(hess[t]/n)*r0-norm(Ut,type = "F")
          if(value>0){
            B0[id,] = 0
            index[t] = 0
          }
        }
        ng1 = which(index>0)
        # cat("lambda",ii,ng-length(ng1),"groups are screened out,\n")
      }else{
        ng1 = c(1:ng)
      }
      
      tobj = rep(0,maxiters)
      tobj[1] = obj(X,Y,B0,b0)+pen(B0,group,dg)*lambda
      
      err = Inf
      iter = 0
      
      start_time <- Sys.time()
      while ((err>=epsilon) & (iter<maxiters)){
        ll1 <- tobj[iter]
        iter <- iter + 1
        Bbold = rbind(B0,t(b0))
        
        for (t in ng1){
          P = pb(X,B0,b0)
          bxi = P-Y
          id = which(group==t)
          Xt = X[,id]
          Ut = t(Xt)%*%bxi%*%t(V)/n
          
          ht = hess[t]
          A1 = ht/2*B0[id,]-Ut
          A1 = as.matrix(A1)
          a = norm(A1,type = "F")
          a1 = max(0,1-lambda*sqrt(dg[t])/a)
          B0[id,] = a1*(B0[id,]-Ut*2/ht)
        }
        P = pb(X,B0,b0)
        bxi = P-Y
        u = V%*%t(bxi)%*%matrix(1,n)/n
        b0 = b0 - 2*u
        tem = obj(X,Y,B0,b0)+ pen(B0,group,dg)*lambda
        
        tobj[iter] <- tem
        ll2 <- tem
        Bb0 = rbind(B0,t(b0))
        err = norm(Bb0-Bbold,type = "F")
        if (display)
          print(paste("Iteration ", iter, ", err: ", err, sep =""))
      }
      end_time <- Sys.time()
      ct[ii] = end_time - start_time
      
      
      ## ----------------------------------------## 
      ## Summary
      ## ----------------------------------------## 
      
      rB[[ii]] = B0
      rb[[ii]] = b0
      robj[[ii]] = tobj[1:iter]
      fitp[[ii]] = pb(X,B0,b0)
      
      if (ii>1){
        rr[ii] = length(ng1)
      }
      
      # gB = gnorm(B0,group,dg)
      # ng2 = which(gB>0)
      # rr[ii] = length(ng1)/length(ng2)
      if (display)
        print(paste("Ref: lambda ", ii, ", err: ", err, sep =""))
    }
    
    list(B=rB, b=rb, obj=robj, fitted=fitp,lambda=lambdas,time = ct, active = rr)
  }
  #source('C:/Users/sli16/Dropbox/Simulation-Fast Multinomial Logistic Regression with Group Sparsity/Code/Simulation-func.R')
  # data&setting --------------------------------------------------------------------
  res_clu <- kmeans(t(testdata$X), 9)
  table(res_clu$cluster)
  group <- as.numeric(res_clu$cluster)
  id_group <- list()
  for (i in 1:max(group)) {
    id_group[[i]] <- which(group == i)
  }
  dg <- NULL
  for (i in 1:max(group)) {
    dg[i] <- length(which(group == i))
  }
  data_recog <- cbind(testdata$X, testdata$Y)
  n <- ncol(data_recog)
  k <- max(testdata$Y)
  
  
  
  id_tr <- group_by(data.frame(testdata$X), data.frame(testdata$Y))
  data_tr <- sample_frac(id_tr, size = 60 / length(testdata$Y))
  data_rest <- setdiff(id_tr, data_tr)
  id_tu <- group_by(data_rest)
  data_tu <- sample_frac(id_tu, size = 60 / length(data_rest$testdata.Y))
  data_te <- setdiff(id_tu, data_tu)
  
  data_tr <- as.matrix(data_tr)
  data_te <- as.matrix(data_te)
  data_tu <- as.matrix(data_tu)
  x_tr <- data_tr[, -n]
  y_tr <- data_tr[, n]
  x_te <- data_te[, -n]
  y_te <- data_te[, n]
  x_tu <- data_tu[, -n]
  y_tu <- data_tu[, n]
  
  Y <- matrix(0, nr = nrow(x_tr), nc = max(testdata$Y))
  for (i in 1:nrow(x_tr)) {
    Y[i, y_tr[i]] <- 1
  }
  X <- x_tr
  epsilon <- 1e-4
  maxiters <- 500
  nlambda <- 100
  display = "F"
  error <- list()
  SF <- list()
  SF_gr <- list()
  # main --------------------------------------------------------------------
  main <- function(para){
    #lambda-sequence
    lam0 <- max(lammaxlasso(X,Y),
                lammax(X,Y,k,group,dg),
                lambda(x_tr, y_tr, alpha = 0, d = 2, lambda.min = 0.001, standardize = F)[1])
    lambda.min.ratio <- 0.1
    lambda <- lam0*seq(1, lambda.min.ratio, length.out=nlambda) * 1.5
    
    
    #MLR+Lasso
    fit_ml <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 1, lambda = lambda)
    err_ml <- NULL
    for (i in 1:length(fit_ml$lambda)) {
      err_ml[i] <- length(which(as.numeric(predict(fit_ml, newx = x_te, s = fit_ml$lambda[i], type = "class")) != y_te)) / length(y_te)
    }
    error[['LMLR']] <- err_ml
    for (i in 1:length(fit_ml$lambda)) {
      beta_lmlr <- 0
      for (j in 1:length(fit_ml$beta)) {
        beta_lmlr <- beta_lmlr + abs(c(fit_ml$beta[[j]][, i]))
      }
      SF[['LMLR']][i] <- length(as.numeric(which(beta_lmlr != 0)))
      grou_sel <- group[as.numeric(which(beta_lmlr != 0))]
      SF_gr[['LMLR']][i] <- length(grou_sel[!duplicated(grou_sel)])
    }
    
    
    
    # #MLR+Ridge
    # fit_mr <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 0, lambda = lambda)
    # err_tu <- NULL
    # for (i in 1:length(fit_mr$lambda)) {
    #   err_tu[i] <- length(which(as.numeric(predict(fit_mr, newx = x_tu, s = fit_mr$lambda[i], type = "class")) != y_tu))
    # }
    # index <- which.min(err_tu)
    # pre_res <- as.numeric(predict(fit_mr, newx = x_te, s = fit_mr$lambda[index], type = "class"))
    # error[['RMLR']] <- length(which(pre_res != y_te)) / length(y_te)
    # beta_rmlr <- 0
    # for (i in 1:length(fit_mr$beta)) {
    #   beta_rmlr <- beta_rmlr + abs(c(fit_mr$beta[[i]][, index]))
    # }
    # SF[['RMLR']] <- as.numeric(which(beta_rmlr != 0))
    
    #MLR+GL
    fit_mg <- fit(x_tr, y_tr, alpha = 0, lambda = lambda, grouping = group, groupWeights = sqrt(dg), intercept = T, standardize = F)
    pre_te <- predict(fit_mg, x_te)
    err_gl <- Err(pre_te, classes = y_te)
    error[['GMLR']] <- err_gl
    for (i in 1:nlambda) {
      SF[['GMLR']][i] <- length(as.numeric(which(colSums(fit_mg$beta[[i]][, -1]) != 0)))
      grou_sel <- group[as.numeric(which(colSums(fit_mg$beta[[i]][, -1]) != 0))]
      SF_gr[['GMLR']][i] <- length(grou_sel[!duplicated(grou_sel)])
    }
    
    #GMLR+path GCD
    fit_gmlr <- gmlr_path(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
    err_gmlr <- NULL
    for (i in 1:nlambda) {
      prob_te <- pb_gmlr(x_te, fit_gmlr$B[[i]], fit_gmlr$b[[i]])
      err_gmlr[i] <- length(which(apply(prob_te, 1, pred) != y_te)) / length(y_te)
      SF[['GMLR_path']][i] <- length(which(rowSums(abs(fit_gmlr$B[[i]])) != 0))
      grou_gmlr <- group[which(rowSums(abs(fit_gmlr$B[[i]])) != 0)]
      SF_gr[['GMLR_path']][i] <- length(grou_gmlr[!duplicated(grou_gmlr)])
    }
    error[['GMLR_path']] <- err_gmlr
    
    
    
    #SGMLR
    fit_gprgg <- gsmlr_gapref(X, Y, group, dg, epsilon, maxiters, nlambda, display, lambda)
    err_sgl <- NULL
    for (i in 1:nlambda) {
      prob_te <- pb(x_te, fit_gprgg$B[[i]], fit_gprgg$b[[i]])
      err_sgl[i] <- length(which(apply(prob_te, 1, pred) != y_te)) / length(y_te)
      SF[['GMLR_gapref']][i] <- length(which(rowSums(abs(fit_gprgg$B[[i]])) != 0))
      grou_sel <- group[which(rowSums(abs(fit_gprgg$B[[i]])) != 0)]
      SF_gr[['GMLR_gapref']][i] <- length(grou_sel[!duplicated(grou_sel)])
    }
    error[['GMLR_gapref']] <- err_sgl
    
    # #MLR+Lasso
    # fit_ml <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 1, type.measure = "class", lambda = lambda, intercept = T)
    # err_tu <- NULL
    # for (i in 1:nlambda) {
    #   err_tu[i] <- length(which(as.numeric(predict(fit_ml, newx = x_tu, s = fit_ml$lambda[i], type = "class")) != y_tu))
    # }
    # index <- which.min(err_tu)
    # pre_res <- as.numeric(predict(fit_ml, newx = x_te, s = fit_ml$lambda[index], type = "class"))
    # error[['LMLR']] <- length(which(pre_res != y_te)) / length(y_te)
    # beta_lmlr <- 0
    # for (i in 1:k) {
    #   beta_lmlr <- beta_lmlr + abs(c(fit_ml$beta[[i]][, index])) 
    # }
    # SF[['LMLR']] <- as.numeric(which(beta_lmlr != 0))
    # grou_sel <- group[as.numeric(which(beta_lmlr != 0))]
    # SF_gr[['LMLR']] <- grou_sel[!duplicated(grou_sel)]
    # 
    # #MLR+Ridge
    # fit_mr <- glmnet(x_tr, y_tr, family = "multinomial", alpha = 0, type.measure = "class", lambda = lambda, intercept = T)
    # err_tu <- NULL
    # for (i in 1:nlambda) {
    #   err_tu[i] <- length(which(as.numeric(predict(fit_mr, newx = x_tu, s = fit_mr$lambda[i], type = "class")) != y_tu))
    # }
    # index <- which.min(err_tu)
    # pre_res <- as.numeric(predict(fit_mr, newx = x_te, s = fit_mr$lambda[index], type = "class"))
    # error[['RMLR']] <- length(which(pre_res != y_te)) / length(y_te)
    # beta_rmlr <- 0
    # for (i in 1:k) {
    #   beta_rmlr <- beta_rmlr + abs(c(fit_mr$beta[[i]][, index])) 
    # }
    # SF[['RMLR']] <- as.numeric(which(beta_rmlr != 0))
    # grou_sel <- group[as.numeric(which(beta_rmlr != 0))]
    # SF_gr[['RMLR']] <- grou_sel[!duplicated(grou_sel)]
    # 
    # #MLR+GL
    # fit_mg <- fit(x_tr, y_tr, alpha = 0, lambda = lambda, grouping = group, groupWeights = sqrt(dg), intercept = T, standardize = F)
    # res_tu <- predict(fit_mg, x_tu)
    # index <- which.min(Err(res_tu, classes = y_tu))
    # res_pre <- predict(fit_mg, x_te)
    # error[['GMLR']] <- length(which(apply(res_pre$response[[index]], 2, pred) != y_te)) / length(y_te)
    # SF[['GMLR']] <- as.numeric(which(colSums(fit_mg$beta[[index]][, -1]) != 0))
    # grou_sel <- group[as.numeric(which(colSums(fit_mg$beta[[index]][, -1]) != 0))]
    # SF_gr[['GMLR']] <- grou_sel[!duplicated(grou_sel)]
    # 
    # #GMLR+path GCD
    # fit_gmlr <- gmlr_path(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
    # err_tu <- NULL
    # for (i in 1:nlambda) {
    #   prob_tu <- pb_gmlr(x_tu, fit_gmlr$B[[i]], fit_gmlr$b[[i]])
    #   err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
    # }
    # index <- which.min(err_tu)
    # error[['GMLR_path']] <- length(which(apply(pb_gmlr(x_te, fit_gmlr$B[[index]], fit_gmlr$b[[index]]), 1, pred) != y_te)) / length(y_te)
    # SF[['GMLR_path']] <- which(rowSums(abs(fit_gmlr$B[[index]])) > 0)
    # grou_sel <- group[which(rowSums(abs(fit_gmlr$B[[index]])) != 0)]
    # SF_gr[['GMLR_path']] <- grou_sel[!duplicated(grou_sel)]
    # 
    # #GSMLR+path refined GAP GCD
    # fit_gprgg <- gsmlr_gapref(X, Y, group, dg, epsilon, maxiters, nlambda, display, lambda)
    # err_tu <- NULL
    # for (i in 1:nlambda) {
    #   prob_tu <- pb(x_tu, fit_gprgg$B[[i]], fit_gprgg$b[[i]])
    #   err_tu[i] <- length(which(apply(prob_tu, 1, pred) != y_tu))
    # }
    # index <- which.min(err_tu)
    # error[['GSMLR_gapref']] <- length(which(apply(pb(x_te, fit_gprgg$B[[index]], fit_gprgg$b[[index]]), 1, pred) != y_te)) / length(y_te)
    # SF[['GSMLR_gapref']] <- which(rowSums(abs(fit_gprgg$B[[index]])) > 0)
    # grou_sel <- group[which(rowSums(abs(fit_gprgg$B[[index]])) != 0)]
    # SF_gr[['GMLR_gapref']] <- grou_sel[!duplicated(grou_sel)]
    
    res <- list(error = error, SF = SF, SF_gr = SF_gr, lambda = lambda)
    return(res)
  }
  fit <- try(res <- main(1), silent=TRUE)
  if('try-error' %in% class(fit)){
    res <- "error"
  }
  return(res)
}
res_all <- list()
for (i in 65:100) {
  res_all[[i]] <- app_func(1)
  save(res_all, file = "res_smlr_ISOLET_lambda.Rdata")
  print(i)
}






cl <- detectCores()
clus <- makeCluster(cl - 2)
res_all <- list()
parm <- cbind(rep(1, 100))
res_all <- parLapply(clus, parm, fun = app_func)
save(res_all, file = "res_smlr_recog.Rdata")
# plot -------------------------------------------------------------------

#boxplot
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(ggpubr)
load("~/sim_GSMLR_1.Rdata")

time_all <- list()
for (i in 1:4) {
  res_time <- matrix(nc = 7, nr = 100)
  for (j in 1:100) {
    for (k in 7:13) {
      unit <- attr(res_all_sim_1[[i]][[j]][["time"]][[k]], "units")
      if(unit == "mins"){
        res_time[j, (k-6)] <- as.numeric(res_all_sim_1[[i]][[j]][["time"]][[k]]) * 60
      }else{
        res_time[j, (k-6)] <- as.numeric(res_all_sim_1[[i]][[j]][["time"]][[k]])
      }
    }
  }
  time_all[[i]] <- res_time
}
p_boxplot <- list()
for (i in 1:4) {
  data <- data.frame(
    name = rep(c("A", "B", "C", "D", "E", "F", "G"), each = 100),
    value = c(time_all[[i]])
  )
  xlab_name <- names(res_all_sim_1[[1]][[1]][["error"]])[c(7:13)]
  xlab_name[4:7] <- c("GMLR_path", "GSMLR_path", "GSMLR_GSS", "GSMLR_RSS")
  p_boxplot[[i]] <- 
    data %>%
    ggplot(aes(x = name, y = value, fill = name)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6, option ="viridis") +
    #theme_ipsum() +
    theme(panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_line(color = "gray", linewidth = 0.3),
          panel.grid.minor = element_blank())+
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      #axis.text.y = element_text(vjust = 0.5),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    ggtitle(paste("d=", i*50+50, sep = "")) +
    xlab("") + ylab("time(s)")+
    scale_x_discrete(labels = xlab_name)
}

ggarrange(p_boxplot[[1]], p_boxplot[[2]], 
             p_boxplot[[3]], p_boxplot[[4]], 
             nrow = 2, ncol = 2)






#radar
library(fmsb)
data <- as.data.frame(matrix(c(9298, 10, 9, 25, 256,
                               4480, 4, 15, 40, 529,
                               2600, 2, 9, 65, 500,
                               900, 15, 9, 60, 617), nr = 4, nc = 5, byrow = T))
rownames(data) <- c("USPS", "Recognition", "MADELON", "ISOLET")
colnames(data) <- c("Sample size", "Classes", "Groups", 
                    "Training data size", "Dimension")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(as.numeric(apply(data, 2, max)), as.numeric(apply(data, 2, min)), data)

# Color vector
colors_border = c(rgb(0.2, 0.5, 0.5, 0.9), rgb(0.8, 0.2, 0.5, 0.9), rgb(0.7, 0.5, 0.1, 0.9), rgb(0.3, 0.7, 0.2, 0.9))
colors_in = c(rgb(0.2, 0.5, 0.5, 0.4), rgb(0.8, 0.2, 0.5, 0.4), rgb(0.7, 0.5, 0.1, 0.4), rgb(0.3, 0.7, 0.2, 0.4))


# plot with default options:
radarchart(data, axistype = 1 , 
           #custom polygon
           pcol = colors_border, pfcol = colors_in, 
           plwd = 4 , plty=1,
           #custom the grid
           cglcol = "grey", cglty = 1, axislabcol = "grey", 
           cglwd = 0.8,
           #custom labels
           vlcex = 0.8 
)

# Add a legend
legend("topright", legend = c("USPS", "Recognition", "MADELON", "ISOLET"), 
       bty = "n", pch = 20 , col = colors_in , text.col = "black", 
       cex = 1.2, pt.cex = 3)



#line
load("~/res_smlr_ISOLET_lambda.Rdata")
res_lambda <- 0
for (i in 1:100) {
  res_lambda <- res_lambda + res_all[[i]][["lambda"]] / 100
}

res_error <- list()
for (i in 1:4) {
  res_error[[i]] <- 0
}
names(res_error) <- names(res_all[[1]][["error"]])

for (i in 1:100) {
  for (j in 1:4) {
    res_error[[j]] <- res_error[[j]] + res_all[[i]][["error"]][[j]] / 100
  }
}
error_seq <- data.frame(cbind(res_lambda, res_error[[1]], 
                                res_error[[3]], res_error[[4]]))
colnames(error_seq) <- c("k", "A", "B", "C")
error_seq <- reshape2::melt(data = error_seq, id.vars = "k")

error_plot <- ggplot(data = error_seq, 
                    aes(x = k, y = value, group = variable, 
                        color = variable, linetype = variable))+
  xlab(expression(lambda))+
  ylab("Error")+
  geom_line(size = 1.2)+
  #geom_point(data = error_seq, aes(x = k, y = value, colour = variable))+
  labs(color="", linetype = "")+
  scale_colour_manual(values = c("#000000", "#7CAE00",  "#00BFC4"), 
                      labels = c("LMLR", "GMLR", "GSMLR"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("LMLR", "GMLR", "GSMLR"))+
  theme(plot.title=element_text(hjust=0.5))+
  guides(color=guide_legend(override.aes = list(size=1.2)))+
  theme(legend.key.width=unit(1.2,'cm'))+
  theme(text = element_text(size = 15))+
  #theme_ipsum()+
  xlim(max(res_lambda), min(res_lambda))+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "gray", linewidth = 0.3),
        panel.grid.minor = element_blank())
error_plot


res_SF <- list()
for (i in 1:4) {
  res_SF[[i]] <- 0
}
names(res_SF) <- names(res_all[[1]][["SF"]])

for (i in 1:100) {
  for (j in 1:4) {
    res_SF[[j]] <- res_SF[[j]] + res_all[[i]][["SF"]][[j]] / 100
  }
}
SF_seq <- data.frame(res_lambda, res_SF[[1]], 
                              res_SF[[3]], res_SF[[4]])
colnames(SF_seq) <- c("k", "A", "B", "C")
SF_seq <- reshape2::melt(data = SF_seq, id.vars = "k")

SF_plot <- ggplot(data = SF_seq, 
                     aes(x = k, y = value, group = variable, 
                         color = variable, linetype = variable))+
  geom_line(size = 1.5)+
  #geom_point(data = SF_seq, aes(x = k, y = value, colour = variable))+
  labs(color="", linetype = "")+
  scale_colour_manual(values = c("#000000", "darkorange",  "#00BFC4"), 
                      labels = c("LMLR", "GMLR", "GSMLR"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("LMLR", "GMLR", "GSMLR"))+
  #labs(title = "Number of features selected by three methods")+
  theme(plot.title=element_text(hjust=0.5))+
  guides(color=guide_legend(override.aes = list(size=1.2)))+
  theme(legend.key.width=unit(1.2,'cm'))+
  theme(text = element_text(size = 15))+
  #theme_ipsum()+
  ylab("Number")+xlab(expression(lambda))+
  xlim(max(res_lambda), min(res_lambda))+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "gray", linewidth = 0.3),
        panel.grid.minor = element_line(color = "gray", linewidth = 0.3))
  # theme(panel.background = element_rect(fill = "transparent"),
  #       panel.grid.major = element_line(color = "gray", linewidth = 0.3),
  #       panel.grid.minor = element_blank())
SF_plot


res_SF_gr <- list()
for (i in 1:4) {
  res_SF_gr[[i]] <- 0
}
names(res_SF_gr) <- names(res_all[[1]][["SF_gr"]])

for (i in 1:100) {
  for (j in 1:4) {
    res_SF_gr[[j]] <- res_SF_gr[[j]] + res_all[[i]][["SF_gr"]][[j]] / 100
  }
}
SF_gr_seq <- data.frame(res_lambda, res_SF_gr[[1]], res_SF_gr[[3]], res_SF_gr[[4]])
colnames(SF_gr_seq) <- c("k", "A", "B", "C")
SF_gr_seq <- reshape2::melt(data = SF_gr_seq, id.vars = "k")

SF_gr_plot <- ggplot(data = SF_gr_seq, 
                  aes(x = k, y = value, group = variable, 
                      color = variable, linetype = variable))+
  xlab(expression(n))+
  ylab("Value")+
  geom_line(size = 1.5)+
  #geom_point(data = SF_gr_seq, aes(x = k, y = value, colour = variable))+
  labs(color="", linetype = "")+
  scale_colour_manual(values = c("#000000", "darkorange",  "#00BFC4"), 
                      labels = c("LMLR", "GMLR", "GSMLR"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("LMLR", "GMLR", "GSMLR"))+
  #ggtitle("Number of groups selsted by three methods") +
  theme(plot.title = element_text(hjust=0.5),
        legend.key.width=unit(1.2,'cm'))+
  guides(color=guide_legend(override.aes = list(size=1.2)))+
  theme(text = element_text(size = 15))+
  #theme_ipsum()+
  ylab("Number")+xlab(expression(lambda))+
  xlim(max(res_lambda), min(res_lambda))+
  theme(panel.background = element_rect(fill = "transparent"),
                                               panel.grid.major = element_line(color = "gray", linewidth = 0.3),
                                               panel.grid.minor = element_line(color = "gray", linewidth = 0.3))
  # theme(panel.background = element_rect(fill = "transparent"),
  #       panel.grid.major = element_line(color = "gray", linewidth = 0.3),
  #       panel.grid.minor = element_blank())
SF_gr_plot

ggarrange(SF_plot, 
          SF_gr_plot, 
          nrow = 1, ncol = 2, common.legend = T)
