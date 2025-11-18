# Packages ----------------------------------------------------------------
library(glmnet)
library(mvtnorm)
library(nnet)
library(e1071)
library(grpnet)
library(dplyr)
library(randomForest)
library(party)
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

#data need to be ontained first
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
            cv.glmnet(x_tr, y_tr)$res$lambda.1se)
lambda.min.ratio <- ifelse(n < d, 0.01, 1e-04)
lambda <- lam0*seq(1, lambda.min.ratio, length.out=nlambda)


time <- list()
error <- list()
SF <- list()

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

#RF
start_time <- Sys.time()
fit_rf <- randomForest(x_tr, as.factor(y_tr))
end_time <- Sys.time()
pre_res <- as.numeric(predict(fit_rf, x_te))
time[['RF']] <- end_time - start_time
error[['RF']] <- length(which(pre_res != y_te)) / length(y_te)
SF[['RF']] <- d

#DT
start_time <- Sys.time()
data_dt <- data.frame(x_tr, as.factor(y_tr))
colnames(data_dt)[ncol(data_dt)] <- "y_tr"
fit_dt <- ctree(y_tr~.,  data_dt)
end_time <- Sys.time()
pre_res <- as.numeric(predict(fit_dt, data.frame(x_te)))
time[['DT']] <- end_time - start_time
error[['DT']] <- length(which(pre_res != y_te)) / length(y_te)
SF[['DT']] <- d

#nB
start_time <- Sys.time()
fit_nb <- naiveBayes(x_tr, as.factor(y_tr))
end_time <- Sys.time()
pre_res <- predict(object = fit_nb, newdata = x_te)
time[['nB']] <- end_time - start_time
error[['nB']] <- length(which(pre_res != y_te)) / length(y_te)
SF[['nB']] <- d

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
SF[['LMLR']] <- as.numeric(which(beta_lmlr > 0.002))


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
SF[['RMLR']] <- as.numeric(which(beta_rmlr > 0.002))

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
