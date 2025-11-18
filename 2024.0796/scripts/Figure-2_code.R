# Working path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Packages ----------------------------------------------------------------
library(parallel)
library(glmnet)
# library(msgl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
# Application study for the ISOLET database -------------------------------
app_func <- function(para){
  load("isolet.Rdata")
  library(glmnet)
  # library(msgl)
  library(dplyr)
  # BasiC functionS ----------------------------------------------------------
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
  # Data&setting --------------------------------------------------------------------
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
  # Main function--------------------------------------------------------------------
  main <- function(para){
    #lambda-sequence
    lam0 <- max(lammaxlasso(X, Y),
                lammax(X, Y, k, group, dg),
                cv.glmnet(x_tr, y_tr)$res$lambda.1se)
    lambda.min.ratio <- 0.1
    lambda <- lam0 * seq(2, lambda.min.ratio, length.out=nlambda)
    
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
    
    #GMLR+path GCD
    fit_gmlr <- gmlr_path(X,Y,group,dg,epsilon,maxiters,nlambda,display,lambda)
    err_gmlr <- NULL
    for (i in 1:nlambda) {
      prob_te <- pb_gmlr(x_te, fit_gmlr$B[[i]], fit_gmlr$b[[i]])
      err_gmlr[i] <- length(which(apply(prob_te, 1, pred) != y_te)) / length(y_te)
      SF[['GMLR']][i] <- length(which(rowSums(abs(fit_gmlr$B[[i]])) != 0))
      grou_gmlr <- group[which(rowSums(abs(fit_gmlr$B[[i]])) != 0)]
      SF_gr[['GMLR']][i] <- length(grou_gmlr[!duplicated(grou_gmlr)])
    }
    error[['GMLR']] <- err_gmlr
    
    
    
    #GSMLR
    fit_gprgg <- gsmlr_gapref(X, Y, group, dg, epsilon, maxiters, nlambda, display, lambda)
    err_sgl <- NULL
    for (i in 1:nlambda) {
      prob_te <- pb(x_te, fit_gprgg$B[[i]], fit_gprgg$b[[i]])
      err_sgl[i] <- length(which(apply(prob_te, 1, pred) != y_te)) / length(y_te)
      SF[['GSMLR']][i] <- length(which(rowSums(abs(fit_gprgg$B[[i]])) != 0))
      grou_sel <- group[which(rowSums(abs(fit_gprgg$B[[i]])) != 0)]
      SF_gr[['GSMLR']][i] <- length(grou_sel[!duplicated(grou_sel)])
    }
    error[['GSMLR']] <- err_sgl
    
    res <- list(error = error, SF = SF, SF_gr = SF_gr, lambda = lambda)
    return(res)
  }
  fit <- try(res <- main(1), silent=TRUE)
  if('try-error' %in% class(fit)){
    res <- "error"
  }
  return(res)
}
cl <- detectCores()
clus <- makeCluster(cl - 2)
res_all <- list()
parm <- cbind(rep(1, 100))
res_all <- parLapply(clus, parm, fun = app_func)
save(res_all, file = "res_smlr_ISOLET_lambda.Rdata")
# Plot --------------------------------------------------------------------
load("res_smlr_ISOLET_lambda.Rdata")
res_lambda <- 0
for (i in 1:100) {
  res_lambda <- res_lambda + res_all[[i]][["lambda"]] / 100
}

res_SF <- list()
for (i in 1:3) {
  res_SF[[i]] <- 0
}
names(res_SF) <- names(res_all[[1]][["SF"]])

for (i in 1:100) {
  for (j in 1:3) {
    res_SF[[j]] <- res_SF[[j]] + res_all[[i]][["SF"]][[j]] / 100
  }
}
SF_seq <- data.frame(res_lambda, res_SF[[1]], res_SF[[2]],
                     res_SF[[3]])
colnames(SF_seq) <- c("k", "A", "B", "C")
SF_seq <- melt(data = SF_seq, id.vars = "k")

SF_plot <- ggplot(data = SF_seq, 
                  aes(x = k, y = value, group = variable, 
                      color = variable))+
  geom_line(size = 1.5)+
  labs(color="", linetype = "")+
  scale_colour_manual(values = c("#000000", "red",  "blue"), 
                      labels = c("LMLR", "GMLR", "GSMLR"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("LMLR", "GMLR", "GSMLR"))+
  theme(plot.title=element_text(hjust=0.5))+
  guides(color=guide_legend(override.aes = list(size=1.2)))+
  theme(legend.key.width=unit(1.2,'cm'))+
  theme(text = element_text(size = 15))+
  ylab("Number")+xlab(expression(lambda))+
  xlim(max(res_lambda), min(res_lambda))+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "gray", linewidth = 0.3),
        panel.grid.minor = element_line(color = "gray", linewidth = 0.3))


res_SF_gr <- list()
for (i in 1:3) {
  res_SF_gr[[i]] <- 0
}
names(res_SF_gr) <- names(res_all[[1]][["SF_gr"]])

for (i in 1:100) {
  for (j in 1:3) {
    res_SF_gr[[j]] <- res_SF_gr[[j]] + res_all[[i]][["SF_gr"]][[j]] / 100
  }
}
SF_gr_seq <- data.frame(res_lambda, res_SF_gr[[1]], res_SF_gr[[2]], res_SF_gr[[3]])
colnames(SF_gr_seq) <- c("k", "A", "B", "C")
SF_gr_seq <- melt(data = SF_gr_seq, id.vars = "k")

SF_gr_plot <- ggplot(data = SF_gr_seq, 
                     aes(x = k, y = value, group = variable, 
                         color = variable))+
  xlab(expression(n))+
  ylab("Value")+
  geom_line(size = 1.5)+
  labs(color="", linetype = "")+
  scale_colour_manual(values = c("#000000", "red",  "blue"), 
                      labels = c("LMLR", "GMLR", "GSMLR"))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), 
                        labels = c("LMLR", "GMLR", "GSMLR"))+
  theme(plot.title = element_text(hjust=0.5),
        legend.key.width=unit(1.2,'cm'))+
  guides(color=guide_legend(override.aes = list(size=1.2)))+
  theme(text = element_text(size = 15))+
  ylab("Number")+xlab(expression(lambda))+
  xlim(max(res_lambda), min(res_lambda))+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "gray", linewidth = 0.3),
        panel.grid.minor = element_line(color = "gray", linewidth = 0.3))

ggarrange(SF_plot, 
          SF_gr_plot, 
          nrow = 1, ncol = 2, common.legend = T)


