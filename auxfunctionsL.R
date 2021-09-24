##############################################################################
# EM - main functions
##############################################################################

#general case
EM.Skew<- function(formFixed, #two sided formula for the fixed effects (e.g. y ~ x)
                   formRandom, #one sided formula for the random effects (e.g. ~ x)
                   data, #object containing the dataset
                   groupVar, #a string containing the name of the grouping variable (in data)
                   distr,#'st' or 'sn' 
                   beta1,sigmae,D1,Deltab,nu=NULL, #initial values
                   lb=2,lu=100, #limit values for nu (if distr=='st')
                   tol=1e-5, #tolerance for EM algorithm
                   estimanuint=TRUE, #should nu be estimated as an integer?
                   estimanuoptim=FALSE, #should nu be estimated as an real number? If TRUE, time required may increase considerably
                   informa=TRUE, #should SE (information matrix) be computed?
                   max.iter=500, #maximum number of iterations
                   showiter=TRUE, #if FALSE nothing will be printed
                   showerroriter=FALSE #if TRUE the current convergence criteria will be printed at each iteration
){
  ti = Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  r<- ncol(Deltab)
  #
  if (all(Deltab==0)) {
    logveroDelta<-function(deltav){
      Deltab <- matrix(deltav,nrow=nrow(D1))
      -logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)
    }
    Deltab <- matrix(NMOF::gridSearch(logveroDelta,lower=rep(-1,length(Deltab)),
                                      upper=rep(1,length(Deltab)),
                                      n=4,printDetail = FALSE)$minlevel,
                     nrow=nrow(D1))
    if (showiter) cat('Initial Delta: ',Deltab,'\n')
  }
  #
  if (distr!='sn'&&estimanuint) if (lb-round(lb)>1e-5) lb = max(2,round(lb))  
  if (distr!='sn'&&estimanuint) if (lu-round(lu)>1e-5) lu = round(lu)
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Deltab,nu)
  
  criterio<-10
  count<-0
  llji = logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")
  
  while((criterio > tol)&&(count<max.iter)){
    
    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emj,y=y, x=x, z=z, beta1=beta1, D1=D1,
                                 Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    sum6 = Reduce("+",res_emj$sum6)
    uhat = unlist(res_emj$uhat,use.names = F)
    
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    D1<-sum4/m
    Deltab<-sum5%*%solve(sum6)
    #
    if (distr!="sn" && estimanuoptim) {
      logvero1<-function(nu){logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)}
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    else if (distr!="sn" && estimanuint) {
      logvero2<-function(nu){-logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)}
      nu <- NMOF::gridSearch(logvero2,levels=list(lb:lu),printDetail = FALSE)$minlevel
    }
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Deltab,nu)
    llj1 <- llji
    llji <- logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)
    if (count>10) criterio <- abs((llji-llj1)/llj1)
    if (showiter&&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  dd<-D1[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,Deltab,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (is.null(colnames(z))) colnames(z) <- paste0("b",0:(ncol(D1)-1))
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("D",1:length(dd)),paste0("Delta",rep(1:q1,r),rep(1:r,each=q1)))
  else names(theta)<- c(colnames(x),"sigma2",paste0("D",1:length(dd)),paste0("Delta",rep(1:q1,r),rep(1:r,each=q1)),'nu')
  res_bi <- tapply(1:N,ind,calcbi,y=y, x=x, z=z, beta1=beta1, D1=D1,
                   Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu)
  res_bi_out <- matrix(unlist(res_bi),ncol=2,byrow = T)
  rownames(res_bi_out) <- names(res_bi)
  colnames(res_bi_out) <- colnames(z) 
  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           dsqrt=dd,D=D1,Delta=Deltab),
                  random.effects=res_bi_out,uhat=uhat)
  if (distr != "sn") obj.out$estimates$nu = nu
  tf = Sys.time()
  if (informa) {
    #Louis
    res_score = tapply(1:N,ind,scorej,y=y, x=x, z=z, beta1=beta1, D1=D1,
                       Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu)
    h2 <- try(solve(Reduce("+",res_score)),silent = T)
    if (class(h2)[1]=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error2=rep(NA,length(theta))
    } else{
      desvios2 <- try(c(diag(matrix.sqrt(h2)),rep(NA,length(nu))),silent = T)
      if (class(desvios2)=="try-error") desvios2=rep(NA,length(theta))
      names(desvios2) <- names(theta)
      obj.out$std.error=desvios2
    }
  }
  obj.out$loglik <-as.numeric(llji)
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

#special case of Sahu, Branco & Dey (2003)
EM.Skew.Sahu<- function(formFixed, #two sided formula for the fixed effects (e.g. y ~ x)
                        formRandom, #one sided formula for the random effects (e.g. ~ x)
                        data, #object containing the dataset
                        groupVar, #a string containing the name of the grouping variable (in data)
                        distr,#'st' or 'sn' 
                        beta1,sigmae,D1,Deltav,nu=NULL, #initial values, Deltav=diag(Deltab) should be a vector
                        lb=2,lu=100, #limit values for nu (if distr=='st')
                        tol=1e-5, #tolerance for EM algorithm
                        estimanuint=TRUE, #should nu be estimated as an integer?
                        estimanuoptim=FALSE, #should nu be estimated as an real number? If TRUE, time required may increase considerably
                        informa=TRUE, #should SE (information matrix) be computed?
                        max.iter=500, #maximum number of iterations
                        showiter=TRUE, #if FALSE nothing will be printed
                        showerroriter=FALSE #if TRUE the current convergence criteria will be printed at each iteration
){
  ti = Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #r<- ncol(Deltab)
  if (length(Deltav)==1) stop('For scalar Delta please use the function EM.Skew')
  Deltab <- diag(Deltav)
  #
  if (all(Deltab==0)) {
    logveroDeltav<-function(deltav){
      -logvero(y, x, z, ind, beta1, sigmae, D1, diag(deltav), distr, nu)
    }
    Deltav <- matrix(NMOF::gridSearch(logveroDeltav,lower=rep(-1,length(Deltav)),
                                      upper=rep(1,length(Deltav)),
                                      n=4,printDetail = FALSE)$minlevel,
                     nrow=nrow(D1)) %>% as.numeric()
    Deltab <- diag(Deltav)
    if (showiter) cat('Initial Delta: ',Deltab,'\n')
  }
  #
  if (distr!='sn'&&estimanuint) if (lb-round(lb)>1e-5) lb = max(2,round(lb))  
  if (distr!='sn'&&estimanuint) if (lu-round(lu)>1e-5) lu = round(lu)
  
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Deltav,nu)
  
  criterio<-10
  count<-0
  llji = logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")
  
  while((criterio > tol)&&(count<max.iter)){
    
    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emj,y=y, x=x, z=z, beta1=beta1, D1=D1,
                                 Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    sum6 = Reduce("+",res_emj$sum6)
    uhat = unlist(res_emj$uhat,use.names = F)
    
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    D1<-sum4/m
    Deltav<-diag(sum5%*%solve(sum6))
    Deltab<- diag(Deltav)
    #
    if (distr!="sn" && estimanuoptim) {
      logvero1<-function(nu){logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)}
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    else if (distr!="sn" && estimanuint) {
      logvero2<-function(nu){-logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)}
      nu <- NMOF::gridSearch(logvero2,levels=list(lb:lu),printDetail = FALSE)$minlevel
    }
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Deltav,nu)
    llj1 <- llji
    llji <- logvero(y, x, z, ind, beta1, sigmae, D1, Deltab, distr, nu)
    if (count>10) criterio <- abs((llji-llj1)/llj1)
    if (showiter&&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") 
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") 
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  dd<-D1[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,Deltav,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (is.null(colnames(z))) colnames(z) <- paste0("b",0:(ncol(D1)-1))
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("D",1:length(dd)),
                                   paste0("Delta",1:q1,1:q1))
  else names(theta)<- c(colnames(x),"sigma2",paste0("D",1:length(dd)),
                        paste0("Delta",1:q1,1:q1),'nu')
  res_bi <- tapply(1:N,ind,calcbi,y=y, x=x, z=z, beta1=beta1, D1=D1,
                   Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu)
  res_bi_out <- matrix(unlist(res_bi),ncol=2,byrow = T)
  rownames(res_bi_out) <- names(res_bi)
  colnames(res_bi_out) <- colnames(z) 
  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           dsqrt=dd,D=D1,Delta=Deltab),
                  random.effects=res_bi_out,uhat=uhat)
  if (distr != "sn") obj.out$estimates$nu = nu
  tf = Sys.time()
  if (informa) {
    #Louis
    res_score = tapply(1:N,ind,scorej,y=y, x=x, z=z, beta1=beta1, D1=D1,
                       Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu)
    h1 <- Reduce("+",res_score)
    indh1<-c(colnames(x),"sigma2",paste0("D",1:length(dd)),paste0("Delta",rep(1:q1,q1),rep(1:q1,each=q1))) %in% names(theta)
    h2 <- try(solve(h1[indh1,indh1]),silent = T)
    if (class(h2)[1]=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error2=rep(NA,length(theta))
    } else{
      desvios2 <- try(c(diag(matrix.sqrt(h2)),rep(NA,length(nu))),silent = T)
      if (class(desvios2)=="try-error") desvios2=rep(NA,length(theta))
      names(desvios2) <- names(theta)
      obj.out$std.error=desvios2
    }
  }
  obj.out$loglik <-as.numeric(llji)
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

################################################################
# creating some auxiliar functions
################################################################

is.wholenumber <- function(x, tol1 = .Machine$double.eps^0.5)  abs(x - round(x)) < tol1

#Root of a symmetric matrix
matrix.sqrt <- function(A)
{
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd e decomposi??o espectral
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}

#trace of a matrix of dim >=1
traceM <- function(Mat){
  if(length(Mat)==1) tr<- as.numeric(Mat)
  else tr<-sum(diag(Mat))
  tr
}

#revert a list
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}

#
Dmatrix <- function(dd) {
  q2 <- length(dd)
  q1 <- -.5+sqrt(1+8*q2)/2
  if (q1%%1 != 0) stop("wrong dimension of dd")
  D1 <- matrix(nrow = q1,ncol = q1)
  D1[upper.tri(D1,diag = T)] <- as.numeric(dd)
  D1[lower.tri(D1)] <- t(D1)[lower.tri(D1)]
  return(D1)
}

# generating a ST sample (for 1 subject)
gen_ind_ST = function(ni,Sig,Di,beta,Delta,#q x r matrix
                        distr="sn",#"sn" or "st"
                        nu=NULL) {
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  Xi = cbind(1,((1:ni)-(1/2+ni/2))/(ni-1)*2) #sim2
  Zi =Xi
  Beta = matrix(beta,ncol=1)
  q1 = nrow(Di)
  r = ncol(Delta)
  si = t(abs(rmvnorm(1,sigma=ui^(-1)*diag(r))))
  bi = t(rmvnorm(1,Delta%*%(c.*matrix(1,ncol=1,nrow=r)+si),sigma=ui^(-1)*Di))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))
  return(data.frame(y=Yi,x=Xi[,2],tempo=1:ni,ui=ui,b0=bi[1],b1=bi[2]))
}

################################################################
#Log-likelihood 
################################################################
ljnormal <-function(j,y,x,z,beta1,D1,Deltab,sigmae){
  c. = -sqrt(2/pi)
  y1=y[j]
  p= ncol(x);q1=ncol(z);r=ncol(Deltab)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab%*%matrix(1,ncol=1,nrow=r)
  njj = length(y1)
  Psi<-(z1)%*%D1%*%t(z1)+sigmae*diag(njj) 
  Sigb <- D1 + Deltab%*%t(Deltab)
  Sigy <- Psi + z1%*%Deltab%*%t(z1%*%Deltab)
  sSigy <- solve(Sigy)
  Lambday <- diag(r) - t(z1%*%Deltab)%*%sSigy%*%z1%*%Deltab
  Lambday <- (Lambday+t(Lambday))/2
  Ajj<-t(z1%*%Deltab)%*%sSigy%*%(y1-med)
  log(2^r*dmvnorm(y1,med,Sigy)*as.numeric(pmvnorm(upper=as.numeric(Ajj),sigma=Lambday)))
}
#
ljt <-function(j,nu,y,x,z,beta1,D1,Deltab,sigmae){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  p= ncol(x);q1=ncol(z);r=ncol(Deltab)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab%*%matrix(1,ncol=1,nrow=r)
  njj = length(y1)
  Psi<-(z1)%*%D1%*%t(z1)+sigmae*diag(njj) 
  Sigb <- D1 + Deltab%*%t(Deltab)
  Sigy <- Psi + z1%*%Deltab%*%t(z1%*%Deltab)
  sSigy <- solve(Sigy)
  Lambday <- diag(r) - t(z1%*%Deltab)%*%sSigy%*%z1%*%Deltab
  Lambday <- (Lambday+t(Lambday))/2
  Ajj<-t(z1%*%Deltab)%*%sSigy%*%(y1-med)*as.numeric(sqrt((nu+njj)/(nu+t(y1-med)%*%sSigy%*%(y1-med))))
  log(2^r*dmvt(y1,delta = med, sigma = Sigy, df = nu,log=F)*
        as.numeric(pmvnormt(upper=as.numeric(Ajj),sigma=Lambday,nu = nu+njj)))
}

logvero = function(y,x,z,ind,beta1,sigmae,D1,Deltab,distr,nu){ #ind = indicadora de individuo
  N<-length(ind)
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormal,y=y,x=x,z=z,beta1=beta1,D1=D1,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljt,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,Deltab=Deltab,sigmae=sigmae))
  lv
}

##############################################################################
# EM auxiliar functions
##############################################################################

calcbi = function(jseq, y, x, z, beta1, D1, Deltab, sigmae,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  #
  y1=y[jseq]
  p= ncol(x);q1=ncol(z);r=ncol(Deltab)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab%*%matrix(1,ncol=1,nrow=r)
  nj = length(y1)
  Psi<-(z1)%*%D1%*%t(z1)+sigmae*diag(nj)
  Sigb <- D1 + Deltab%*%t(Deltab)
  Sigy <- Psi + z1%*%Deltab%*%t(z1%*%Deltab)
  sSigy <- solve(Sigy)
  dj<-as.numeric(t(y1-med)%*%sSigy%*%(y1-med))
  Lambday <- diag(r) - t(z1%*%Deltab)%*%sSigy%*%z1%*%Deltab
  Lambday <- (Lambday+t(Lambday))/2
  Ajj<-t(z1%*%Deltab)%*%sSigy%*%(y1-med)
  sD1 <- solve(D1)
  #
  if  (distr=="sn"){
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),Sigma = Lambday,dist = 'normal')
    shat <- ttrunc$mean
  }
  
  if (distr=="st"){
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),
                         Sigma = (nu+dj)/(nu+nj)*Lambday,dist = 't',nu = nu+nj)
    shat <- ttrunc$mean
  }
  
  Bmat <- solve(sD1+t(z1)%*%z1/sigmae) #Tbj<-solve(solve(Gammab)+t(z1)%*%z1/sigmae)
  temp1<-c.*Deltab%*%matrix(1,nrow=r)+Bmat%*%t(z1)%*%(y1-med)/sigmae
  ubhat <- temp1+Bmat%*%sD1%*%Deltab%*%shat
  return(ubhat)
}


emj = function(jseq, y, x, z, beta1, D1, Deltab, sigmae,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  #
  y1=y[jseq]
  p= ncol(x);q1=ncol(z);r=ncol(Deltab)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab%*%matrix(1,ncol=1,nrow=r)
  nj = length(y1)
  Psi<-(z1)%*%D1%*%t(z1)+sigmae*diag(nj)
  Sigb <- D1 + Deltab%*%t(Deltab)
  Sigy <- Psi + z1%*%Deltab%*%t(z1%*%Deltab)
  sSigy <- solve(Sigy)
  dj<-as.numeric(t(y1-med)%*%sSigy%*%(y1-med))
  Lambday <- diag(r) - t(z1%*%Deltab)%*%sSigy%*%z1%*%Deltab
  Lambday <- (Lambday+t(Lambday))/2
  Ajj<-t(z1%*%Deltab)%*%sSigy%*%(y1-med)
  sD1 <- solve(D1)
  #
  if  (distr=="sn"){
    uhat<-1
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),Sigma = Lambday,dist = 'normal')
    ushat <- ttrunc$mean*uhat
    us2hat <- ttrunc$EYY*uhat
  }
  if (distr=="st"){
    uhat <- (nu+nj)/(nu+dj)*
      as.numeric(pmvnormt(upper=as.numeric(Ajj)*sqrt((nu+nj+2)/(nu+dj)),sigma=Lambday,nu = nu+nj+2))/
      as.numeric(pmvnormt(upper=as.numeric(Ajj)*sqrt((nu+nj)/(nu+dj)),sigma=Lambday,nu = nu+nj))
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),Sigma = (nu+dj)/(nu+nj+2)*Lambday,dist = 't',nu = nu+nj+2)
    ushat <- ttrunc$mean*uhat
    us2hat <- ttrunc$EYY*uhat
  }
  
  Bmat <- solve(sD1+t(z1)%*%z1/sigmae) 
  temp1<-c.*Deltab%*%matrix(1,nrow=r)+Bmat%*%t(z1)%*%(y1-med)/sigmae
  ubhat <- temp1*uhat+Bmat%*%sD1%*%Deltab%*%ushat
  ubshat <- temp1%*%t(ushat) + Bmat%*%sD1%*%Deltab%*%us2hat
  ub2hat <- Bmat+ubhat%*%t(temp1)+ubshat%*%t(Deltab)%*%sD1%*%Bmat
  ub2hat <- (ub2hat+t(ub2hat))/2
  #
  #M step
  sum1<-uhat*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uhat*y1-z1%*%ubhat)) #num beta
  sum3<-uhat*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ubhat-
    t(ubhat)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2hat%*%t(z1)%*%z1) #sum of sig2
  sum4<-ub2hat + Deltab%*%us2hat%*%t(Deltab)+ uhat*c.^2*Deltab%*%matrix(1,nrow=r,ncol=r)%*%t(Deltab)- 
    c.*(Deltab%*%matrix(1,nrow=r)%*%t(ubhat) + ubhat%*%matrix(1,ncol=r)%*%t(Deltab))+
    c.*(Deltab%*%matrix(1,nrow=r)%*%t(ushat)%*%t(Deltab)+ Deltab%*%ushat%*%matrix(1,ncol=r)%*%t(Deltab))-
    Deltab%*%t(ubshat)-ubshat%*%t(Deltab)#sum of D1
  sum5<- c.*ubhat%*%matrix(1,ncol=r)+ubshat #sum of Delta
  sum6<- us2hat+c.^2*uhat*matrix(1,ncol=r,nrow=r)+c.*(matrix(1,nrow=r)%*%t(ushat)+ushat%*%matrix(1,ncol=r)) #den do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,sum6=sum6,uhat=uhat)
  return(obj.out)
}

##############################################################################
# functions for standard errors (SE)
##############################################################################

logveroVec = function(theta,y,x.,z,ind,distr,r,nu=NULL){ #ind = indicator of subject
  ##theta = c(beta1,sigmae,dd,Deltab,nu)
  N<-length(ind)
  p<-dim(x.)[2]
  q1<-dim(z)[2]
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  dd <- theta[(p+2):(p+1+q2)]
  D1 <- Dmatrix(dd)
  Deltab <-matrix(theta[(p+2+q2):(p+1+q2+r*q1)],ncol=r)
  #
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormal,y=y,x=x.,z=z,beta1=beta1,D1=D1,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljt,nu=nu,y=y,x=x.,z=z,beta1=beta1,D1=D1,Deltab=Deltab,sigmae=sigmae))
  lv
}
scorej = function(jseq, y, x, z, beta1, D1, Deltab, sigmae,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  #
  y1=y[jseq]
  p= ncol(x);q1=ncol(z);r=ncol(Deltab)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab%*%matrix(1,ncol=1,nrow=r)
  nj = length(y1)
  Psi<-(z1)%*%D1%*%t(z1)+sigmae*diag(nj)
  Sigb <- D1 + Deltab%*%t(Deltab)
  Sigy <- Psi + z1%*%Deltab%*%t(z1%*%Deltab)
  sSigy <- solve(Sigy)
  dj<-as.numeric(t(y1-med)%*%sSigy%*%(y1-med))
  Lambday <- diag(r) - t(z1%*%Deltab)%*%sSigy%*%z1%*%Deltab
  Lambday <- (Lambday+t(Lambday))/2
  Ajj<-t(z1%*%Deltab)%*%sSigy%*%(y1-med)
  sD1 <- solve(D1)
  #
  if  (distr=="sn"){
    uhat<-1
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),Sigma = Lambday,dist = 'normal')
    ushat <- ttrunc$mean*uhat
    us2hat <- ttrunc$EYY*uhat
  }
  if (distr=="st"){
    uhat <- (nu+nj)/(nu+dj)*
      as.numeric(pmvnormt(upper=as.numeric(Ajj)*sqrt((nu+nj+2)/(nu+dj)),sigma=Lambday,nu = nu+nj+2))/
      as.numeric(pmvnormt(upper=as.numeric(Ajj)*sqrt((nu+nj)/(nu+dj)),sigma=Lambday,nu = nu+nj))
    ttrunc <- meanvarTMD(lower=rep(0,r),mu = as.numeric(Ajj),Sigma = (nu+dj)/(nu+nj+2)*Lambday,dist = 't',nu = nu+nj+2)
    ushat <- ttrunc$mean*uhat
    us2hat <- ttrunc$EYY*uhat
  }
  
  Bmat <- solve(sD1+t(z1)%*%z1/sigmae)
  temp1<-c.*Deltab%*%matrix(1,nrow=r)+Bmat%*%t(z1)%*%(y1-med)/sigmae
  ubhat <- temp1*uhat+Bmat%*%sD1%*%Deltab%*%ushat
  ubshat <- temp1%*%t(ushat) + Bmat%*%sD1%*%Deltab%*%us2hat
  ub2hat <- Bmat+ubhat%*%t(temp1)+ubshat%*%t(Deltab)%*%sD1%*%Bmat
  ub2hat <- (ub2hat+t(ub2hat))/2
  #
  #score
  scbeta<- 1/sigmae*t(x1)%*%(uhat*(y1-x1%*%beta1)-z1%*%ubhat)
  scsigmae<- -nj/(2*sigmae)+1/(2*sigmae^2)*(uhat*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ubhat-
                                              t(ubhat)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2hat%*%t(z1)%*%z1)) 
  scD1<- -1/2*sD1+1/2*sD1%*%(ub2hat + Deltab%*%us2hat%*%t(Deltab)+ uhat*c.^2*Deltab%*%matrix(1,nrow=r,ncol=r)%*%t(Deltab)- 
                               c.*(Deltab%*%matrix(1,nrow=r)%*%t(ubhat) + ubhat%*%matrix(1,ncol=r)%*%t(Deltab))+
                               c.*(Deltab%*%matrix(1,nrow=r)%*%t(ushat)%*%t(Deltab)+ Deltab%*%ushat%*%matrix(1,ncol=r)%*%t(Deltab))-
                               Deltab%*%t(ubshat)-ubshat%*%t(Deltab))%*%sD1
  scDelta<- sD1%*%(c.*ubhat%*%matrix(1,ncol=r)+ubshat) -sD1%*%Deltab%*%(
    us2hat+c.^2*uhat*matrix(1,ncol=r,nrow=r)+c.*(matrix(1,nrow=r)%*%t(ushat)+ushat%*%matrix(1,ncol=r)))
  ##teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Deltab,nu)
  scvec <- c(scbeta,scsigmae,scD1[upper.tri(D1, diag = T)],scDelta)
  return(matrix(scvec)%*%t(matrix(scvec)))
}
