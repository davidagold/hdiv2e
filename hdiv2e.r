##two exogenous variables###

library(purrr)
library(mvtnorm)
library(glmnet)
library(dplyr)
source("estimation.r")
source("utils.r")

##########config###############

trial <- function() {
  args = commandArgs(trailingOnly=TRUE)
  config_id_ <- args[1] %>% as.numeric
  trial_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  configs <- read.csv("configs.csv")
  config <- filter(configs, config_id == config_id_)
  T <- config$T; J <- config$J; L <- config$L
  set.seed(trial_id)
  nn=1000
  TJ=T*J

  sigmaxi=0.5
  xi=rnorm(TJ)*sigmaxi
  betatrue=c(-1,1) # true mean tastes
  K=length(betatrue) # # of attributes
  covX=matrix(c(1,0.5,0.5,1),K,K)
  x=matrix(rnorm(TJ*K),TJ,K)%*%(covX) # products attributes TJ by dim(beta)
  delta=x%*%betatrue+xi #TJ by 1

  v=matrix(rnorm(nn),K,nn) #K by nn
  covrc=diag(.5,K,K)
  rc=t(chol(covrc))%*%v # K by nn
  MU=x%*%rc #TJ by nn

  numer=exp(kronecker(matrix(1,1,nn),delta)+MU) #numerator is a TJ by nn matrix
  tem=cbind(seq(TJ),kronecker(seq(J),matrix(1,T,1)),kronecker(matrix(1,J,1),seq(T))) #TJ by dim(beta)
  datasort=tem[order(tem[,3],tem[,2])]
  a=numer[datasort,]
  temp=apply(a,2,cumsum) #TJ by nn
  temp1=temp[c((1:T)*J),] #T by nn
  temp1[c(2:nrow(temp1)),]=diff(temp1)

  denom=1+temp1 #T by nn
  temp2=kronecker(matrix(1,J,1),denom)
  s=rowMeans(numer/temp2) #T*J by 1: averaging over nn
  s0=rowMeans(1/denom) #T by 1: averaging over nn

  sum(s)+sum(s0)
  min(s)
  max(s)

  R=log(s)-log(kronecker(matrix(1,J,1),s0))
  sMatrix=matrix(s,nrow=T,ncol=J)
  x1Mat<-t(matrix(x[,1],nrow=T,ncol=J))
  x2Mat<-t(matrix(x[,2],nrow=T,ncol=J))
  share<-t(sMatrix) #J by T


  dd1=array(NA,c(J,T,J))
  dz1=dd1
  a1=matrix(NA,J,T)
  D1=matrix(0,(J*J*T),1)
  Dz1=D1

  dd2=array(NA,c(J,T,J))
  dz2=dd2
  a2=matrix(NA,J,T)
  D2=matrix(0,(J*J*T),1)
  Dz2=D2


  ss=array(NA,c(J,T,J))
  sss=matrix(NA,J,T)
  S=matrix(0,(J*J*T),1)

  S0matrix=t(kronecker(matrix(1,1,J),s0))

  ss0=array(NA,c(J,T,J))
  sss0=matrix(NA,J,T)
  S0=matrix(0,(J*J*T),1)

  for (j in 1:J){
    for (t in 1:T){
      sss[,t]=share[j,t]
      ss[,,j]=sss
      S[((j-1)*TJ+1):(j*(TJ))]=as.vector(ss[,,j])

      sss0[,t]=S0matrix[j,t]
      ss0[,,j]=sss0
      S0[((j-1)*TJ+1):(j*(TJ))]=as.vector(ss[,,j])

      a1[,t]=(x1Mat[j,t]-x1Mat[,t])/share[j,t]
      dd1[,,j]=a1
      D1[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd1[,,j])

      a1[,t]=(x1Mat[j,t]-x1Mat[,t])/share[j,t]
      dd1[,,j]=a1
      dz1[,,j]=(x1Mat[j,t]-x1Mat[,t])
      D1[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd1[,,j])
      Dz1[((j-1)*TJ+1):(j*(TJ))]=as.vector(dz1[,,j])


      a2[,t]=x2Mat[j,t]-x2Mat[,t]/share[j,t]
      dd2[,,j]=a2
      dz2[,,j]=(x2Mat[j,t]-x2Mat[,t])
      D2[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd2[,,j])
      Dz2[((j-1)*TJ+1):(j*(TJ))]=as.vector(dz2[,,j])

    }
  }

  D=cbind(S*Dz1,S*Dz2) #TJ*J by 1


  Dz=cbind(Dz1,Dz2)

  dd<-polym(D[,1],D[,2],degree=L,raw=TRUE)
  ddz<-polym(Dz[,1],Dz[,2],degree=L,raw=TRUE)


  basisD=matrix(NA,TJ,length(dd[1,]))
  basisDz=matrix(NA,TJ,length(ddz[1,]))

  for (t in 1:TJ){
    basisD[t,]=colSums(dd[(t*J-J+1):(t*J),]) #basisX is TJ by L
    basisDz[t,]=colSums(ddz[(t*J-J+1):(t*J),]) #basisX is TJ by L
  }


  ######
  dd1=array(NA,c(J,T,J))
  a1=matrix(NA,J,T)
  bX1=matrix(0,(J*J*T),1)
  bXz1=bX1

  dd2=array(NA,c(J,T,J))
  dd1z=dd1
  dd2z=dd2
  a2=matrix(NA,J,T)
  bX2=matrix(0,(J*J*T),1)
  bXz2=bX2

  for (j in 1:J){
    for (t in 1:T){
      a1[,t]=(x1Mat[,t])/s0[t]
      a1[,t]=x1Mat[,t]
      dd1[,,j]=a1
      dd1z[,,j]=x1Mat[,t]
      bX1[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd1[,,j])
      bXz1[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd1z[,,j])

      a2[,t]=x2Mat[,t]/s0[t]
      a2[,t]=x2Mat[,t]
      dd2[,,j]=a2
      dd2z[,,j]=x2Mat[,t]
      bX2[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd2[,,j])
      bXz2[((j-1)*TJ+1):(j*(TJ))]=as.vector(dd2z[,,j])
    }
  }


  bX=cbind(S0*bXz1,S0*bXz2) #TJ*J by 2
  bXz=cbind(bXz1,bXz2)

  xx<-polym(bX[,1],bX[,2],degree=L,raw=TRUE)
  zz<-polym(bXz[,1],bXz[,2],degree=L,raw=TRUE)

  basisX=matrix(NA,TJ,length(xx[1,]))
  basisXz=matrix(NA,TJ,length(zz[1,]))


  for (l in 1:TJ){
    basisX[l,]=colSums(xx[(l*J-J+1):(l*J),]) #basisX is TJ by L
    basisXz[l,]=colSums(zz[(l*J-J+1):(l*J),])
  }


  B=cbind(basisD,basisX)
  IV=cbind(basisDz,basisXz)

  y=R
  X=B
  Z=cbind(x,IV)
  X=cbind(x,B)

  #############IV 1####################
  ginv_PP <- ginv(t(Z)%*%Z)
  pz<- Z%*%ginv_PP%*%t(Z)
  xhat<-pz%*%X

  lm1<-lm(R~xhat -1)
  fit2<-predict(lm1)
  mean((R-fit2)^2)
  coefficients(lm1)[1:2]

  ##########IV 2#######################

  ginv_PP <- ginv(t(IV)%*%IV)
  pz<- IV%*%ginv_PP%*%t(IV)
  Bhat<-pz%*%B
  Xhat=cbind(x,Bhat)
  lm2<-lm(R~Xhat -1)
  fit2<-predict(lm2)
  mean((R-fit2)^2)
  beta0=coefficients(lm2)[1:2]

  lambda_j.Alpha0hat <- .lambda_j.Alpha0hat(B, IV, sigma0_v, tune_type = "CV")
  lambda_j <- lambda_j.Alpha0hat$lambda_j; Alpha0hat <- lambda_j.Alpha0hat$Alpha0hat
  Dhat <- IV %*% Alpha0hat
  XX=cbind(x,Dhat)
  lambda.beta_Lasso <- .lambda.beta_Lasso(y, XX, sigma0_h, no_pen_ids = (1:2))
  lambda <- lambda.beta_Lasso$lambda
  beta_Lasso_XX <- lambda.beta_Lasso$beta_Lasso
  beta1=beta_Lasso_XX[1:2]
  #projected estimator
  beta2=solve(t(x)%*%x)%*%t(x)%*%(y-XX[,3:ncol(XX)]%*%beta_Lasso_XX[3:length(beta_Lasso_XX)])
  #delta_hat

  xihat=y-XX%*%beta_Lasso_XX
  deltahat=x%*%beta1+xihat
  Deltahat=mean(delta-deltahat)

  Betas=cbind(t(beta1),t(beta2))

  A=cbind(Betas,Deltahat)

  data.frame(
    config_id = config_id_,
    trial_id = trial_id,
    beta1_1 = beta1[1],
    beta1_2 = beta1[2],
    beta2_1 = beta2[1],
    beta2_2 = beta2[1],
    Delta.hat = Deltahat
  ) %>%
    write.csv(file = paste("res/res_", config_id_, "_", trial_id, ".csv", sep=''))
}

trial()
