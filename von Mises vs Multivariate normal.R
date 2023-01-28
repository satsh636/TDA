## Two sample problem Mixture of von Mises and Multivarite normal
library(TDA)
library(mvtnorm)


## number of replications
r<-10 

## dimension of data
d<-3

## sample size
sample_size<-c(20,seq(50,200,by=50))

## multivariate normal parameters

mu<-replicate(d,0)
sigma<-cbind(c(1,0.5,0.5),c(0.5,1,0.5),c(0.5,0.5,1))

## parameters of multivariate von mises
theta<-c(3, 4, 5) * rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))

## scale the data by euclidean norm
Enorm<-function(v)
{
  w<-(sum(v^2))^(1/2)
  return(v/w)
  
}

## Computing Betti numbers of Rips complex
## D is ripsdiag$diagram

sample_Betti<-function(D,e)
{
  index<-which(D[,3]>=e)
  Betti_num<-c()
  for (i in 0:(d-1)) 
  {
    B<-c()
    for (j in index) 
    {
      if(D[j,1]==i)
      {
        B<-append(B,j,after = length(B))
      }
    }
    Betti_num<-append(Betti_num,length(B))
    
  }
  return(Betti_num)
}

## Perform Monte Carlo study

Power<-c()
for (n in sample_size) 
  
{
  
  maxscale<-(1/n)^(1/d) 
  e<-maxscale
   
  set.seed(1306)
  Data_1<-replicate(r,Enorm(rmovMF(n,theta)),simplify = "array")
  Data_2<-replicate(r,Enorm(rmovMF(n,theta)),simplify = "array")
  Data_Alt<-replicate(r,Enorm(rmvnorm(n,mu,sigma)),simplify = "array")
  
  
  TV<-c()
  TVA<-c()
  for (l in 1:r) 
  {
    
    #### computing Persistent homology
    Dia_1<-ripsDiag(Data_1[,,l], maxdimension = 1, maxscale , library = "GUDHI")$diagram
    Dia_2<-ripsDiag(Data_2[,,l], maxdimension = 1, maxscale , library = "GUDHI")$diagram
    Dia_Alt<-ripsDiag(Data_Alt[,,l], maxdimension = 1, maxscale , library = "GUDHI")$diagram
    
    ### Test Statistic
    
    TS<-sum(abs(sample_Betti(Dia_1,e)-sample_Betti(Dia_2,e)))  ### L1 Norm Distance
    
    TV<-append(TV,TS,after = length(TV))
    
    TSA<-sum(abs(sample_Betti(Dia_1,e)-sample_Betti(Dia_Alt,e)))
    
    TVA<-append(TVA,TSA,after = length(TVA))
    
  }
  
  Critical_Value <-quantile(TV,(0.95)/2)
  MPower<-length(which(TVA>=Critical_Value))/r
  Power<- append(Power, MPower,after=length(Power))
  
}

Power
plot(sample_size,Power,xlab = "Sample Size",ylab = "Power",type = "l")  
