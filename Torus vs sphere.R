library(TDA)

## define different parameters
## number of replications r is fixed to 100

r<-100
a<-1
b<-2
sample_size<-c(20,seq(50,500,by=50))

## generate data from 3-D Torus(1,2), inner radius 1 and outer radius 2

set.seed(1739)
Data<-torusUnif(20,a,b)
DataDim<-ncol(Data)
Pop_Bettinumber<-c(1,2,1) ### population Betti number 

## Monte Carlo study of Uniform distribution on 3-D torus.

Power<-c()

for (n in sample_size) 
  
{
  set.seed(1557)
  X<-replicate(r,torusUnif(n,a,b),simplify = "array")
  
  ### calculating r values of Test statistic and storing in TV
  
  TV<-c()
  for (l in 1:r)
  {
    #### computing Persistent homology 
    
    maxscale<-1.5
    Dia<-ripsDiag(X[,,l], maxdimension = 2, maxscale , library = "GUDHI")
    
    #### Compute Betti number from Persistent Homology 
    
    index<-which(Dia$diagram[,3]>=maxscale)
    sample_Bettinumber<-c()
    for (i in 0:(DataDim-1)) 
    {
      B<-c()
      for (j in index) 
      {
        if(Dia$diagram[j,1]==i)
        {
          B<-append(B,j,after = length(B))
        }
      }
      sample_Bettinumber<-append(sample_Bettinumber,length(B))
      
    }
    
    
    ### Test Statistic
    
    TS<-sum(abs(sample_Bettinumber-Pop_Bettinumber)) ### L1 Norm Distance
    
    TV<-append(TV,TS,after = length(TV))
  }
  
  ##### Estimation of critical values
  
  Critical_Value <-quantile(TV,0.95) 
  
  ### Monte Carlo study for Uniform distribution on 3-D Sphere  
  
  set.seed(1557)
  XA<-replicate(r,sphereUnif(n, 2, r = 3),simplify = "array")  
  
  ##Calculating test statistic under various alternatives and storing values in TVA
  
  TVA<-c()
  for (l in 1:r) 
  {
    
    maxscale<-2      
    
    #### computing Persistent homology 
    
    Dia<-ripsDiag(XA[,,l], maxdimension = 2, maxscale , library = "GUDHI")
    
    #### Compute Betti number from Persistent Homology 
    
    index<-which(Dia$diagram[,3]>=maxscale)
    Alt_Bettinumber<-c()
    for (i in 0:(DataDim-1)) 
    {
      B<-c()
      for (j in index) 
      {
        if(Dia$diagram[j,1]==i)
        {
          B<-append(B,j,after = length(B))
        }
      }
      Alt_Bettinumber<-append(Alt_Bettinumber,length(B))
      
    }
    
    ### Test Statistic under alternative hypothesis
    
    TSA<-sum(abs(Alt_Bettinumber-Pop_Bettinumber))
    TVA<-append(TVA,TSA,after = length(TVA))
  }
  MPower<-length(which(TVA>=Critical_Value))/r
  Power<- append(Power, MPower,after=length(Power))
}

Power
plot(sample_size,Power,xlab = "Sample Size",ylab = "Power",pch=20,type = "o")


