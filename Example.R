rm(list=ls())
library(CDVineCopulaConditional)
library(foreach)
library(parallel)
library(fitdistrplus)
library(VineCopula)
library(pracma)
library(doParallel)
library(fGarch)
set.seed(123)

path0<-getwd()
setwd(path0)

year_month<-vector()
for(year1 in 1950:2022){
  for(month1 in 1:12)
  {
    year_month<-rbind(year_month,cbind(year1,month1))
  }
}
year_month<-as.data.frame(year_month)
colnames(year_month)<-c("year","month")

norm1<-function(x){
  p1<-fitdist(x,"norm",method="mle")
  x2<-pnorm(x,p1$estimate[1],p1$estimate[2])
  return(x2)
}

source("t.R")
source("gaussian.R")
source("gumbel.R")
source("clayton.R")
source("frank.R")
source("DynamicVineCopula.R")


spei1<-read.table("spei.txt")
sti1<-read.table("sti.txt")
ssi1<-read.table("ssi.txt")

family1<-matrix(0,3,3)
par1<-matrix(0,3,3)
par2<-matrix(0,3,3)
matrix1<-matrix(0,3,3)

for(monthi in 3:5)
{
  k<-which(year_month$month==monthi)
  
  data1<-cbind(norm1(ssi1$V1[k]),norm1(spei1$V1[k]),norm1(sti1$V1[k]))
  colnames(data1)<-c("ssi","spei","sti")
  
  vine1<-CDVineCondFit(data1,familyset=c(-0,-6,-7,-8,-9,-10,-13,-16,-18,-20,-23,-24,-26,-27,-28,-29,-30,-33,-34,
                                         -36,-37,-38,-39,-40,-104,-114,-124,-134,-204,-214,-224,-234),Nx=2,treecrit="BIC",selectioncrit="AIC")
  
  family1[1,monthi-2]<-vine1$family[3,1]
  family1[2,monthi-2]<-vine1$family[3,2]
  family1[3,monthi-2]<-vine1$family[2,1]
  
  par1[1,monthi-2]<-vine1$par[3,1]
  par1[2,monthi-2]<-vine1$par[3,2]
  par1[3,monthi-2]<-vine1$par[2,1]
  
  par2[1,monthi-2]<-vine1$par2[3,1]
  par2[2,monthi-2]<-vine1$par2[3,2]
  par2[3,monthi-2]<-vine1$par2[2,1]
  
  matrix1[1,monthi-2]<-as.numeric(paste(vine1$Matrix[1,1],vine1$Matrix[3,1],sep=""))
  matrix1[2,monthi-2]<-as.numeric(paste(vine1$Matrix[2,2],vine1$Matrix[3,2],sep=""))
  matrix1[3,monthi-2]<-as.numeric(paste(vine1$Matrix[1,1],vine1$Matrix[2,1],vine1$Matrix[3,1],sep=""))
  
}

for(monthi in 3:5)
{
  result1<-DynamicVineCopula_main(data1,family1,matrix1,par1,par2)
  
  TVPDEP1<-result1$TVPDEP
  
  for(TVPDEPi in 1:3){
    path1<-paste("TVPDEP-",monthi,"-",TVPDEPi,".txt",sep="")
    value_t<-TVPDEP1[[TVPDEPi]]
    write.table(value_t,path1,row.names = F,col.names = F)
  }
}


