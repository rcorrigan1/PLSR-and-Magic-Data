##############
# Rachel Corrigan
# 2018 FCR PLSR 
# May 2020
##############

library(pls)  #Load the pls package
library(lubridate)
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
pathD<-"/Users/rachelcorrigan/Dropbox/fcr/MagicData/" #Specify folder where data is located
setwd("/Users/rachelcorrigan/Dropbox/fcr/MagicData/")

################################################################################
#Read in FCR WQ data
################################################################################

WQ<-"allWQoverlap2018.csv"  #Specify file where data is located  
dataWQ<-read.table(file=paste(pathD,WQ,sep=""),sep=",",header = TRUE)  #Import data as .csv file
#dataWQ=dataWQ[,-1]
WQparam <- c("TN_ugL","TP_ugL","NH4_ugL","NO3NO2_ugL","SRP_ugL","DOC_mgL","TFe_mgL","TMn_mgL","SFe_mgL","SMn_mgL")   
dataWQ <- dataWQ[,(-1)]

################################################################################
#### Reading of  FingerPrint (FP) file corresponding to lab concentrations for 
#### calibration hence dataCalFP
#### 
FPcaldata_name<-"allFP_WQ_Overlap_FINAL.csv"
dataCalFP<-read.delim(file=paste(pathD,FPcaldata_name,sep=""),sep=",")  #Import data as .csv file
#dataCalFP=dataCalFP[,-1]
#dataCalFP$ID <- seq.int(nrow(dataCalFP))
#dataCalFP <- dataCalFP[,c(218,1:217)]
colnames(dataCalFP)<-c("ID","Date/Time","status",seq(200,730,2.5)) #Add column names
timesCalFP<-cbind(data.matrix(dataCalFP[,2]),dataCalFP[,1]) #date from scan and port number. I removed port 2 here because there is no corresponding depth for wq samples
dataCalFP<-dataCalFP[,-1:-3] #Remove NO3-N values and NAs at high wavelengths
dataCalFP<-data.matrix(dataCalFP) #Convert to data matrix
dataCalFP<-dataCalFP[,(-214)]

################################################################################
#### This replaces the ID and Date from the original dataWQ with the exact values
#### from the SCAN so that manual values can be plotted later on in the TS plots
dataWQ$ID<-timesCalFP[,2]
dataWQ$DateTime<-timesCalFP[,1]
#dataWQ<-dataWQ[,(-13)]

################################################################################
#### Reading of  FingerPrint (FP) file corresponding to the entire time series (TS) 
#### 
TimeSeriesFP_name<-"FP_with_Port_ALL2018.csv"
TS_FP<-read.table(file=paste(pathD,TimeSeriesFP_name,sep=""),sep=",", skip=1)  #Import Time Series data as .csv file
colnames(TS_FP)<-c("port","Date","status",seq(200,730,2.5)) #Add column names
TS_FP$Date = as.POSIXct(TS_FP$Date, format = "%m/%d/%y %H:%M")
TS_FP<-TS_FP[!(TS_FP$port==11),]
Dat<-strptime(TS_FP$Date, format = "%Y-%m-%d %H:%M:%S") #Create record of date and time

################################################################################
####  Create matrix to store calculated concentrationss:TS_conc
TS_conc<-as.data.frame(matrix(0,dim(TS_FP)[1],12))  #Create data frame for date/time and predicted NO3-N values
TS_conc[,1]<-TS_FP$port
TS_conc[,2]<-as.character(Dat, "%Y-%m-%d %H:%M:%S")
colnames(TS_conc)<-c("port","DateTime",WQparam) #Add column names
TS_FP<-TS_FP[,(-1:-3)]
TS_FP<-data.matrix(TS_FP) #Convert spectrometer output to matrix

#######################
#Specify number of components for wq param
ncomp=7 #7 for TP and TN

################################################################################
#### function which does the calibration and then calculates conc from the TS for
#### a given chemical parameter (param).  It does the calibration for a given 
#### number of components (ncomp)
################################################################################

PLSR_SCAN<-function(param,dataCalFP,dataWQ,TS_FP,ncomp,yesplot=FALSE){
  
  WQ<-data.matrix(subset(dataWQ,select=param)) #Make matrix of the param values
  temp<-cbind(dataCalFP,WQ) #combines FP and WQ columns to remove the rows containing NAs 
  temp<-temp[complete.cases(temp),] #removes the rows containing NAs
  WQ<-data.matrix(subset(temp,select=param)) # recreate a data matrix from the WQ vector minus the NAs
  dataFP<-temp[,-dim(temp)[2]]  # redefines the FP matrix rid off the NA values of missing WQ
  
  
  fit<-plsr(WQ~data.matrix(dataFP),ncomp=ncomp,validation="CV")  #PLSR model to predict param with cross validation
  summary(fit)  #See summary of PLSR model to choose number of components
  Pfit<-predict(fit,dataFP,ncomp=ncomp,type=c("response")) #Predict concentrations based on PLSR model
  #x11()
  WQP<-as.data.frame(matrix(0,1,dim(dataFP)[1])) #Create data frame for predicted param values
  WQP<-as.data.frame(Pfit[1:length(Pfit)])  #Insert predicted param values into data frame
  
  Pfit_TS<-predict(fit,TS_FP,ncomp=ncomp,type=c("response"))
  WQP_TS<-as.data.frame(matrix(0,1,dim(TS_FP)[1])) #Create data frame for predicted Time series values
  WQP_TS<-as.data.frame(Pfit_TS[1:length(Pfit_TS)])  #Insert predicted param values into data frame
  
  if (yesplot==TRUE){
    plot(WQ,as.matrix(WQP),
         xlab=paste("measured",param,"?g/L",sep=" "),
         ylab=c("PLSR_predicted")) #Compare predicted and lab values of param
    
    fit2<-lm(WQ~as.matrix(WQP)) #Linear regression of predicted and lab NO3-N values
    abline(fit2)
    summary(fit2)
    tp_resid <- resid(fit2)
  }
  
  assign("WQP_TS",WQP_TS,env=.GlobalEnv)
}



################################################################################
#### Example for running the function for one parameter only
################################################################################

param<-"TP_ugL"
ncomp=7
PLSR_SCAN(param,dataCalFP,dataWQ,TS_FP,ncomp, yesplot=TRUE)

#plot residuals
tp_resid <- fit$residuals[,,7]
plot(as.matrix(WQP), fit$residuals[,,7], ylab="Residuals", xlab="Predicted", main="Spread of residuals 2018 ")
abline(0, 0)
sd(tp_resid)


TS_conc[,(4)]<-WQP_TS #for TP

TS_conc$uncerTP_max <- NA
TS_conc$uncerTP_min <- NA
TS_conc$uncerTP_max <- WQP_TS + 1.96*sd(as.numeric(fit$residuals[,,7]))
TS_conc$uncerTP_min <- WQP_TS - 1.96*sd(as.numeric(fit$residuals[,,7]))
TS_conc$uncerTP_max <- unlist(TS_conc$uncerTP_max)
TS_conc$uncerTP_min <- unlist(TS_conc$uncerTP_min)

TS_conc$DateTime <- as.POSIXct(TS_conc$DateTime, format="%Y-%m-%d %H:%M:%S")
dataWQ$DateTime <- as.POSIXct(dataWQ$DateTime, format="%m/%d/%y %H:%M")



TS_conc$Depth = NA
TS_conc$Depth[TS_conc$port == 1] = 0.1
TS_conc$Depth[TS_conc$port == 2] = 1.0
TS_conc$Depth[TS_conc$port == 3] = 1.6
TS_conc$Depth[TS_conc$port == 4] = 3.8
TS_conc$Depth[TS_conc$port == 5] = 5.0
TS_conc$Depth[TS_conc$port == 6] = 6.2
TS_conc$Depth[TS_conc$port == 7] = 8.0
TS_conc$Depth[TS_conc$port == 8] = 9.0
TS_conc$Depth[TS_conc$port == 9] = 9.5

#TS_conc_crop = TS_conc %>%
#  filter(DateTime )
#dataWQ<- dataWQ[,(-1)]
#colnames(dataWQ)[2] <- "DateTime"
colnames(dataWQ)[3] <- "Depth"

tp_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=TP_ugL, group=factor(Depth)), size=0.5) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerTP_min, ymax=uncerTP_max, x=DateTime, fill = "band"), alpha = 0.3)+
  geom_point(data=dataWQ, aes(x=DateTime, y=TP_ugL), colour="blue") +
  labs(x="Date", y = expression(paste("Total P (", mu, "g/L)"))) +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  facet_wrap(~Depth, ncol=1) +
  theme(legend.position="none")

plot(tp_test)
