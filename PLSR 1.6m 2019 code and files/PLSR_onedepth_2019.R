################
# PLSR for 1.6m MUX, 2019
# Rachel Corrigan
# June 2020
###############

library(pls)  #Load the pls package
library(lubridate)
library(dplyr)
library(stringr)
library(scales)
pathD<-"/Users/rachelcorrigan/Dropbox/fcr/MagicData/" #Specify folder where data is located
setwd("/Users/rachelcorrigan/Dropbox/fcr/MagicData/")

################################################################################
#Read in FCR WQ data
################################################################################

WQ<-"onedepth2019WQoverlaps.csv" #Specify file where data is located  
dataWQ<-read.table(file=paste(pathD,WQ,sep=""),sep=",",header = TRUE)  #Import data as .csv file
dataWQ=dataWQ[,-c(16:24)]
dataWQ=dataWQ[,-c(1:4)]
#dataWQ = dataWQ[,-10]
WQparam <- c("TN_ugL","TP_ugL","NH4_ugL","NO3NO2_ugL","SRP_ugL","DOC_mgL","DIC_mgL", "DC_mgL", "DN_mgL")   


################################################################################
#### Reading of  FingerPrint (FP) file corresponding to lab concentrations for 
#### calibration hence dataCalFP
#### 
FPcaldata_name<-"onedepth2019FPoverlaps.csv"
dataCalFP<-read.delim(file=paste(pathD,FPcaldata_name,sep=""),sep=",")  #Import data as .csv file
#dataCalFP=dataCalFP[,-c(2)]
colnames(dataCalFP)<-c("ID","Date/Time","status",seq(200,730,2.5)) #Add column names

#subset if doing depth specific model
#dataCalFP <- dataCalFP %>%
#  filter(ID == 1 | ID == 2 | ID == 3) #%>%
#  filter(ID == 2) %>%
#  filter(ID == 3)
#
timesCalFP<-dataCalFP[,1:2] 
dataCalFP<-dataCalFP[,-c(1:3)]
#timesCalFP<-data.matrix(timesCalFP)#date from scan and port number. I removed port 2 here because there is no corresponding depth for wq samples
dataCalFP<-dataCalFP[,-214:-216] #Remove NO3-N values and NAs at high wavelengths
dataCalFP<-data.matrix(dataCalFP) #Convert to data matrix

################################################################################
#### This replaces the ID and Date from the original dataWQ with the exact values
#### from the SCAN so that manual values can be plotted later on in the TS plots
dataWQ$ID<-timesCalFP[,1]
dataWQ$DateTime<-timesCalFP[,2]
#dataWQ<-dataWQ[,(-13)]


################################################################################
#### Reading of  FingerPrint (FP) file corresponding to the entire time series (TS) 
#### 
TimeSeriesFP_name<-"onedepth2019FullTS.csv"
TS_FP<-read.table(file=paste(pathD,TimeSeriesFP_name,sep=""),sep=",", skip=1)  #Import Time Series data as .csv file
colnames(TS_FP)<-c("port","Date","status",seq(200,730,2.5)) #Add column names
TS_FP$Date = as.POSIXct(TS_FP$Date, format = "%m/%d/%y %H:%M")
#TS_FP <- TS_FP %>%
#  filter(port != 10 | port !=11 | port !=12) #%>%
#  filter(port !=11) %>%
#  filter(port !=12)

Dat<-strptime(TS_FP$Date, format = "%Y-%m-%d %H:%M:%S") #Create record of date and time

################################################################################
####  Create matrix to store calculated concentrationss:TS_conc
TS_conc<-as.data.frame(matrix(0,dim(TS_FP)[1],11))  #Create data frame for date/time and predicted NO3-N values
TS_conc[,1]<-TS_FP$port
TS_conc[,2]<-as.character(Dat, "%Y-%m-%d %H:%M:%S")
colnames(TS_conc)<-c("port","DateTime",WQparam) #Add column names
TS_FP<-TS_FP[,(-1:-3)]
TS_FP<-data.matrix(TS_FP) #Convert spectrometer output to matrix


#######################
#Specify number of components for wq param
ncomp=9 #7 for TP and TN

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

param<-"DOC_mgL"
ncomp=6
PLSR_SCAN(param,dataCalFP,dataWQ,TS_FP,ncomp, yesplot=TRUE)

TS_conc$DateTime <- as.POSIXct(TS_conc$DateTime, format="%Y-%m-%d %H:%M:%S")
dataWQ$DateTime <- as.POSIXct(dataWQ$DateTime, format="%m/%d/%y %H:%M")
colnames(dataWQ)[2] <- "Depth"
TS_conc$Depth[TS_conc$port == 2] = 1.6
#WQP_TS <- 10^(WQP_TS)
TS_conc$uncerDOC_max <- NA
TS_conc$uncerDOC_min <- NA
TS_conc$uncerDOC_max <- WQP_TS + 1.96*sd(as.numeric(fit$residuals[,,6]))
TS_conc$uncerDOC_min <- WQP_TS - 1.96*sd(as.numeric(fit$residuals[,,6]))
TS_conc$uncerDOC_max <- unlist(TS_conc$uncerDOC_max)
TS_conc$uncerDOC_min <- unlist(TS_conc$uncerDOC_min)

#sd(as.numeric(fit$residuals[,,6]))
#TS_conc[,(7)]<-WQP_TS #for SRP
TS_conc[,(8)]<-WQP_TS #for DOC
#TS_conc[,(5)]<-10^(WQP_TS) #for NH4, used LOG10 in PLSR!
#TS_conc[,(6)]<-10^(WQP_TS) #for N03N02, used LOG10 in PLSR!
#TS_conc[,(9)]<-WQP_TS #for DIC
#TS_conc[,(10)]<-WQP_TS #for DC
#TS_conc[,(11)]<-WQP_TS #for DN

#TS_conc <- read.csv("TS_PREDICTIONS_2019_16m.csv")

srp_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=SRP_ugL), size=0.8, alpha=0.8) +
 # scale_colour_manual(name='', values=c('Predict'='black', 'Observe'='blue'), guide=guide_legend(),
 #                     labels=c("Predict", "Observe"))+
  geom_ribbon(data=TS_conc, aes(ymin=uncerSRP_min, ymax=uncerSRP_max, x=DateTime, fill = "purple"), alpha = 0.2)+
  #scale_fill_identity(name = '', labels = c('uncert'), guide=guide_legend()) +
   geom_point(data=dataWQ, aes(x=DateTime, y=SRP_ugL), colour="blue") +
  #scale_shape_manual(values=c(16),guide=guide_legend(override.aes=list(colour="blue",size=c(2))))+
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("SRP (", mu, "g/L)")), title = "Predicted vs. Obs SRP at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
srp_test


doc_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=DOC_mgL), size=0.8, alpha=0.8) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerDOC_min, ymax=uncerDOC_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=DOC_mgL), colour="blue") +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("DOC mg/L)")), title = "Predicted vs. Obs DOC at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
doc_test

NH4_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=NH4_ugL), size=0.8, alpha=0.8) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerNH4_min, ymax=uncerNH4_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=NH4_ugL), colour="blue", size=1.5) +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("NH4 (", mu, "g/L)")), title = "Predicted vs. Obs NH4 at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position = "topright")
NH4_test

NO3_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=NO3NO2_ugL), size=0.8, alpha=0.8) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerNO3_min, ymax=uncerNO3_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=NO3NO2_ugL), colour="blue") +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("NO3NO2 (", mu, "g/L)")), title = "Predicted vs. Obs NO3NO2 at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
NO3_test


DIC_test <- ggplot() +
  geom_line(data=TS_conc, aes(x=DateTime,y=DIC_mgL)) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerDIC_min, ymax=uncerDIC_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=DIC_mgL), colour="blue") +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("DIC mg/L)")), title = "Predicted vs. Obs DIC at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
DIC_test

DC_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=DC_mgL), size=0.8, alpha=0.8) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerDC_min, ymax=uncerDC_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=DC_mgL), colour="blue") +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("DC mg/L)")), title = "Predicted vs. Obs DC at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
DC_test

DN_test <- ggplot() +
  geom_point(data=TS_conc, aes(x=DateTime,y=DN_mgL), size=0.8, alpha=0.8) +
  geom_ribbon(data=TS_conc, aes(ymin=uncerDN_min, ymax=uncerDN_max, x=DateTime, fill = "band"), alpha = 0.2)+
  geom_point(data=dataWQ, aes(x=DateTime, y=DN_mgL), colour="blue") +
  #ylim(-5, 50)+
  labs(x="Date", y = expression(paste("DN mg/L)")), title = "Predicted vs. Obs DN at 1.6m, 2019") +
  scale_x_datetime(labels = date_format("%Y-%m-%d"))+
  theme(legend.position="none")
DN_test

#save TS_conc with predictions filled in! 
write.csv(TS_conc, "TS_PREDICTIONS_2019_16m.csv")
TS_conc_chem <- TS_conc[,-13:-28]
write.csv(TS_conc_chem, "TS_PREDICTIONS2019_16m_chem_only.csv")
TS_conc$Reservoir<- "FCR"
TS_conc$Site <- "50"
TS_conc$Depth <- 1.6

TS_format <- TS_conc[,c(29,30,2,12,5:11)]
write.csv(TS_format, "TS_PREDICTIONS2019_FORMATTED.csv")

all_chem_plot <- ggarrange(doc_test, srp_test, NH4_test, NO3_test, DC_test, DN_test, ncol=2, nrow=3)
all_chem_plot
