## SAMTechReportBootstrapped 
# Version 15/08/16 
# this is a reworked version of the CarlyNewBoot.R


rm(list=ls(all=TRUE))
## SET THE WORKING AND RESULTS DIRECTORIES
wkdir <- "D:/Fisheries Research/Abalone"
setwd(wkdir)

# loctaion for graph outputs
resdir <- "D:/Fisheries Research/Abalone/SAM"


## Load raw csv file
infile <- "DataForExport.csv"

library(MASS)
library(boot)
library(car)
library(dplyr)
source("D:/GitCode/AbResearch/SAM_utils_TechReport.R")

samdata <- read.csv(infile,header=T)
dim(samdata)
head(samdata)

#Option - remove Wineglass Bay Site for Season * region analysis
pick <- which(samdata$Site_Name == "Wineglass Bay")
samdata <- samdata[-pick,]
samdata <- droplevels(samdata)

# Simplify Maturity classes
samdata$Mat <- NA
samdata$Mat <- ifelse(samdata$Sex=="I", c("I"), c("M"))

# Characterize Sites
Sites <- unique(samdata$Site_Name); Sites
NS <- length(Sites)

#recode date

samdata$NewDate <- as.character(samdata$Sample_Date)
samdata$NewDate <- substr(samdata$NewDate, 1, 10)

samdata$FishMonth<-sapply(strsplit(samdata$NewDate, "/"), "[", 2)
samdata$FishMonth<-as.numeric(samdata$FishMonth)
samdata$FishQtr<-samdata$FishMonth
samdata$FishQtr[samdata$FishQtr %in% c("1", "2", "3")] <- "Q1"
samdata$FishQtr[samdata$FishQtr %in% c("4", "5", "6")] <- "Q2"
samdata$FishQtr[samdata$FishQtr %in% c("7", "8", "9")] <- "Q3"
samdata$FishQtr[samdata$FishQtr %in% c("10", "11", "12")] <- "Q4"
samdata$FishQtr<-as.factor(samdata$FishQtr)
summary(samdata)

#PORT davey data

# PDavey<-read.csv('SAM_Port_davey_2016.csv',header=T)
# summary(PDavey)
# # Simplify Maturity classes
# PDavey$Mat <- NA
# PDavey$Mat <- ifelse(PDavey$Sex=="I", c("I"), c("M"))
# # Characterize Sites
# PDavey$Region<-PDavey$Location
#PDavey$Site_Name<-PDavey$Location
# Sites <- unique(PDavey$Region); Sites
# NS <- length(Sites)
# samdata<-PDavey


# Code Season and Site
samdata$Season <- samdata$FishQtr
samdata$Region <- samdata$Site_Name

## Outlier Removal: Unusual small mature animals
pick <- which((samdata$Length <60) & (samdata$Mat=="M"))
outlier <- rbind(outlier,samdata[pick,])
samdata <- samdata[-pick,]

## Outlier Removal: Unusual large immature animals
pick <- which((samdata$Length >139) & (samdata$Mat=="I"))
outlier <- samdata[pick,]
samdata <- samdata[-pick,]

## Allocate into Size Classes

samdata$SizeC <- NA
pick <- which(samdata$Length <= 90)
samdata$SizeC[pick] <- "Small"
pick <- which((samdata$Length > 90) & (samdata$Length <=130))
samdata$SizeC[pick] <- "Medium"
pick <- which(samdata$Length > 130)
samdata$SizeC[pick] <- "Large"


#source("D:\\Students\\Carly Giosio/carly_utils.R")

## Calculate Base cases for each site.
starttime <- unclass(Sys.time())
sizeCl <- c("Small","Medium","Large")
as.character(Sites)


# 1: George III Rock April
# 2: George III Rock January
# 3: St Helens April
# 4: St Helens January
# 5: Wineglass Bay January

##
# ## Step 1: Choose a site based on the number index above
# pSite <- 1            #  CHANGE THIS TO SELECT WHICH Location TO USE
# #  pSite <- 1  will analyse gardens
# #  pSite <- 3  will analyse GIII
# Sites
# 
# Site1 <- Sites[pSite]
# pick <- which(samdata$Site_Name == Site1)
# if (length(pick) > 0) {
#  subdata <- samdata[pick,]
# } else { "Oops something has gone wrong"
# }
Regions<-unique(samdata$Region)

## Calculate Base cases for each site.
columns <- c("LM50", "LCI50", "UCI50", "LM75", "LCI75", "UCI75", "LM95","LCI95", "UCI95", "IQ","a","b","smallN","mediumN","largeN","totalN")
BaseResults <- matrix(0,nrow=NS,ncol=length(columns),dimnames=list(Regions,columns))

for (pSite in 1:NS) {
    Site <- Regions[pSite]
    pick <- which(samdata$Site_Name == Site)
    if (length(pick) > 0) {
      subdata <- samdata[pick,]
    } else { "Oops something has gone wrong"
    }

##
## Step 3: Run from here to the very bottom
reps <- 1000
 columns <- c("LM50","IQ","a","b","smallN","mediumN","largeN","totalN")
 bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))

sampledDat <- split(subdata,subdata$SizeC)
catSizeLarge <- nrow(sampledDat[[1]])
catSizeMedium <- nrow(sampledDat[[2]])
catSizeSmall <- nrow(sampledDat[[3]])
pickfrom1 <- seq(1,catSizeLarge,1)
pickfrom2 <- seq(1,catSizeMedium,1)
pickfrom3 <- seq(1,catSizeSmall,1)

for (iter in 1:reps) {
 ## do the random sampling
 adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
 adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
 adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
 
 subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
 subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
#   pickfrom <- seq(1,catSize,1)
#    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
#   subdata2 <-subdata[picked,]
 ## Create required fields for logistic regression
 SizeMat <- table(subdata2$Length, subdata2$Mat)
 SizeMat <- as.data.frame(rbind(SizeMat))
 SizeMat$Length <- as.numeric(rownames(SizeMat))
 head(SizeMat)
 SizeMat$Total <- NA
 SizeMat$Total <- SizeMat$I + SizeMat$M
 SizeMat$MatRatio <- NA
 SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
 out <- doLogistic(SizeMat)
 bootResults[iter,1] <- out$LM50
 bootResults[iter,2] <- out$IQ
 model <- out$Model
 bootResults[iter,3:4] <- model$coef
 scN <- numeric(3)
 pick <- which(SizeMat$Length <= 90)
 scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
 pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
 scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
 pick <- which(SizeMat$Length > 130)
 scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
 bootResults[iter,5:8] <- c(scN,sum(scN))
 #  plotgraph(SizeMat,out$LM50,Site1,scN,savefile=T)
}


droplevels(Sites[pSite])
p <- c(2.5, 50, 97.5)/100
bootResults <- as.data.frame(bootResults)
BaseResults[pSite,1]<-quantile(bootResults$LM50,0.5)
BaseResults[pSite,2]<-quantile(bootResults$LM50,0.05)
BaseResults[pSite,3]<-quantile(bootResults$LM50,0.95)
BaseResults[pSite,10] <- mean(out$IQ)
BaseResults[pSite,11] <- mean(bootResults[pSite,3])
BaseResults[pSite,12] <- mean(bootResults[pSite,4])
BaseResults[pSite,13] <- mean(bootResults[pSite,5])
BaseResults[pSite,14] <- mean(bootResults[pSite,6])
BaseResults[pSite,15] <- mean(bootResults[pSite,7])
BaseResults[pSite,16] <- mean(bootResults[pSite,8])

# # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
# LM50boot<-as.data.frame(quantile(bootResults$LM50, p))

##

#####LM75


## Step 3: Run from here to the very bottom
reps <- 1000
columns <- c("LM75","IQ","a","b","smallN","mediumN","largeN","totalN")
bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))

sampledDat <- split(subdata,subdata$SizeC)
catSizeLarge <- nrow(sampledDat[[1]])
catSizeMedium <- nrow(sampledDat[[2]])
catSizeSmall <- nrow(sampledDat[[3]])
pickfrom1 <- seq(1,catSizeLarge,1)
pickfrom2 <- seq(1,catSizeMedium,1)
pickfrom3 <- seq(1,catSizeSmall,1)

for (iter in 1:reps) {
  ## do the random sampling
  
  
  adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
  adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
  adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
  
  subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
  subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
  #   pickfrom <- seq(1,catSize,1)
  #    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
  #   subdata2 <-subdata[picked,]
  ## Create required fields for logistic regression
  SizeMat <- table(subdata2$Length, subdata2$Mat)
  SizeMat <- as.data.frame(rbind(SizeMat))
  SizeMat$Length <- as.numeric(rownames(SizeMat))
  head(SizeMat)
  SizeMat$Total <- NA
  SizeMat$Total <- SizeMat$I + SizeMat$M
  SizeMat$MatRatio <- NA
  SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
  out <- doLogistic(SizeMat)
  bootResults[iter,1] <- out$LM75
  bootResults[iter,2] <- out$IQ
  model <- out$Model
  bootResults[iter,3:4] <- model$coef
  scN <- numeric(3)
  pick <- which(SizeMat$Length <= 90)
  scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
  pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
  scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
  pick <- which(SizeMat$Length > 130)
  scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
  bootResults[iter,5:8] <- c(scN,sum(scN))
  #  plotgraph(SizeMat,out$LM50,Site1,scN,savefile=T)
}


droplevels(Sites[pSite])
p <- c(2.5, 50, 97.5)/100
bootResults <- as.data.frame(bootResults)
BaseResults[pSite,4]<-quantile(bootResults$LM75,0.5)
BaseResults[pSite,5]<-quantile(bootResults$LM75,0.05)
BaseResults[pSite,6]<-quantile(bootResults$LM75,0.95)
# # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
# LM75boot<-as.data.frame(quantile(bootResults$LM75, p))

##
## Step 3: Run from here to the very bottom
reps <- 1000
columns <- c("LM95","IQ","a","b","smallN","mediumN","largeN","totalN")
bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))

sampledDat <- split(subdata,subdata$SizeC)
catSizeLarge <- nrow(sampledDat[[1]])
catSizeMedium <- nrow(sampledDat[[2]])
catSizeSmall <- nrow(sampledDat[[3]])
pickfrom1 <- seq(1,catSizeLarge,1)
pickfrom2 <- seq(1,catSizeMedium,1)
pickfrom3 <- seq(1,catSizeSmall,1)

for (iter in 1:reps) {
  ## do the random sampling
  
  
  adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
  adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
  adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
  
  subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
  subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
  #   pickfrom <- seq(1,catSize,1)
  #    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
  #   subdata2 <-subdata[picked,]
  ## Create required fields for logistic regression
  SizeMat <- table(subdata2$Length, subdata2$Mat)
  SizeMat <- as.data.frame(rbind(SizeMat))
  SizeMat$Length <- as.numeric(rownames(SizeMat))
  head(SizeMat)
  SizeMat$Total <- NA
  SizeMat$Total <- SizeMat$I + SizeMat$M
  SizeMat$MatRatio <- NA
  SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
  out <- doLogistic(SizeMat)
  bootResults[iter,1] <- out$LM95
  bootResults[iter,2] <- out$IQ
  model <- out$Model
  bootResults[iter,3:4] <- model$coef
  scN <- numeric(3)
  pick <- which(SizeMat$Length <= 90)
  scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
  pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
  scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
  pick <- which(SizeMat$Length > 130)
  scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
  bootResults[iter,5:8] <- c(scN,sum(scN))
  #  plotgraph(SizeMat,out$LM50,Site1,scN,savefile=F)
}


droplevels(Sites[pSite])
p <- c(2.5, 50, 97.5)/100
bootResults <- as.data.frame(bootResults)
BaseResults[pSite,7]<-quantile(bootResults$LM95,0.5)
BaseResults[pSite,8]<-quantile(bootResults$LM95,0.05)
BaseResults[pSite,9]<-quantile(bootResults$LM95,0.95)

plotgraph(SizeMat,out$LM50,out$LM75, out$LM95,Site,scN,savefile=F)

# # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
# LM95boot<-as.data.frame(quantile(bootResults$LM95, p))
# 
}
####
#### HI #############
####



HIdata<-read.csv("HI_SAM2014.csv", header=T)
summary(HIdata)

# Simplify Maturity classes
HIdata$Mat <- NA
HIdata$Mat <- ifelse(HIdata$Sex=="I", c("I"), c("M"))

# Characterize Sites
Sites <- unique(HIdata$Region_Name); Sites
NS <- length(Sites)

## Allocate into Size Classes

HIdata$SizeC <- NA
pick <- which(HIdata$Length <= 80)
HIdata$SizeC[pick] <- "Small"
pick <- which((HIdata$Length > 80) & (HIdata$Length <=110))
HIdata$SizeC[pick] <- "Medium"
pick <- which(HIdata$Length > 110)
HIdata$SizeC[pick] <- "Large"


#source("D:\\Students\\Carly Giosio/carly_utils.R")

## Calculate Base cases for each site.
starttime <- unclass(Sys.time())
sizeCl <- c("Small","Medium","Large")
as.character(Sites)


Regions <- unique(HIdata$Region); 
NS <-  length(Regions)

#Regions<-unique(HIdata$Region_Name)


## Calculate Base cases for each site.
columns <- c("LM50", "LCI50", "UCI50", "LM75", "LCI75", "UCI75", "LM95","LCI95", "UCI95", "IQ","a","b","smallN","mediumN","largeN","totalN")
BaseResults <- matrix(0,nrow=NS,ncol=length(columns),dimnames=list(Regions,columns))

for (pSite in 1:NS) {
  Site <- Regions[pSite]
  pick <- which(HIdata$Region == Site)
  if (length(pick) > 0) {
    subdata <- HIdata[pick,]
  } else { "Oops something has gone wrong"
  }

  ## Step 3: Run from here to the very bottom
  reps <- 1000
  columns <- c("LM50","IQ","a","b","smallN","mediumN","largeN","totalN")
  bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))
  
  sampledDat <- split(subdata,subdata$SizeC)
  catSizeLarge <- nrow(sampledDat[[1]])
  catSizeMedium <- nrow(sampledDat[[2]])
  catSizeSmall <- nrow(sampledDat[[3]])
  pickfrom1 <- seq(1,catSizeLarge,1)
  pickfrom2 <- seq(1,catSizeMedium,1)
  pickfrom3 <- seq(1,catSizeSmall,1)
  
  for (iter in 1:reps) {
    ## do the random sampling
    adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
    adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
    adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
    
    subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
    subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
    #   pickfrom <- seq(1,catSize,1)
    #    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
    #   subdata2 <-subdata[picked,]
    ## Create required fields for logistic regression
    SizeMat <- table(subdata2$Length, subdata2$Mat)
    SizeMat <- as.data.frame(rbind(SizeMat))
    SizeMat$Length <- as.numeric(rownames(SizeMat))
    head(SizeMat)
    SizeMat$Total <- NA
    SizeMat$Total <- SizeMat$I + SizeMat$M
    SizeMat$MatRatio <- NA
    SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
    out <- doLogistic(SizeMat)
    bootResults[iter,1] <- out$LM50
    bootResults[iter,2] <- out$IQ
    model <- out$Model
    bootResults[iter,3:4] <- model$coef
    scN <- numeric(3)
    pick <- which(SizeMat$Length <= 90)
    scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
    scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which(SizeMat$Length > 130)
    scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
    bootResults[iter,5:8] <- c(scN,sum(scN))
    #  plotgraph(SizeMat,out$LM95,Site,scN,savefile=F)
  }
  
  
  #droplevels(Sites[pSite])
  p <- c(2.5, 50, 97.5)/100
  bootResults <- as.data.frame(bootResults)
  BaseResults[pSite,1]<-quantile(bootResults$LM50,0.5)
  BaseResults[pSite,2]<-quantile(bootResults$LM50,0.05)
  BaseResults[pSite,3]<-quantile(bootResults$LM50,0.95)
  BaseResults[pSite,10] <- mean(out$IQ)
  BaseResults[pSite,11] <- mean(bootResults[pSite,3])
  BaseResults[pSite,12] <- mean(bootResults[pSite,4])
  BaseResults[pSite,13] <- mean(bootResults[pSite,5])
  BaseResults[pSite,14] <- mean(bootResults[pSite,6])
  BaseResults[pSite,15] <- mean(bootResults[pSite,7])
  BaseResults[pSite,16] <- mean(bootResults[pSite,8])
  
  # # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
  # LM50boot<-as.data.frame(quantile(bootResults$LM50, p))
  
  ##
  
  #####LM75
  
  
  ## Step 3: Run from here to the very bottom
  reps <- 1000
  columns <- c("LM75","IQ","a","b","smallN","mediumN","largeN","totalN")
  bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))
  
  sampledDat <- split(subdata,subdata$SizeC)
  catSizeLarge <- nrow(sampledDat[[1]])
  catSizeMedium <- nrow(sampledDat[[2]])
  catSizeSmall <- nrow(sampledDat[[3]])
  pickfrom1 <- seq(1,catSizeLarge,1)
  pickfrom2 <- seq(1,catSizeMedium,1)
  pickfrom3 <- seq(1,catSizeSmall,1)
  
  for (iter in 1:reps) {
    ## do the random sampling
    
    
    adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
    adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
    adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
    
    subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
    subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
    #   pickfrom <- seq(1,catSize,1)
    #    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
    #   subdata2 <-subdata[picked,]
    ## Create required fields for logistic regression
    SizeMat <- table(subdata2$Length, subdata2$Mat)
    SizeMat <- as.data.frame(rbind(SizeMat))
    SizeMat$Length <- as.numeric(rownames(SizeMat))
    head(SizeMat)
    SizeMat$Total <- NA
    SizeMat$Total <- SizeMat$I + SizeMat$M
    SizeMat$MatRatio <- NA
    SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
    out <- doLogistic(SizeMat)
    bootResults[iter,1] <- out$LM75
    bootResults[iter,2] <- out$IQ
    model <- out$Model
    bootResults[iter,3:4] <- model$coef
    scN <- numeric(3)
    pick <- which(SizeMat$Length <= 90)
    scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
    scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which(SizeMat$Length > 130)
    scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
    bootResults[iter,5:8] <- c(scN,sum(scN))
    #  plotgraph(SizeMat,out$LM50,Site1,scN,savefile=T)
  }
  
  
  #droplevels(Sites[pSite])
  p <- c(2.5, 50, 97.5)/100
  bootResults <- as.data.frame(bootResults)
  BaseResults[pSite,4]<-quantile(bootResults$LM75,0.5)
  BaseResults[pSite,5]<-quantile(bootResults$LM75,0.05)
  BaseResults[pSite,6]<-quantile(bootResults$LM75,0.95)
  # # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
  # LM75boot<-as.data.frame(quantile(bootResults$LM75, p))
  
  ##
  ## Step 3: Run from here to the very bottom
  reps <- 1000
  columns <- c("LM95","IQ","a","b","smallN","mediumN","largeN","totalN")
  bootResults <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(seq(1,reps,1),columns))
  
  sampledDat <- split(subdata,subdata$SizeC)
  catSizeLarge <- nrow(sampledDat[[1]])
  catSizeMedium <- nrow(sampledDat[[2]])
  catSizeSmall <- nrow(sampledDat[[3]])
  pickfrom1 <- seq(1,catSizeLarge,1)
  pickfrom2 <- seq(1,catSizeMedium,1)
  pickfrom3 <- seq(1,catSizeSmall,1)
  
  for (iter in 1:reps) {
    ## do the random sampling
    
    
    adddata1 <- sort(sample(pickfrom1,catSizeLarge,replace=TRUE))
    adddata2 <- sort(sample(pickfrom2,catSizeMedium,replace=TRUE))
    adddata3 <- sort(sample(pickfrom3,catSizeSmall,replace=TRUE))
    
    subdata2 <- rbind(sampledDat[[1]][adddata1,],sampledDat[[2]][adddata2,])    # combine the bits
    subdata2 <- rbind(subdata2,sampledDat[[3]][adddata3,])    # combine the bits
    #   pickfrom <- seq(1,catSize,1)
    #    picked <- sort(sample(pickfrom,catSize,replace=TRUE))
    #   subdata2 <-subdata[picked,]
    ## Create required fields for logistic regression
    SizeMat <- table(subdata2$Length, subdata2$Mat)
    SizeMat <- as.data.frame(rbind(SizeMat))
    SizeMat$Length <- as.numeric(rownames(SizeMat))
    head(SizeMat)
    SizeMat$Total <- NA
    SizeMat$Total <- SizeMat$I + SizeMat$M
    SizeMat$MatRatio <- NA
    SizeMat$MatRatio <- SizeMat$M/SizeMat$Total
    out <- doLogistic(SizeMat)
    bootResults[iter,1] <- out$LM95
    bootResults[iter,2] <- out$IQ
    model <- out$Model
    bootResults[iter,3:4] <- model$coef
    scN <- numeric(3)
    pick <- which(SizeMat$Length <= 90)
    scN[1] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which((SizeMat$Length > 90) & (SizeMat$Length <= 130))
    scN[2] <- sum(SizeMat$Total[pick],na.rm=T)
    pick <- which(SizeMat$Length > 130)
    scN[3] <- sum(SizeMat$Total[pick],na.rm=T)
    bootResults[iter,5:8] <- c(scN,sum(scN))
    #  plotgraph(SizeMat,out$LM50,Site1,scN,savefile=T)
  }
  
  
  #droplevels(Sites[pSite])
  p <- c(2.5, 50, 97.5)/100
  bootResults <- as.data.frame(bootResults)
  BaseResults[pSite,7]<-quantile(bootResults$LM95,0.5)
  BaseResults[pSite,8]<-quantile(bootResults$LM95,0.05)
  BaseResults[pSite,9]<-quantile(bootResults$LM95,0.95)
  # # Obtain the 95% CI of Lm50 - Lower CI=2.5% & upper CI=97.5% 
  # LM95boot<-as.data.frame(quantile(bootResults$LM95, p))
  # 
  # boots<-cbind(LM50boot, LM75boot, LM95boot)
  # colnames(boots) <- c("L50%", "L75%", "L95%")
  # boots <- add_rownames(boots, "C.I.")
  # 
  # boots$Site<-s
  # 
  # 
  # if (exists("BootsOut")){BootsOut <- rbind(BootsOut, boots)} else { BootsOut <- boots}
  plotgraph(SizeMat,out$LM50,out$LM75, out$LM95,Site,scN,savefile=F)
}