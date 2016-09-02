# Version 05/06/2014  

rm(list=ls(all=TRUE))
## SET THE WORKING AND RESULTS DIRECTORIES
wkdir <- "D:/R_Stuff/SAM"
setwd(wkdir)

## Load raw csv file
#infile <- "Basic Sam Query.csv"

library(car)
library(MASS)
library(boot)
library(dplyr)
library(plyr)
library(gdata)
library(ggplot2)
library(multcompView)
library(devtools)
#sources
#source("D:/GitCode/AbResearch/SAM_utils_TechReport.R")
source("D:/GitCode/AbResearch/Grwth_matrix.r")

# # load D:\R_Stuff\SAM    SAM2016_April  SAM_TechReport230816.RData
# keep(SamResults, sure =T)

#rename columns
#SamResults<-rename(SamResults, c("SIT_Latitude"="Latitude", "SIT_Longitude"="Longitude", "SIT_StatBlock"="BlockNo", "SIT_SubBlock"="SubBlockNo"))
names(SamResults)[names(SamResults)=='SIT_Latitude']<-"Latitude"
names(SamResults)[names(SamResults)=='SIT_Longitude']<-"Longitude"
names(SamResults)[names(SamResults)=='SIT_StatBlock']<-"BlockNo"
names(SamResults)[names(SamResults)=='SIT_SubBlock']<-"SubBlockNo"

#Add columns for sample size (n), and % of sample <L05 <L50 and >L95
SamResults$n<-SamResults$I+SamResults$M
SamResults$PctL05<-SamResults$N.underLD05/SamResults$n*100
SamResults$PctL50<-SamResults$N.underLD50/SamResults$n*100
SamResults$PctL95<-SamResults$N.overLD95/SamResults$n*100

#add column for year of sample
SamResults$FishYear<-format(SamResults$SAM_Date,'%Y')
SamResults$FishYear<-as.factor(SamResults$FishYear)

#Reformat SublockNo
SamResults$SubBlockNo<-paste(SamResults$BlockNo,SamResults$SubBlockNo, sep="")

#Add zone 
SamResults$Zone[SamResults$BlockNo %in%  c(seq(14,30,1)) | SamResults$SubBlockNo %in% c("13C", "13D", "13E", "31A")] <- "E"
SamResults$Zone[SamResults$BlockNo %in%  c(seq(7,12,1)) | SamResults$SubBlockNo %in% c("13A", "13B", "06D", "6D")] <- "W"
SamResults$Zone[SamResults$SubBlockNo %in% c("6A","6B", "6C")] <- "CW"
SamResults$Zone[SamResults$SubBlockNo %in% c("5A", "5B", "5C")] <- "N"
SamResults$Zone[SamResults$BlockNo %in%  c(1, 2, 3, 4,47, 48, 49,39, 40) | SamResults$SubBlockNo %in% c("31B")] <- "N" 
SamResults$Zone[SamResults$BlockNo %in% c(seq(32, 38,1),seq(41,46,1), seq(50,57,1))] <- "BS"

#Add LML 
SamResults$LML[SamResults$BlockNo %in%  c(14:26,29:30) | SamResults$SubBlockNo %in% c("13C", "13D", "13E", "31A")] <- 138
SamResults$LML[SamResults$BlockNo %in%  c(seq(7,12,1)) | SamResults$SubBlockNo %in% c("13A", "13B", "06D", "6D")] <- 140
SamResults$LML[SamResults$SubBlockNo %in% c("6A","6B", "6C", "5D")] <- 132
SamResults$LML[SamResults$BlockNo %in%  c(1, 2, 3, 4,39, 40) |SamResults$SubBlockNo %in% c("5A", "5B", "5C", "31B", "49D")] <- 127
SamResults$LML[SamResults$BlockNo %in% c("27", "28")] <- 145
SamResults$LML[SamResults$BlockNo %in% c(seq(32, 38,1),seq(41,46,1), seq(50,57,1))] <- 110
SamResults$LML[SamResults$BlockNo %in%  c(47,48) |SamResults$SubBlockNo %in% c("49A", "49B","49C")] <- 120


#Ld bootstrap range
SamResults$Ld50BootRange<-SamResults$Ld50BootU95-SamResults$Ld50BootL95
SamResults$Ld95BootRange<-SamResults$Ld95BootU95-SamResults$Ld95BootL95

#remove files which break metarules of <5% or >95% mature
pick <- which(SamResults$PctL05 <5)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)
pick <- which(SamResults$PctL95 <5)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)

#DT Seasonal SAM MSc used the same sites for SAM work between 1999-2001. There is some concern that by resampling the same sites continuously for 
#for 2.5 years there might be effects on the L50 outside of those found at all other sites. Therefore all but the 1st samples in 1999 and 2000 are removed.

DT_MSc<-c(171, 172)
pick <- which(SamResults$SIT_Id == DT_MSc)
DT_MSc<-SamResults[pick,]

pick <- which(DT_MSc$FishYear == 2001)
DT_MSc<-DT_MSc[-pick,]
DT_MSc[order(as.Date(DT_MSc$SAM_Date, format="%d/%m/%Y")),]


Sites<-unique(DT_MSc$SIT_Name)
if (exists("DT_out")) 
  rm(DT_out)

for(d in Sites){
  choice<-subset(DT_MSc, SIT_Name == d)
  pick<-ddply(choice,.(FishYear),function(x) head(x,1))
  if (exists("DT_out"))
    DT_out <- rbind(DT_out, pick)
  else
    DT_out <- pick
}

#Remove all DT_MSc sites and replace with DT_out

SamFilter<-droplevels(subset(SamResults, SIT_Id != 171))
SamFilter<-droplevels(subset(SamFilter, SIT_Id != 172))

SamFilter<-rbind(SamFilter, DT_out)

# save(SamResults, file="SamResults.Rdata")
# write.csv(SamResults, file='SamResultsLatLong.csv')


#remove duplicate records

SamFilter<-SamFilter[!duplicated(SamFilter[,1]),]

#
###############################################
#Calculate the Size transition matrix and eLML
###############################################
#
source("D:/GitCode/AbResearch/Grwth_matrix.r")


Sites<-unique(SamResults$SiteCode)

#####
#     L50%
#####
if (exists("eLMLResults")) 
  rm(eLMLResults)

for(i in Sites){
  choice<-subset(SamResults, SiteCode == i)
  choice$L50<-1.1539*choice$LD50-15.335
  choice$L95<-1.0862*choice$LD50+32.461
  choice$MaxDL<-0.46095*choice$L95-0.46856*choice$L50+5.58943
  choice$SigMax<-choice$MaxDL/(1+exp((log(19)*(choice$LD50-choice$L50)/(choice$L95-choice$L50))))
  
  param <- c(choice$MaxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
  Lm50 <- choice$LD50 # estimated size at 50% maturity
  #eLML from L50
  midpts <- seq(2,210,2)
  G <- STM(param,midpts)
  Nt <- numeric(105)
  Nt[trunc(Lm50/2)] <- 1000
  Nt1 <- G %*% (G %*% Nt)
  choice$eLML<-(findmedL(Nt1))

  pick<-choice[,c(1,6,33:37)]  
  if (exists("eLMLResults"))
    eLMLResults <- rbind(eLMLResults, pick)
  else
    eLMLResults <- pick
}
#pass back to GwthResults
GwthResults<-eLMLResults



#####
#     L50% UCI as base for eLML
#####
if (exists("eLMLResults")) 
  rm(eLMLResults)

for(i in Sites){
  choice<-subset(SamResults, SiteCode == i)
  if (is.na(choice$Ld50BootU95)) next
  choice$L50<-1.1539*choice$Ld50BootU95-15.335
  choice$L95<-1.0862*choice$Ld50BootU95+32.461
  choice$MaxDL<-0.46095*choice$L95-0.46856*choice$Ld50BootU95+5.58943
  choice$SigMax<-choice$MaxDL/(1+exp((log(19)*(choice$Ld50BootU95-choice$L50)/(choice$L95-choice$L50))))
  
  param <- c(choice$MaxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
  #eLML from Lower CI L50
  Lm50 <- choice$Ld50BootU95 # estimated size at 50% maturity
  midpts <- seq(2,210,2)
  G <- STM(param,midpts)
  Nt <- numeric(105)
  Nt[trunc(Lm50/2)] <- 1000
  Nt1 <- G %*% (G %*% Nt)
  choice$eLMLbootU95<-(findmedL(Nt1))
  
  pick<-choice[,c(1,37)]  
  if (exists("eLMLResults"))
    eLMLResults <- rbind(eLMLResults, pick)
  else
    eLMLResults <- pick
}
#join to GwthResults
GwthResults<-join(GwthResults, eLMLResults, by ="SiteCode", type ="left" )

#
##############################
#           L95% elml
##############################
#
if (exists("eLMLResults")) 
  rm(eLMLResults)

for(i in Sites){
  choice<-subset(SamResults, SiteCode == i)
  choice$L50<-1.1539*choice$LD95-15.335
  choice$L95<-1.0862*choice$LD95+32.461
  choice$MaxDL<-0.46095*choice$L95-0.46856*choice$L50+5.58943
  choice$SigMax<-choice$MaxDL/(1+exp((log(19)*(choice$LD95-choice$L50)/(choice$L95-choice$L50))))
  
  param <- c(choice$MaxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
  Lm50 <- choice$LD95 # estimated size at 50% maturity
  #eLML from L50
  midpts <- seq(2,210,2)
  G <- STM(param,midpts)
  Nt <- numeric(105)
  Nt[trunc(Lm50/2)] <- 1000
  Nt1 <- G %*% (G %*% Nt)
  choice$L95eLML<-(findmedL(Nt1))
  
  pick<-choice[,c(1,37)]  
  if (exists("eLMLResults"))
    eLMLResults <- rbind(eLMLResults, pick)
  else
    eLMLResults <- pick
}
#pass back to GwthResults
GwthResults<-join(GwthResults, eLMLResults, by ="SiteCode", type ="left" )

#####
#     L95% UCI
#####

if (exists("eLMLResults")) 
  rm(eLMLResults)
if (exists("choice")) 
  rm(choice)

for(i in Sites){
  choice<-subset(SamResults, SiteCode == i)
  if (is.na(choice$Ld95BootU95)) next
  choice$L50<-1.1539*choice$Ld95BootU95-15.335
  choice$L95<-1.0862*choice$Ld95BootU95+32.461
  choice$MaxDL<-0.46095*choice$L95-0.46856*choice$Ld95BootU95+5.58943
  choice$SigMax<-choice$MaxDL/(1+exp((log(19)*(choice$Ld95BootU95-choice$L50)/(choice$L95-choice$L50))))
  
  param <- c(choice$MaxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
  #eLML from Lower CI L50
  Lm50 <- choice$Ld95BootU95 # estimated size at 50% maturity
  midpts <- seq(2,210,2)
  G <- STM(param,midpts)
  Nt <- numeric(105)
  Nt[trunc(Lm50/2)] <- 1000
  Nt1 <- G %*% (G %*% Nt)
  choice$L95eLMLbootU95<-(findmedL(Nt1))
  
  pick<-choice[,c(1,37)]  
  if (exists("eLMLResults"))
    eLMLResults <- rbind(eLMLResults, pick)
  else
    eLMLResults <- pick
}
#join to GwthResults
GwthResults<-join(GwthResults, eLMLResults, by ="SiteCode", type ="left" )

#pass back to SamResults
SamResults<-join(SamResults, GwthResults[,c(1,7:10)], by ="SiteCode", type ="left" )

#difference between l50 and LML
SamResults$LMLDiff<-SamResults$LML-SamResults$eLML


BlockSumStats<-ddply(SamResults,.(BlockNo, Zone), summarize,  n = length(SiteCode), 
                        mn.L50 = mean(LD50, na.rm=T), mn.LCI50 = mean(Ld50BootL95, na.rm=T), mn.UCI50 = mean(Ld50BootU95, na.rm=T),
                        mn.L95 = mean(LD95, na.rm=T), mn.LCI95 = mean(Ld95BootL95, na.rm=T), mn.UCI95 = mean(Ld95BootU95, na.rm=T),
                        mn.IQR = mean(IQR, na.rm=T), sd.IQR = sd(IQR, na.rm=T),
                        mn.pct.L50 = mean(PctL50, na.rm=T), sd.pct.L50 = mean(PctL50, na.rm=T),
                        mn.Bootrange.L50 = mean(Ld50BootRange, na.rm=T), sd.Bootrange.L50 = mean(Ld50BootRange, na.rm=T),
                        mn.eLML = mean(eLML, na.rm=T) , sd.eLML = sd(eLML, na.rm=T),
                        mn.eLMLbootL95 = mean(eLMLbootL95, na.rm=T) , sd.eLMLbootL95 = sd(eLMLbootL95, na.rm=T),
                        mn.eLMLbootU95 = mean(eLMLbootU95, na.rm=T) , sd.eLMLbootU95 = sd(eLMLbootU95, na.rm=T),
                        diffLML = mean(LMLDiff, na.rm=T))
#write.csv(BlockSumStats, file='SAMBlockSummaryStats.csv')

ZoneSumStats<-ddply(SamResults,.(Zone), summarize,  n = length(SiteCode), 
                    mn.L50 = mean(LD50, na.rm=T), mn.LCI50 = mean(Ld50BootL95, na.rm=T), mn.UCI50 = mean(Ld50BootU95, na.rm=T),
                    mn.L95 = mean(LD95, na.rm=T), mn.LCI95 = mean(Ld95BootL95, na.rm=T), mn.UCI95 = mean(Ld95BootU95, na.rm=T),
                    mn.IQR = mean(IQR, na.rm=T), sd.IQR = sd(IQR, na.rm=T),
                    mn.pct.L50 = mean(PctL50, na.rm=T), sd.pct.L50 = mean(PctL50, na.rm=T),
                    mn.Bootrange.L50 = mean(Ld50BootRange, na.rm=T), sd.Bootrange.L50 = mean(Ld50BootRange, na.rm=T),
                    mn.eLML = mean(eLML, na.rm=T) , sd.eLML = sd(eLML, na.rm=T),
                    mn.eLMLbootL95 = mean(eLMLbootL95, na.rm=T) , sd.eLMLbootL95 = sd(eLMLbootL95, na.rm=T),
                    mn.eLMLbootU95 = mean(eLMLbootU95, na.rm=T) , sd.eLMLbootU95 = sd(eLMLbootU95, na.rm=T))

write.csv(SamResults, file='SAMResultsLatLong.csv')

rm(choice,DT_out,pick)
