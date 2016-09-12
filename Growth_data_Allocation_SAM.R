

# 
# choice<-subset(SamFilterIL, SiteCode == "125_1988_8")
# choice<-droplevels(choice)
# choice$SigMax<-choice$maxDL/(1+exp((log(19)*(choice$LD50-choice$L50)/(choice$L95-choice$L50))))
# 
# param <- c(choice$maxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
# Lm50 <- choice$LD50 # estimated size at 50% maturity
# #eLML from L50
# midpts <- seq(2,210,2)
# G <- STM(param,midpts)
# Nt <- numeric(105)
# Nt[trunc(Lm50/2)] <- 1000
# Nt1 <- G %*% (G %*% Nt)
# choice$eLML<-(findmedL(Nt1))
# Nt1df<-as.data.frame(Nt1)
# Nt1df<-add_rownames(Nt1df, "Length")
# pick<-which(Nt1df$V1 > 0)
# Nt1df <- Nt1df[pick,]
# pick<-which(Nt1df$Length >= 138)
# U.LML <- Nt1df[pick,]
# choice$PctU.LML<-sum(U.LML$V1/10)
# #pick<-choice[,c(1,6,33:37)]  
# 
# 
# ggplot(Nt1df, aes(x=Length, y = V1/10)) + 
#   xlab('Length (mm)')+ 
#   ylab('Frequency')+
#   geom_density()+
#   theme_bw()+#white background
#   #scale_fill_grey(start = 0.9, end = 0.1)+
#   theme(legend.position="none",
#         axis.text.x  = element_text(size=14, angle = 90, vjust=0.4),
#         axis.title.y = element_text(size=14),
#         axis.text.y  = element_text(size=14),
#         axis.title.x=element_blank())+
#   scale_x_continuous(breaks=seq(as.numeric(min(Nt1df$Length)), as.numeric(max(Nt1df$Length)), 10))
# ggsave("Mod_DenLML_plot.tiff", width = 9, height = 9, units = "cm")

#use SamFilter as precusror to this file 
# use SamFilter030816.RData 


library(car)
library(MASS)
library(boot)
library(dplyr)
library(tibble)
library(plyr)
library(gdata)
library(ggplot2)
library(multcompView)
library(devtools)
keep(SamFilter, SamResults, sure=T)

#recode database subblock errors
SamFilter$SubBlockNo[SamFilter$SubBlockNo==11] <- "11A"

SamFilterIL<-SamFilter

setwd("D:/Fisheries Research/Abalone/SAM")
IL.info<-read.csv("Inv.Log.Data.csv", header = T)
summary(IL.info)

names(IL.info)[names(IL.info)=='Site_code']<-"Growth_Id"
names(IL.info)[names(IL.info)=='Site.ID']<-"SIT_Id"
names(IL.info)[names(IL.info)=='Site']<-"Site_Growth"
#Add zone 
IL.info$Zone[IL.info$BlockNo %in%  c(seq(14,30,1)) | IL.info$SubBlockNo %in% c("13C", "13D", "13E", "31A")] <- "E"
IL.info$Zone[IL.info$BlockNo %in%  c(seq(7,12,1)) | IL.info$SubBlockNo %in% c("13A", "13B", "06D", "6D")] <- "W"
IL.info$Zone[IL.info$SubBlockNo %in% c("6A","6B", "6C")] <- "CW"
IL.info$Zone[IL.info$SubBlockNo %in% c("5A", "5B", "5C")] <- "N"
IL.info$Zone[IL.info$BlockNo %in%  c(1, 2, 3, 4,47, 48, 49,39, 40) | IL.info$SubBlockNo %in% c("31B")] <- "N" 
IL.info$Zone[IL.info$BlockNo %in% c(seq(32, 38,1),seq(41,46,1), seq(50,57,1))] <- "BS"


#match SAM data to Growth data by site ID
SAM.IL<-left_join(SamFilter,IL.info[,4:9], by = 'SIT_Id')
Unpick<-subset(SAM.IL, is.na(MaxDL))
SIT_ID<-subset(SAM.IL, !is.na(MaxDL))

#anova diferences L%
boxcox(SIT_ID$L50^2.5~SIT_ID$LD50)
fit<-lm(L50^2.5~LD50, data=SIT_ID)
summary(fit)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))


ggplot(data = SIT_ID, aes(x=LD50,  y=L50)) + 
  geom_point()+
  xlab(bquote(''~LM['50%']~'(mm)')) + ylab(bquote(''~L['50%']~'(mm)'))+
  geom_smooth(method=lm, se=F, fill='Black', fullrange=F, size=1.2, color='black')+
  #ggtitle(paste(dum$SubBlockNo, FishYear))+
  #labs(title= Yeardum$SubBlockNo, size=10)+
  #geom_histogram(binwidth=50)+
  theme_bw()+
  scale_color_identity()+ #this makes sure the color follows the color argument above in aes()
  theme(legend.position=c(0.9, 0.8))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14))+
  theme(axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))









East.IL<-droplevels(subset(IL.info, Zone == 'E'))
East.SamF<-droplevels(subset(SamFilter, Zone == 'E'))

East.SAM.IL<-left_join(East.SamF,East.IL[,4:9], by = 'SIT_Id')
# look at unallocated SiteCodes
East.Unpick<-subset(East.SAM.IL, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

East_UnMatched<-subset(East.SAM.IL, is.na(MaxDL)) # used in output 5 onwards

#       OUTPUT  1
East_SIT_ID<-subset(East.SAM.IL, !is.na(MaxDL))

#++++++++++++++++++++++++++++++++++++++++
#  Dealing with mutliple growth data allocations in subblocks
#++++++++++++++++++++++++++++++++++++++++
SubBlockDupes<-c('13D', '13E', '14A', '16D')
East_SBDupes_IL<-droplevels(subset(East.IL, SubBlockNo %in% SubBlockDupes))

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from East.Unpick and extract all obs from SubBlockDupes
East_SBDupes_SAM<-droplevels(subset(East.Unpick, SubBlockNo %in% SubBlockDupes))

#new column with start of site name
East_SBDupes_SAM$Part_Name<-sapply(strsplit(East_SBDupes_SAM$SIT_Name, "\\ "), `[[`, 1)
East_SBDupes_IL$Part_Name<-sapply(strsplit(as.character(East_SBDupes_IL$Site_Growth), "\\ "), `[[`, 1)

East_Sitebind<-left_join(East_SBDupes_SAM,East_SBDupes_IL[,c(4:5,7:9,12)], by = 'Part_Name')
East.Unpick<-subset(East_Sitebind, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  2
East_Sitebind<-subset(East_Sitebind, !(is.na(MaxDL)))

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from East.Unpick and extract 13D and 13E for remaining SAM in those subblocks
Site_IDs<-c(815, 879) # this selects appropriate IL from 13D and 13E for remaining SAM in those subblocks
East_Site_ID_IL<-droplevels(subset(East.IL, SIT_Id %in% Site_IDs))

East_SiteID<-left_join(East.Unpick,East_Site_ID_IL[,c(3:5, 7:9)], by = 'SubBlockNo')

East.Unpick<-subset(East_SiteID, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  3
East_ID<-subset(East_SiteID, !(is.na(MaxDL)))
names(East_ID)[names(East_ID)=='SIT_Id.x']<-"SIT_Id"


#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from blocks 14 and 16 and use average grwoth in those subblocks
###########
sb<-c(14,16)
pick <- which(East_SBDupes_IL$BlockNo %in% sb)
East.BLK<-East_SBDupes_IL[pick,]
SB_mean<-ddply(East.BLK,.(BlockNo), summarize,  MaxDL = mean(MaxDL, na.rm=T), L50 = mean(L50), L95 = mean(L95))
SB_mean$Growth_Id<-'BlockAvg'
SB_mean$Site_Growth<-'Site_Growth'

#       OUTPUT  4
East_dupes<-left_join(East.Unpick,SB_mean, by = 'BlockNo')


#rm(SB_mean,East.Unpick,East_Site_ID_IL, Site_IDs)

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from east unmatched 2ith unique IL data but no site_ID match

pick <- which(East_UnMatched$SubBlockNo %in% SubBlockDupes)
East_UnMatched<-East_UnMatched[-pick,] # remove all those subblocks already dealt with above
East_UnMatched<-East_UnMatched[,1:32]

#match Growth and SAM by subblock
East.SBMatch<-left_join(East_UnMatched,East.IL[,3:9], by = 'SubBlockNo')
East.Unpick<-subset(East.SBMatch, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  5  - matched subblocks
East.SBMatch<-subset(East.SBMatch, !(is.na(MaxDL)))
names(East.SBMatch)[names(East.SBMatch)=='SIT_Id.x']<-"SIT_Id"

#match Growth and SAM by block
East.BlkMatch<-left_join(East.Unpick,East.IL[,c(2,4:9)], by = 'BlockNo')
East.Unpick<-subset(East.BlkMatch, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  6- matched blocks
East.BlkMatch<-subset(East.BlkMatch, !(is.na(MaxDL)))
names(East.BlkMatch)[names(East.BlkMatch)=='SIT_Id.x']<-"SIT_Id"

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from east East.Unpick and match by block +1 for closest growth data
East.Unpick$Blk1<-East.Unpick$BlockNo+1
East.BlkPlus1<-merge(East.Unpick,East.IL[,c(2,4:9)], by.x = c("Blk1"), by.y = c('BlockNo'))

#       OUTPUT  7
East.BlkPlus1<-subset(East.BlkPlus1, !(is.na(MaxDL)))
names(East.BlkPlus1)[names(East.BlkPlus1)=='SIT_Id.x']<-"SIT_Id"

#######RBIND ALL SUBSETS (function below for unmatch columns)

rbind.match.columns <- function(input1, input2) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}

SAM_Out<-rbind.match.columns(East_SIT_ID, East_Sitebind)
SAM_Out<-rbind.match.columns(SAM_Out, East_ID)
SAM_Out<-rbind.match.columns(SAM_Out, East_dupes)
SAM_Out<-rbind.match.columns(SAM_Out, East.SBMatch)
SAM_Out<-rbind.match.columns(SAM_Out, East.BlkMatch)
SAM_Out<-rbind.match.columns(SAM_Out, East.BlkPlus1)


SamOut<-rbind(East.BlkPlus1[,2:39], East.BlkMatch, East.SBMatch, East_dupes, East_ID, East_Sitebind[,-33], East_SIT_ID)
SamOutEast<-SAM_Out[!duplicated(SAM_Out[,1]),]

keep(SamFilter, SamResults, SamOutEast, IL.info, sure=T)

#--------------------------------------------------END EAST

#================================================================
#
#                           WEST
#
#================================================================

W.IL<-droplevels(subset(IL.info, Zone == 'W'))
W.SamF<-droplevels(subset(SamFilter, Zone == 'W'))

W.SAM.IL<-left_join(W.SamF,W.IL[,4:9], by = 'SIT_Id')
# look at unallocated SiteCodes
W.Unpick<-subset(W.SAM.IL, is.na(MaxDL))
W.Unpick<-W.Unpick[,1:32]

W_UnMatched<-subset(W.SAM.IL, is.na(MaxDL)) # used in output 5 onwards

#       OUTPUT  1
W_SIT_ID<-subset(W.SAM.IL, !is.na(MaxDL))

#++++++++++++++++++++++++++++++++++++++++
#  Dealing with mutliple growth data allocations in subblocks
#++++++++++++++++++++++++++++++++++++++++
SubBlockDupes<-c('11A')
W_SBDupes_IL<-droplevels(subset(W.IL, SubBlockNo %in% SubBlockDupes))

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from W.Unpick and extract all obs from SubBlockDupes
W_SBDupes_SAM<-droplevels(subset(W.Unpick, SubBlockNo %in% SubBlockDupes))

#new column with start of site name
W_SBDupes_SAM$Part_Name<-sapply(strsplit(W_SBDupes_SAM$SIT_Name, "\\ "), `[[`, 1)
W_SBDupes_IL$Part_Name<-sapply(strsplit(as.character(W_SBDupes_IL$Site_Growth), "\\ "), `[[`, 1)

W_Sitebind<-left_join(W_SBDupes_SAM,W_SBDupes_IL[,c(4:5,7:9,12)], by = 'Part_Name')
W.Unpick<-subset(W_Sitebind, is.na(MaxDL))
W.Unpick<-W.Unpick[,1:32]

#       OUTPUT  2
W_Sitebind<-subset(W_Sitebind, !(is.na(MaxDL)))

#++++++++++++++++++++++++++++++++++++++++
# Manually assign growth within subblock 
W.Unpick$Growth_Id<-c("GI", "GI", "GI","TR")

#       OUTPUT  3
W_SiteID<-left_join(W.Unpick,W_SBDupes_IL[,c(4:5, 7:9)], by = 'Growth_Id')

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from W unmatched 2ith unique IL data but no site_ID match

pick <- which(W_UnMatched$SubBlockNo %in% SubBlockDupes)
W_UnMatched<-W_UnMatched[-pick,] # remove all those subblocks already dealt with above
W_UnMatched<-W_UnMatched[,1:32]

#match Growth and SAM by subblock
W.SBMatch<-left_join(W_UnMatched,W.IL[,3:9], by = 'SubBlockNo')
W.Unpick<-subset(W.SBMatch, is.na(MaxDL))
W.Unpick<-W.Unpick[,1:32]

#       OUTPUT  5  - matched subblocks
W.SBMatch<-subset(W.SBMatch, !(is.na(MaxDL)))
names(W.SBMatch)[names(W.SBMatch)=='SIT_Id.x']<-"SIT_Id"


# allocate 11E grwoth data from 11c
W.pick<-subset(W.Unpick, SubBlockNo == '11E')
W.pickIL<-subset(W.IL, SubBlockNo == '11C')

W.SB11E<-left_join(W.pick,W.pickIL[,c(2,4:9)], by = 'BlockNo')


W.pick<-subset(W.Unpick, SubBlockNo != '11E')



#match Growth and SAM by block
W.BlkMatch<-left_join(W.Unpick,W.IL[,c(2,4:9)], by = 'BlockNo')
W.Unpick<-subset(W.BlkMatch, is.na(MaxDL))
W.Unpick<-W.Unpick[,1:32]

#       OUTPUT  6- matched blocks
W.BlkMatch<-subset(W.BlkMatch, !(is.na(MaxDL)))
names(W.BlkMatch)[names(W.BlkMatch)=='SIT_Id.x']<-"SIT_Id"

#++++++++++++++++++++++++++++++++++++++++
# Take the unmatched obs from W W.Unpick and match by block +1 for closest growth data
W.Unpick$Blk1<-W.Unpick$BlockNo+1
W.BlkPlus1<-merge(W.Unpick,W.IL[,c(2,4:9)], by.x = c("Blk1"), by.y = c('BlockNo'))

#       OUTPUT  7
W.BlkPlus1<-subset(W.BlkPlus1, !(is.na(MaxDL)))
names(W.BlkPlus1)[names(W.BlkPlus1)=='SIT_Id.x']<-"SIT_Id"

#######RBIND ALL SUBSETS
IL_out<-rbind(Out1, Out2, Out3, Out4, Out5, Out6, Out7)

rbind.match.columns <- function(input1, input2) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}

SAM_Out<-rbind.match.columns(W_SIT_ID, W_Sitebind)
SAM_Out<-rbind.match.columns(SAM_Out, W_ID)
SAM_Out<-rbind.match.columns(SAM_Out, W_dupes)
SAM_Out<-rbind.match.columns(SAM_Out, W.SBMatch)
SAM_Out<-rbind.match.columns(SAM_Out, W.BlkMatch)
SAM_Out<-rbind.match.columns(SAM_Out, W.BlkPlus1)


SamOut<-rbind(W.BlkPlus1[,2:39], W.BlkMatch, W.SBMatch, W_dupes, W_ID, W_Sitebind[,-33], W_SIT_ID)
SamOutW<-SAM_Out[!duplicated(SAM_Out[,1]),]

keep(SamFilter, SamResults, SamOutW, sure=T)

#================================================================
#
#                           NORTH
#
#================================================================

N.IL<-droplevels(subset(IL.info, Zone == 'N'))
N.SamF<-droplevels(subset(SamFilter, Zone == 'N'))

N.SAM.IL<-left_join(N.SamF,N.IL[,4:9], by = 'SIT_Id')
# look at unallocated SiteCodes
N.Unpick<-subset(N.SAM.IL, is.na(MaxDL))
N.Unpick<-N.Unpick[,1:32]

N_UnMatched<-subset(N.SAM.IL, is.na(MaxDL)) # used in output 5 onwards

#       OUTPUT  1
N_SIT_ID<-subset(N.SAM.IL, !is.na(MaxDL))
