choice<-subset(SamFilterIL, SiteCode == "125_1988_8")
choice<-droplevels(choice)
choice$SigMax<-choice$maxDL/(1+exp((log(19)*(choice$LD50-choice$L50)/(choice$L95-choice$L50))))

param <- c(choice$maxDL,choice$L50,choice$L95,choice$SigMax) # MaxDL, L50, L95, SigMax
Lm50 <- choice$LD50 # estimated size at 50% maturity
#eLML from L50
midpts <- seq(2,210,2)
G <- STM(param,midpts)
Nt <- numeric(105)
Nt[trunc(Lm50/2)] <- 1000
Nt1 <- G %*% (G %*% Nt)
choice$eLML<-(findmedL(Nt1))
Nt1df<-as.data.frame(Nt1)
Nt1df<-add_rownames(Nt1df, "Length")
pick<-which(Nt1df$V1 > 0)
Nt1df <- Nt1df[pick,]
pick<-which(Nt1df$Length >= 138)
U.LML <- Nt1df[pick,]
choice$PctU.LML<-sum(U.LML$V1/10)
#pick<-choice[,c(1,6,33:37)]  


ggplot(Nt1df, aes(x=Length, y = V1/10)) + 
  xlab('Length (mm)')+ 
  ylab('Frequency')+
  geom_density()+
  theme_bw()+#white background
  #scale_fill_grey(start = 0.9, end = 0.1)+
  theme(legend.position="none",
        axis.text.x  = element_text(size=14, angle = 90, vjust=0.4),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14),
        axis.title.x=element_blank())+
  scale_x_continuous(breaks=seq(as.numeric(min(Nt1df$Length)), as.numeric(max(Nt1df$Length)), 10))
ggsave("Mod_DenLML_plot.tiff", width = 9, height = 9, units = "cm")




names(IL.info)[names(IL.info)=='Site_code']<-"Growth_Id"
names(IL.info)[names(IL.info)=='Site.ID']<-"SIT_Id"
names(IL.info)[names(IL.info)=='Site ']<-"Site_Growth"
#Add zone 
IL.info$Zone[IL.info$BlockNo %in%  c(seq(14,30,1)) | IL.info$SubBlockNo %in% c("13C", "13D", "13E", "31A")] <- "E"
IL.info$Zone[IL.info$BlockNo %in%  c(seq(7,12,1)) | IL.info$SubBlockNo %in% c("13A", "13B", "06D", "6D")] <- "W"
IL.info$Zone[IL.info$SubBlockNo %in% c("6A","6B", "6C")] <- "CW"
IL.info$Zone[IL.info$SubBlockNo %in% c("5A", "5B", "5C")] <- "N"
IL.info$Zone[IL.info$BlockNo %in%  c(1, 2, 3, 4,47, 48, 49,39, 40) | IL.info$SubBlockNo %in% c("31B")] <- "N" 
IL.info$Zone[IL.info$BlockNo %in% c(seq(32, 38,1),seq(41,46,1), seq(50,57,1))] <- "BS"




East.IL<-droplevels(subset(IL.info, Zone == 'E'))
East.SamF<-droplevels(subset(SamFilter, Zone == 'E'))

East.SAM.IL<-left_join(East.SamF,East.IL[,4:9], by = 'SIT_Id')
# look at unallocated SiteCodes
East.Unpick<-subset(East.SAM.IL, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

East_UnMatched<-subset(East.SAM.IL, is.na(MaxDL))

#       OUTPUT  1
East_SIT_ID<-subset(East.SAM.IL, !is.na(MaxDL))
Out1<-subset(East_SIT_ID, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))

#++++++++++++++++++++++++++++++++++++++++
#  Dealing with mutliple growth data allocations in subblocks
#++++++++++++++++++++++++++++++++++++++++
SubBlockDupes<-c('13D', '13E', '14A', '16D')
East_SBDupes_IL<-droplevels(subset(East.IL, SubBlockNo %in% SubBlockDupes))
East_SBDupes_SAM<-droplevels(subset(East.Unpick, SubBlockNo %in% SubBlockDupes))

#new column with start of site name
East_SBDupes_SAM$Part_Name<-sapply(strsplit(East_SBDupes_SAM$SIT_Name, "\\ "), `[[`, 1)
East_SBDupes_IL$Part_Name<-sapply(strsplit(as.character(East_SBDupes_IL$Site_Growth), "\\ "), `[[`, 1)

East_Sitebind<-left_join(East_SBDupes_SAM,East_SBDupes_IL[,c(4:5,7:9,12)], by = 'Part_Name')
East.Unpick<-subset(East_Sitebind, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  2
East_Sitebind<-subset(East_Sitebind, !(is.na(MaxDL)))
Out2<-subset(East_Sitebind, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))


Site_IDs<-c(815, 879) # this selects appropriate IL from 13D and 13E for remaining SAM in those subblocks
East_Site_ID_IL<-droplevels(subset(East.IL, SIT_Id %in% Site_IDs))

East_ID<-left_join(East.Unpick,East_Site_ID_IL[,3:9], by = 'SubBlockNo')

East.Unpick<-subset(East_ID, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  3
East_ID<-subset(East_ID, !(is.na(MaxDL)))
Out3<-subset(East_ID, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))

###########
sb<-c(14,16)
pick <- which(East_SBDupes_IL$BlockNo  %in% sb)
East.Siteq<-East_SBDupes_IL[pick,]
SB_mean<-ddply(East.Siteq,.(BlockNo), summarize,  MaxDL = mean(MaxDL, na.rm=T), L50 = mean(L50), L95 = mean(L95))
SB_mean$Growth_Id<-'BlockAvg'
SB_mean$Site_Growth<-'Site_Growth'

#       OUTPUT  4
East_dupes<-left_join(East.Unpick,SB_mean, by = 'BlockNo')
Out4<-subset(East_dupes, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))


#rm(SB_mean,East.Unpick,East_Site_ID_IL, Site_IDs)




# Datasets with unique IL data but no site_ID match

pick <- which(East_UnMatched$SubBlockNo  %in% SubBlockDupes)
East_UnMatched<-East_UnMatched[-pick,]
East_UnMatched<-East_UnMatched[,1:32]

East.SBMatch<-left_join(East_UnMatched,East.IL[,3:9], by = 'SubBlockNo')
East.Unpick<-subset(East.SBMatch, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  5
East.SBMatch<-subset(East.SBMatch, !(is.na(MaxDL)))
Out5<-subset(East_dupes, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))


East.BlkMatch<-left_join(East.Unpick,East.IL[,c(2,4:9)], by = 'BlockNo')
East.Unpick<-subset(East.BlkMatch, is.na(MaxDL))
East.Unpick<-East.Unpick[,1:32]

#       OUTPUT  6
East.BlkMatch<-subset(East.BlkMatch, !(is.na(MaxDL)))
Out6<-subset(East.BlkMatch, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))


East.Unpick$Blk1<-East.Unpick$BlockNo+1
East.BlkPlus1<-merge(East.Unpick,East.IL[,c(2,4:9)], by.x = c("Blk1"), by.y = c('BlockNo'))

#       OUTPUT  7
East.BlkPlus1<-subset(East.BlkPlus1, !(is.na(MaxDL)))
Out7<-subset(East.BlkPlus1, select=c("SiteCode", 'Site_Growth', "Growth_Id", 'MaxDL', 'L50', 'L95'))

#######RBIND ALL SUBSETS
IL_out<-rbind(Out1, Out2, Out3, Out4, Out5, Out6, Out7)
