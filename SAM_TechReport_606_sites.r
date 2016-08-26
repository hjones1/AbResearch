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
#sources
source("D:/GitCode/r-AbSpatialAnalyses/GraphsUtils.r")

source("D:/GitCode/AbResearch/SAM_utils_TechReport.R")


# load D:\R_Stuff\SAM    SAM2016_April  SAM_TechReport230816.RData
keep(SamResults, sure =T)

colnames(SamResults)[21] <- "BlockNo"
colnames(SamResults)[22] <- "SubBlockNo"
SamResults$n<-SamResults$I+SamResults$M
SamResults$PctL50<-SamResults$N.underLD50/SamResults$n*100

SamResults$FishYear<-format(SamResults$SAM_Date,'%Y')
SamResults$FishYear<-as.factor(SamResults$FishYear)

SamResults$SubBlockNo<-paste(SamResults$BlockNo,SamResults$SubBlockNo, sep="")


n.per.sample<-ddply(SamResults,.(SiteCode, BlockNo, SubBlockNo), summarize,  n = n)

subdata<-droplevels(subset(n.per.sample, n>=100))
#recode date


# samdata$FishMonth<-sapply(strsplit(samdata$NewDate, "/"), "[", 2)
# samdata$FishMonth<-as.numeric(samdata$FishMonth)
# samdata$FishQtr<-samdata$FishMonth
# samdata$FishQtr[samdata$FishQtr %in% c("1", "2", "3")] <- "Q1"
# samdata$FishQtr[samdata$FishQtr %in% c("4", "5", "6")] <- "Q2"
# samdata$FishQtr[samdata$FishQtr %in% c("7", "8", "9")] <- "Q3"
# samdata$FishQtr[samdata$FishQtr %in% c("10", "11", "12")] <- "Q4"
# samdata$FishQtr<-as.factor(samdata$FishQtr)
# summary(samdata)


SamResults$Zone[SamResults$BlockNo %in%  c(seq(14,30,1)) | SamResults$SubBlockNo %in% c("13C", "13D", "13E", "31A")] <- "E"
SamResults$Zone[SamResults$BlockNo %in%  c(seq(7,12,1)) | SamResults$SubBlockNo %in% c("13A", "13B", "06D", "6D")] <- "W"
SamResults$Zone[SamResults$SubBlockNo %in% c("6A","6B", "6C")] <- "CW"
SamResults$Zone[SamResults$SubBlockNo %in% c("5A", "5B", "5C")] <- "N"
SamResults$Zone[SamResults$BlockNo %in%  c(1, 2, 3, 4,47, 48, 49,39, 40) | SamResults$SubBlockNo %in% c("31B")] <- "N" 
SamResults$Zone[SamResults$BlockNo %in% c(seq(32, 38,1),seq(41,46,1), seq(50,57,1))] <- "BS"

#Ld bootstrap range
SamResults$Ld50BootRange<-SamResults$Ld50BootU95-SamResults$Ld50BootL95
SamResults$Ld95BootRange<-SamResults$Ld95BootU95-SamResults$Ld95BootL95

#remove outliers
pick <- which(SamResults$Ld50BootRange > 50)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)

pick <- which(SamResults$Ld95BootRange > 50)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)



save(SamResults, file="SamResults.Rdata")



SubBlockSumStats<-ddply(SamResults,.(BlockNo, Zone), summarize,  n = length(SiteCode), 
                        mn.L50 = mean(LD50, na.rm=T), mn.LCI50 = mean(Ld50BootL95, na.rm=T), mn.UCI50 = mean(Ld50BootU95, na.rm=T),
                        mn.L95 = mean(LD95, na.rm=T), mn.LCI95 = mean(Ld95BootL95, na.rm=T), mn.UCI95 = mean(Ld95BootU95, na.rm=T),
                        mn.IQR = mean(IQR, na.rm=T), sd.IQR = sd(IQR, na.rm=T),
                        mn.pct.L50 = mean(PctL50, na.rm=T), sd.pct.L50 = mean(PctL50, na.rm=T),
                        mn.Bootrange.L50 = mean(Ld50BootRange, na.rm=T), sd.Bootrange.L50 = mean(Ld50BootRange, na.rm=T),
                        mn.MLStc = mean(MLStc, na.rm=T) , sd.MLStc = sd(MLStc, na.rm=T))
write.csv(SubBlockSumStats, file='SAMBlockSummaryStats.csv')


####RESULTS ANALYSIS
##
#
#

#Option - remove Wineglass Bay Site for Season * region analysis
pick <- which(SamResults$PctL50 >= 78 & SamResults$L50Rci >= 5)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)

pick <- which(SamResults$PctL50 < 10 & SamResults$L50Rci >= 19)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)

pick <- which(SamResults$L50Rci >= 15)
SamResults <- SamResults[-pick,]
SamResults <- droplevels(SamResults)


ggplot(data = SamResults, aes(x=PctL50,  y=Ld50BootRange)) + #color=ifelse(L50Rci>12, 'red', 'black')
  geom_point(color="grey20")+
  xlab(bquote('% Sample <'~L['50%']~'.')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
  geom_smooth(method=loess, se=F, fill='Black', fullrange=F, size=1.2, color='black')+
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

fitd <- loess(L50Rci~PctL50, data=SamResults)
summary(fitd)

pick <- which(SamResults$PctL50 >= 30)
Minus30 <- SamResults[-pick,]
mean(Minus30$Ld50BootRange, na.rm=T)
sd(Minus30$Ld50BootRange, na.rm=T)


pick <- which(SamResults$PctL50 < 30)
Plus30 <- SamResults[-pick,]
#Plus30 <- droplevels(Plus30)
mean(Plus30$Ld50BootRange, na.rm=T)
sd(Plus30$Ld50BootRange, na.rm=T)


##reformatting the data for bw plots
Plus30<-as.data.frame(Plus30$Ld50BootRange)
colnames(Plus30)[1] <- "L50Rci"
Plus30$Op<-">=30%"
Minus30<-as.data.frame(Minus30$Ld50BootRange)
colnames(Minus30)[1] <- "L50Rci"
Minus30$Op<-"<30%"
PctSample<-rbind(Plus30, Minus30)

#anova diferences L50CI%
boxcox(PctSample$L50Rci^-0.7~PctSample$Op)
fit<-aov(L50Rci^-0.7~Op, data=PctSample)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))

#Boxplot by subblock
ggplot(PctSample, aes(x=Op, y=L50Rci)) + 
  xlab(bquote('% Sample <'~L['50%']~'.')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
  geom_boxplot(outlier.colour = "black", outlier.size = 3)+
  theme_bw()+#white background
  theme(legend.position="none",
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))



rangeCI95<-as.data.frame(SamResults$Ld95BootRange)
colnames(rangeCI95)[1] <- "CIRange"
rangeCI95$LM<-"L95%"

rangeCI50<-as.data.frame(SamResults$Ld50BootRange)
colnames(rangeCI50)[1] <- "CIRange"
rangeCI50$LM<-"L50%"
CIRange<-rbind(rangeCI50, rangeCI95)
rm(rangeCI50, rangeCI95)




#anova diferences L%
boxcox(CIRange$CIRange^-0.3~CIRange$LM)
fit<-aov(CIRange^-0.3~LM, data=CIRange)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))
tHSDlm<- TukeyHSD(fit, ordered = FALSE, conf.level = 0.95)
tHSDlm

#set working dataframe for Tukey label function
ASM<-CIRange


#Boxplot by maturity L%
ggplot(CIRange, aes(x=LM, y=CIRange)) + 
  xlab("Maturity Estimate") + labs(y=expression(paste("C.I. Range (mm)")))+
  geom_boxplot(outlier.colour = "black", outlier.size = 3)+
  theme_bw()+#white background
  theme(legend.position="none",
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))
  #geom_text(data = generate_label_df(tHSDlm, "LM"), aes(x = plot.labels, y = 0.2, label = labels))
dev.off()
ddply(CIRange,.(LM), summarize,  M = median(CIRange, na.rm=T))


#Range of maturity onset (l95-L05) against CIraangeL50%
SamResults1$Mrange<-SamResults1$LD95-SamResults1$LD05
plot(SamResults1$Mrange, SamResults1$Ld50BootRange)
#anova diferences L%
boxcox(SamResults1$Ld50BootRange^-0.5~SamResults1$Mrange)
fit<-lm(Ld50BootRange^-0.5~Mrange, data=SamResults)
summary(fit)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))


ggplot(data = SamResults, aes(x=Mrange,  y=Ld50BootRange, color=ifelse(Ld50BootRange>25, 'red', 'black'))) + 
  geom_point()+
  xlab(bquote(~M['range']~'(mm)')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
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


#Range of IQ against CIrangeL50%
plot(SamResults$IQ, SamResults$Ld50BootRange)


#remove outliers from dataset and run lm
pick <- which(SamResults$Ld50BootRange >25)
SamResults1 <- SamResults[-pick,]
SamResults1 <- droplevels(SamResults1)

#anova diferences L%
boxcox(SamResults1$Ld50BootRange~SamResults1$IQR)
fit<-lm(Ld50BootRange~IQR, data=SamResults1)
summary(fit)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))


ggplot(data = SamResults, aes(x=IQR,  y=Ld50BootRange, color=ifelse(Ld50BootRange>25, 'red', 'black'))) + 
  geom_point()+
  xlab(bquote(~IQ['range']~'(mm)')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
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


####Histogram of CIrangeL50
ggplot(data = SamResults, aes(x=Ld50BootRange)) + 
  geom_histogram(bins = 30)+
  xlab(bquote(~CI['range']~'L'['50%']))+
  ylim(0,120)+
  theme_bw()+
  #scale_fill_identity()+ #this makes sure the color follows the color argument above in aes()
  theme(legend.position=c(0.9, 0.8))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14))+
  theme(axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))

#anova L50 by Zone
boxcox(SamResults$LD50^1.5~SamResults$Zone)
fit<-aov(LD50^1.5~Zone, data=SamResults)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))
tHSDlm<- TukeyHSD(fit, ordered = FALSE, conf.level = 0.95)
tHSDlm

#set working dataframe for Tukey label function
ASM<-SamResults


#Boxplot by maturity L%
ggplot(SamResults, aes(x=Zone, y=LD50)) + 
  xlab("Zone") +  ylab(bquote(~L['50%']~'(mm)'))+
  geom_boxplot(outlier.colour = "black", outlier.size = 3)+
  theme_bw()+#white background
  theme(legend.position="none",
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))+
geom_text(data = generate_label_df(tHSDlm, "Zone"), aes(x = plot.labels, y = 60, label = labels))


#Boxplot by MLStc by zone
ggplot(SamResults, aes(x=Zone, y=MLStc)) + 
  xlab("Zone") +  ylab("MLStc")+
  geom_boxplot(outlier.colour = "black", outlier.size = 3)+
  theme_bw()+#white background
  theme(legend.position="none",
        axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))+
  #geom_text(data = generate_label_df(tHSDlm, "Zone"), aes(x = plot.labels, y = 60, label = labels))