
####RESULTS ANALYSIS
##
pick <- which(SamResults$Ld50BootRange>=9.4)
picked <- SamResults[-pick,]
picked<-droplevels(picked)

ggplot(data = picked, aes(x=PctL50,  y=Ld50BootRange, color=ifelse(Ld50BootRange>9.5, 'red', 'black')))+
  geom_point()+
  geom_point(data = SamResults, aes(x=PctL50,  y=Ld50BootRange, color=ifelse(Ld50BootRange>9.5, 'red', 'black')))+
  xlab(bquote('% Sample <'~L['50%']~'.')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
  #geom_smooth(method=lm, se=F, fill='Black', fullrange=F, size=1.2, color='black')+
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

fitd <- lm(Ld50BootRange~PctL50, data=SamResults)
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


# ##reformatting the data for bw plots
# Plus30<-as.data.frame(Plus30$Ld50BootRange)
# colnames(Plus30)[1] <- "L50Rci"
# Plus30$Op<-">=30%"
# Minus30<-as.data.frame(Minus30$Ld50BootRange)
# colnames(Minus30)[1] <- "L50Rci"
# Minus30$Op<-"<30%"
# PctSample<-rbind(Plus30, Minus30)
# 
# #anova diferences L50CI%
# boxcox(PctSample$L50Rci^-0.7~PctSample$Op)
# fit<-aov(L50Rci^-0.7~Op, data=PctSample)
# anova(fit)
# par(mfrow = c(2,2))
# plot(fit)
# par(mfrow = c(1,1))
# 
# #Boxplot by subblock
# ggplot(PctSample, aes(x=Op, y=L50Rci)) + 
#   xlab(bquote('% Sample <'~L['50%']~'.')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
#   geom_boxplot(outlier.colour = "black", outlier.size = 3)+
#   theme_bw()+#white background
#   theme(legend.position="none",
#         axis.title.x = element_text(size=14),
#         axis.text.x  = element_text(size=14),
#         axis.title.y = element_text(size=14),
#         axis.text.y  = element_text(size=14))



rangeCI95<-as.data.frame(SamResults$Ld95BootRange)
colnames(rangeCI95)[1] <- "CIRange"
rangeCI95$LM<-"L95%"

rangeCI50<-as.data.frame(SamResults$Ld50BootRange)
colnames(rangeCI50)[1] <- "CIRange"
rangeCI50$LM<-"L50%"
CIRange<-rbind(rangeCI50, rangeCI95)
rm(rangeCI50, rangeCI95)




#anova diferences L%
boxcox(log(CIRange$CIRange)~CIRange$LM)
fit<-aov(log(CIRange)~LM, data=CIRange)
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

# 
# #remove outliers from dataset and run lm
# pick <- which(SamResults$Ld50BootRange >25)
# SamResults1 <- SamResults[-pick,]
# SamResults1 <- droplevels(SamResults1)

#anova diferences L%
boxcox(SamResults$Ld50BootRange^-0.6~SamResults$IQR)
fit<-lm(Ld50BootRange^-0.6~IQR, data=SamResults)
summary(fit)
anova(fit)
par(mfrow = c(2,2))
plot(fit)
par(mfrow = c(1,1))


ggplot(data = SamResults, aes(x=IQR,  y=Ld50BootRange, color=ifelse(Ld50BootRange>12, 'red', 'black'))) + 
  geom_point()+
  xlab(bquote(~IQ['range']~'(mm)')) + ylab(bquote('C.I. Range'~L['50%']~'(mm)'))+
  #geom_smooth(method=lm, se=F, fill='Black', fullrange=F, size=1.2, color='black')+
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
  #ylim(0,120)+
  theme_bw()+
  #scale_fill_identity()+ #this makes sure the color follows the color argument above in aes()
  theme(legend.position=c(0.9, 0.8))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14))+
  theme(axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))

#% below thresholds
pick <- which(SamResults$Ld50BootRange <10)
ci5less <- SamResults[pick,]
ci5less <- droplevels(ci5less)

dim(ci5less)
296/305*100
max(ci5less$b)

#anova L50 by Zone
boxcox(SamResults$LD50^3.6~SamResults$Zone)
fit<-aov(LD50^3.6~Zone, data=SamResults)
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




plot(SamResults$b, SamResults$Ld50BootRange)

ggplot(data = SamResults, aes(x=b)) + 
  geom_histogram(bins = 20)+
  xlab("Slope")+
  #ylim(0,120)+
  theme_bw()+
  #scale_fill_identity()+ #this makes sure the color follows the color argument above in aes()
  theme(legend.position=c(0.9, 0.8))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=14))+
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14))+
  theme(axis.title.y = element_text(size=14),
        axis.text.y  = element_text(size=14))
#
#########################
#           Plots for eLML L50%
#########################
#
#
doPlot = function(LFPlot) {
  dum = subset(BlockSumStats, Zone == LFPlot)
  ggobj = ggplot(data = dum, aes(y=mn.eLML, x=as.factor(BlockNo))) + 
    xlab("BlockNo") +
    ylab("eLML (mm)") +
    labs(title= dum$Zone, size=10)+
    #ylim(min(dum$sd.eLML-10), max(dum$sd.eLML+10))+
    geom_point(position=position_dodge(), stat="identity", size =3) +
    geom_errorbar(aes(ymin=dum$mn.eLMLbootL95, ymax=dum$mn.eLMLbootU95),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=14),
          legend.position="none")
  #ggsave(sprintf("%s_LFplot.tiff", LFPlot))
  print(ggobj)
}
lapply(unique(BlockSumStats$Zone), doPlot)

doPlot = function(LFPlot) {
  dum = subset(SamResults, Zone == LFPlot)
  ggobj = ggplot(data = dum, aes(y=eLML, x=as.factor(BlockNo))) + 
    xlab("BlockNo") +
    ylab("eLML (mm)") +
    labs(title= dum$Zone, size=10)+
    #ylim(min(dum$sd.eLML-10), max(dum$sd.eLML+10))+
    geom_boxplot(outlier.colour = "black", outlier.size = 3)+
        theme_bw()+
    theme(legend.title=element_blank(),
          legend.text = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=14),
          legend.position="none")
  #ggsave(sprintf("%s_LFplot.tiff", LFPlot))
  print(ggobj)
}
lapply(unique(BlockSumStats$Zone), doPlot)


#eLML boxplot from the BootU95 of L50%
doPlot = function(LFPlot) {
  dum = subset(SamResults, Zone == LFPlot)
  ggobj = ggplot(data = dum, aes(y=eLMLbootU95, x=as.factor(BlockNo))) + 
    xlab("BlockNo") +
    ylab("eLML (UCI) (mm)") +
    labs(title= dum$Zone, size=10)+
    #ylim(min(dum$sd.eLML-10), max(dum$sd.eLML+10))+
    geom_boxplot(outlier.colour = "black", outlier.size = 3)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=14),
          legend.position="none")
  #ggsave(sprintf("%s_LFplot.tiff", LFPlot))
  print(ggobj)
}
lapply(unique(BlockSumStats$Zone), doPlot)


########################################
#                   eLML L95%
#######################################

doPlot = function(LFPlot) {
  dum = subset(SamResults, Zone == LFPlot)
  ggobj = ggplot(data = dum, aes(y=L95eLML, x=as.factor(BlockNo))) + 
    xlab("BlockNo") +
    ylab("L95% eLML (mm)") +
    labs(title= dum$Zone, size=10)+
    #ylim(min(dum$sd.eLML-10), max(dum$sd.eLML+10))+
    geom_boxplot(outlier.colour = "black", outlier.size = 3)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=14),
          legend.position="none")
  #ggsave(sprintf("%s_LFplot.tiff", LFPlot))
  print(ggobj)
}
lapply(unique(BlockSumStats$Zone), doPlot)


#eLML boxplot from the BootU95 of L50%
doPlot = function(LFPlot) {
  dum = subset(SamResults, Zone == LFPlot)
  ggobj = ggplot(data = dum, aes(y=eLMLbootU95, x=as.factor(BlockNo))) + 
    xlab("BlockNo") +
    ylab("eLML (UCI) (mm)") +
    labs(title= dum$Zone, size=10)+
    #ylim(min(dum$sd.eLML-10), max(dum$sd.eLML+10))+
    geom_boxplot(outlier.colour = "black", outlier.size = 3)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y  = element_text(size=14),
          legend.position="none")
  #ggsave(sprintf("%s_LFplot.tiff", LFPlot))
  print(ggobj)
}
lapply(unique(BlockSumStats$Zone), doPlot)