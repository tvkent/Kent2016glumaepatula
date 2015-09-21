library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
setwd('~/RiceResults/Thetas/')
#read in files
allo <- fread('allo_Diversity.thetas.gz.pestPG',sep='\t')
sympno276 <- fread('sympno276_Diversity.thetas.gz.pestPG',sep='\t')
ogno276 <- fread('ogno276_Diversity.thetas.gz.pestPG',sep='\t')
og <- fread('og_Diversity.thetas.gz.pestPG',sep='\t')
symp <- fread('symp_Diversity.thetas.gz.pestPG',sep='\t')
indica <- fread('indica_Diversity.thetas.gz.pestPG',sep='\t')

#add population code
allo <- allo %>% mutate(pop='allo')
sympno276 <- sympno276 %>% mutate(pop='sympno276')
ogno276 <- ogno276 %>% mutate(pop='ogno276')
og <- og %>% mutate(pop='og')
symp <- symp %>% mutate(pop='symp')
indica <- indica %>% mutate(pop='indica')

#mean values
#Wattersons
allo.watterson <- allo %>% group_by(Chr) %>% summarise(mean=mean(tW))
sympno276.watterson <- sympno276 %>% group_by(Chr) %>% summarise(mean=mean(tW))
ogno276.watterson <- ogno276 %>% group_by(Chr) %>% summarise(mean=mean(tW))
indica.watterson <- indica %>% group_by(Chr) %>% summarise(mean=mean(tW))

#TajD
allo.tajd <- allo %>% group_by(Chr) %>% summarise(mean=mean(Tajima))
sympno276.tajd <- sympno276 %>% group_by(Chr) %>% summarise(mean=mean(Tajima))
ogno276.tajd <- ogno276 %>% group_by(Chr) %>% summarise(mean=mean(Tajima))
indica.tajd <- indica %>% group_by(Chr) %>% summarise(mean=mean(Tajima))

#split by chr
#Chr 1
allo.1 <- subset(allo, Chr==1)
sympno276.1 <- subset(sympno276, Chr==1)
ogno276.1 <- subset(ogno276, Chr==1)
og.1 <- subset(og, Chr==1)
symp.1 <- subset(symp, Chr==1)
indica.1 <- subset(indica, Chr==1)
#Chr 2
allo.2 <- subset(allo, Chr==2)
sympno276.2 <- subset(sympno276, Chr==2)
ogno276.2 <- subset(ogno276, Chr==2)
og.2 <- subset(og, Chr==2)
symp.2 <- subset(symp, Chr==2)
indica.2 <- subset(indica, Chr==2)
#Chr 3
allo.3 <- subset(allo, Chr==3)
sympno276.3 <- subset(sympno276, Chr==3)
ogno276.3 <- subset(ogno276, Chr==3)
og.3 <- subset(og, Chr==3)
symp.3 <- subset(symp, Chr==3)
indica.3 <- subset(indica, Chr==3)
#Chr 4
allo.4 <- subset(allo, Chr==4)
sympno276.4 <- subset(sympno276, Chr==4)
ogno276.4 <- subset(ogno276, Chr==4)
og.4 <- subset(og, Chr==4)
symp.4 <- subset(symp, Chr==4)
indica.4 <- subset(indica, Chr==4)
#Chr 5
allo.5 <- subset(allo, Chr==5)
sympno276.5 <- subset(sympno276, Chr==5)
ogno276.5 <- subset(ogno276, Chr==5)
og.5 <- subset(og, Chr==5)
symp.5 <- subset(symp, Chr==5)
indica.5 <- subset(indica, Chr==5)
#Chr 6
allo.6 <- subset(allo, Chr==6)
sympno276.6 <- subset(sympno276, Chr==6)
ogno276.6 <- subset(ogno276, Chr==6)
og.6 <- subset(og, Chr==6)
symp.6 <- subset(symp, Chr==6)
indica.6 <- subset(indica, Chr==6)
#Chr 7
allo.7 <- subset(allo, Chr==7)
sympno276.7 <- subset(sympno276, Chr==7)
ogno276.7 <- subset(ogno276, Chr==7)
og.7 <- subset(og, Chr==7)
symp.7 <- subset(symp, Chr==7)
indica.7 <- subset(indica, Chr==7)
#Chr 8
allo.8 <- subset(allo, Chr==8)
sympno276.8 <- subset(sympno276, Chr==8)
ogno276.8 <- subset(ogno276, Chr==8)
og.8 <- subset(og, Chr==8)
symp.8 <- subset(symp, Chr==8)
indica.8 <- subset(indica, Chr==8)
#Chr 9
allo.9 <- subset(allo, Chr==9)
sympno276.9 <- subset(sympno276, Chr==9)
ogno276.9 <- subset(ogno276, Chr==9)
og.9 <- subset(og, Chr==9)
symp.9 <- subset(symp, Chr==9)
indica.9 <- subset(indica, Chr==9)
#Chr 10
allo.10 <- subset(allo, Chr==10)
sympno276.10 <- subset(sympno276, Chr==10)
ogno276.10 <- subset(ogno276, Chr==10)
og.10 <- subset(og, Chr==10)
symp.10 <- subset(symp, Chr==10)
indica.10 <- subset(indica, Chr==10)
#Chr 11
allo.11 <- subset(allo, Chr==11)
sympno276.11 <- subset(sympno276, Chr==11)
ogno276.11 <- subset(ogno276, Chr==11)
og.11 <- subset(og, Chr==11)
symp.11 <- subset(symp, Chr==11)
indica.11 <- subset(indica, Chr==11)
#Chr 12
allo.12 <- subset(allo, Chr==12)
sympno276.12 <- subset(sympno276, Chr==12)
ogno276.12 <- subset(ogno276, Chr==12)
og.12 <- subset(og, Chr==12)
symp.12 <- subset(symp, Chr==12)
indica.12 <- subset(indica, Chr==12)

#drop sites with tW & tP == 0 for tajima plot
#Chr 1
allo.1.taj <- subset(allo.1, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.1.taj <- subset(sympno276.1, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.1.taj <- subset(ogno276.1, tW!=0 & tP!=0) %>% subset(nSites>150)
og.1.taj <- subset(og.1, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.1.taj <- subset(symp.1, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.1.taj <- subset(indica.1, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 2
allo.2.taj <- subset(allo.2, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.2.taj <- subset(sympno276.2, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.2.taj <- subset(ogno276.2, tW!=0 & tP!=0) %>% subset(nSites>150)
og.2.taj <- subset(og.2, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.2.taj <- subset(symp.2, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.2.taj <- subset(indica.2, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 3
allo.3.taj <- subset(allo.3, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.3.taj <- subset(sympno276.3, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.3.taj <- subset(ogno276.3, tW!=0 & tP!=0) %>% subset(nSites>150)
og.3.taj <- subset(og.3, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.3.taj <- subset(symp.3, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.3.taj <- subset(indica.3, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 4
allo.4.taj <- subset(allo.4, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.4.taj <- subset(sympno276.4, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.4.taj <- subset(ogno276.4, tW!=0 & tP!=0) %>% subset(nSites>150)
og.4.taj <- subset(og.4, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.4.taj <- subset(symp.4, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.4.taj <- subset(indica.4, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 5
allo.5.taj <- subset(allo.5, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.5.taj <- subset(sympno276.5, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.5.taj <- subset(ogno276.5, tW!=0 & tP!=0) %>% subset(nSites>150)
og.5.taj <- subset(og.5, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.5.taj <- subset(symp.5, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.5.taj <- subset(indica.5, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 6
allo.6.taj <- subset(allo.6, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.6.taj <- subset(sympno276.6, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.6.taj <- subset(ogno276.6, tW!=0 & tP!=0) %>% subset(nSites>150)
og.6.taj <- subset(og.6, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.6.taj <- subset(symp.6, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.6.taj <- subset(indica.6, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 7
allo.7.taj <- subset(allo.7, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.7.taj <- subset(sympno276.7, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.7.taj <- subset(ogno276.7, tW!=0 & tP!=0) %>% subset(nSites>150)
og.7.taj <- subset(og.7, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.7.taj <- subset(symp.7, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.7.taj <- subset(indica.7, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 8
allo.8.taj <- subset(allo.8, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.8.taj <- subset(sympno276.8, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.8.taj <- subset(ogno276.8, tW!=0 & tP!=0) %>% subset(nSites>150)
og.8.taj <- subset(og.8, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.8.taj <- subset(symp.8, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.8.taj <- subset(indica.8, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 9
allo.9.taj <- subset(allo.9, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.9.taj <- subset(sympno276.9, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.9.taj <- subset(ogno276.9, tW!=0 & tP!=0) %>% subset(nSites>150)
og.9.taj <- subset(og.9, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.9.taj <- subset(symp.9, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.9.taj <- subset(indica.9, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 10
allo.10.taj <- subset(allo.10, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.10.taj <- subset(sympno276.10, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.10.taj <- subset(ogno276.10, tW!=0 & tP!=0) %>% subset(nSites>150)
og.10.taj <- subset(og.10, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.10.taj <- subset(symp.10, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.10.taj <- subset(indica.10, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 11
allo.11.taj <- subset(allo.11, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.11.taj <- subset(sympno276.11, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.11.taj <- subset(ogno276.11, tW!=0 & tP!=0) %>% subset(nSites>150)
og.11.taj <- subset(og.11, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.11.taj <- subset(symp.11, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.11.taj <- subset(indica.11, tW!=0 & tP!=0) %>% subset(nSites>150)
#Chr 12
allo.12.taj <- subset(allo.12, tW!=0 & tP!=0) %>% subset(nSites>150)
sympno276.12.taj <- subset(sympno276.12, tW!=0 & tP!=0) %>% subset(nSites>150)
ogno276.12.taj <- subset(ogno276.12, tW!=0 & tP!=0) %>% subset(nSites>150)
og.12.taj <- subset(og.12, tW!=0 & tP!=0) %>% subset(nSites>150)
symp.12.taj <- subset(symp.12, tW!=0 & tP!=0) %>% subset(nSites>150)
indica.12.taj <- subset(indica.12, tW!=0 & tP!=0) %>% subset(nSites>150)

#merge and get spearman

#plot shit
#Watterson's
#Chr 1
ggplot(allo.1, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.1, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.1, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.1, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 1")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black')) 
#Chr 2
ggplot(allo.2, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.2, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.2, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.2, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 2")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 3
ggplot(allo.3, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.3, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.3, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.3, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 3")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 4
ggplot(allo.4, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.4, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.4, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.4, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 4")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 5
ggplot(allo.5, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.5, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.5, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.5, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 5")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 6
ggplot(allo.6, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.6, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.6, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.6, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 6")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 7
ggplot(allo.7, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.7, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.7, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.7, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 7")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 8
ggplot(allo.8, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.8, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.8, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.8, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 8")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 9
ggplot(allo.9, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.9, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.9, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.9, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 9")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 10
ggplot(allo.10, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.10, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.10, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.10, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 10")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 11
ggplot(allo.11, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.11, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.11, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.11, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 11")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 12
ggplot(allo.12, aes(x=WinCenter, y=tW))+
  stat_smooth(data = allo.12, aes(x=WinCenter, y=tW, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.12, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.12, aes(x=WinCenter, y=tW, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Watterson's Theta along chromosome 12")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))

#TajD
#Chr 1
ggplot(allo.1.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.1.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.1.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.1.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 1")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black')) 
#Chr 2
ggplot(allo.2.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.2.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.2.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.2.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 2")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 3
ggplot(allo.3.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.3.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.3.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.3.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 3")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 4
ggplot(allo.4.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.4.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.4.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.4.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 4")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 5
ggplot(allo.5.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.5.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.5.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.5.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D Theta along chromosome 5")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 6
ggplot(allo.6.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.6.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.6.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.6.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D Theta along chromosome 6")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 7
ggplot(allo.7.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.7.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.7.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.7.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 7")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 8
ggplot(allo.8.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.8.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.8.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.8.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 8")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 9
ggplot(allo.9.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.9.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.9.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.9.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 9")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 10
ggplot(allo.10.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.10.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.10.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.10.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 10")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 11
ggplot(allo.11.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.11.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.11.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.11.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 11")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 12
ggplot(allo.12.taj, aes(x=WinCenter, y=Tajima))+
  stat_smooth(data = allo.12.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.12.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = indica.12.taj, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  theme_bw()+
  ggtitle("Tajima's D along chromosome 12")+
  scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))

#tajima histograms
#allo
hist(allo.1.taj$Tajima)
hist(allo.2.taj$Tajima)
hist(allo.3.taj$Tajima)
hist(allo.4.taj$Tajima)
hist(allo.5.taj$Tajima)
hist(allo.6.taj$Tajima)
hist(allo.7.taj$Tajima)
hist(allo.8.taj$Tajima)
hist(allo.9.taj$Tajima)
hist(allo.10.taj$Tajima)
hist(allo.11.taj$Tajima)
hist(allo.12.taj$Tajima)
#sympno276
hist(sympno276.1.taj$Tajima)
hist(sympno276.2.taj$Tajima)
hist(sympno276.3.taj$Tajima)
hist(sympno276.4.taj$Tajima)
hist(sympno276.5.taj$Tajima)
hist(sympno276.6.taj$Tajima)
hist(sympno276.7.taj$Tajima)
hist(sympno276.8.taj$Tajima)
hist(sympno276.9.taj$Tajima)
hist(sympno276.10.taj$Tajima)
hist(sympno276.11.taj$Tajima)
hist(sympno276.12.taj$Tajima)
#indica
hist(indica.1.taj$Tajima)
hist(indica.2.taj$Tajima)
hist(indica.3.taj$Tajima)
hist(indica.4.taj$Tajima)
hist(indica.5.taj$Tajima)
hist(indica.6.taj$Tajima)
hist(indica.7.taj$Tajima)
hist(indica.8.taj$Tajima)
hist(indica.9.taj$Tajima)
hist(indica.10.taj$Tajima)
hist(indica.11.taj$Tajima)
hist(indica.12.taj$Tajima)

    #Pi
    #Chr 1
    ggplot(allo.1, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.1, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.1, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.1, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 1")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black')) 
    #Chr 2
    ggplot(allo.2, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.2, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.2, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.2, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 2")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 3
    ggplot(allo.3, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.3, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.3, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.3, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 3")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 4
    ggplot(allo.4, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.4, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.4, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.4, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 4")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 5
    ggplot(allo.5, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.5, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.5, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.5, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 5")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 6
    ggplot(allo.6, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.6, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.6, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.6, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 6")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 7
    ggplot(allo.7, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.7, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.7, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.7, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 7")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 8
    ggplot(allo.8, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.8, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.8, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.8, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 8")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 9
    ggplot(allo.9, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.9, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.9, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.9, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 9")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 10
    ggplot(allo.10, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.10, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.10, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.10, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 10")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 11
    ggplot(allo.11, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.11, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.11, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.11, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 11")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))
    #Chr 12
    ggplot(allo.12, aes(x=WinCenter, y=tP/nSites))+
      stat_smooth(data = allo.12, aes(x=WinCenter, y=tP/nSites, colour=pop, fill = pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = sympno276.12, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      stat_smooth(data = indica.12, aes(x=WinCenter, y=tP/nSites, colour=pop, fill=pop),
                  method="gam", formula=y~s(x,k=80))+
      theme_bw()+
      ggtitle(expression(paste(pi, " along chromosome 12")))+
      scale_colour_manual(name='Population', values=c('darkred','darkblue','darkgreen'),
                          guide = guide_legend(override.aes=aes(fill=NA)))+
      scale_fill_manual(values = c('red','blue','green'),guide=FALSE)+
      theme(legend.key = element_rect(colour = 'white'),
            legend.background=element_rect(colour='black'))

#Zeng's E doesn't work
#Chr 1
ggplot(ogno276.1, aes(x=WinCenter, y=zeng))+
  stat_smooth(data = ogno276.1, aes(x=WinCenter, y=zeng, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=40))+
  theme_bw()+
  ggtitle("Zeng's E along chromosome 1")+
  scale_colour_manual(name='Population', values='orange',
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = 'yellow',guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black')) 
#Chr 2
ggplot(allo.2, aes(x=WinCenter, y=Tajima)+
         stat_smooth(data = allo.2, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.2, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 2")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 3
ggplot(allo.3, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.3, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.3, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 3")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 4
ggplot(allo.4, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.4, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.4, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 4")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 5
ggplot(allo.5, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.5, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.5, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D Theta along chromosome 5")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 6
ggplot(allo.6, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.6, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.6, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D Theta along chromosome 6")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 7
ggplot(allo.7, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.7, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.7, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 7")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 8
ggplot(allo.8, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.8, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.8, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 8")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 9
ggplot(allo.9, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.9, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.9, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 9")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 10
ggplot(allo.10, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.10, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.10, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 10")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 11
ggplot(allo.11, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.11, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.11, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Tajima's D along chromosome 11")+
         scale_colour_manual(name='Population', values=c('darkred','darkblue'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('red','blue'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))
#Chr 12
ggplot(allo.12, aes(x=WinCenter, y=Tajima))+
         stat_smooth(data = allo.12, aes(x=WinCenter, y=Tajima, colour=pop, fill = pop),
                     method="gam", formula=y~s(x,k=80))+
         stat_smooth(data = sympno276.12, aes(x=WinCenter, y=Tajima, colour=pop, fill=pop),
                     method="gam", formula=y~s(x,k=80))+
         theme_bw()+
         ggtitle("Zheng's E along chromosome 12")+
         scale_colour_manual(name='Population', values=c('orange'),
                             guide = guide_legend(override.aes=aes(fill=NA)))+
         scale_fill_manual(values = c('yellow'),guide=FALSE)+
         theme(legend.key = element_rect(colour = 'white'),
               legend.background=element_rect(colour='black'))