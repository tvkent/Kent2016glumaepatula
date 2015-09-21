library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
setwd('~/RiceResults/Fst/')
#set column names
column.names <- c('A','AB','f','FST','Pvar','Chr','bp')

#read in files and set headers
allo.sympno276 <- fread('allo.sympno276.fst')
intersect <- fread('intersect.allo.sympno276_intergenic.txt')
allo.sympno276 <- cbind(allo.sympno276,intersect)
setnames(allo.sympno276,column.names)

allo.indica <- fread('allo.indica.fst')
intersect <- fread('intersect.allo.indica_intergenic.txt')
allo.indica <- cbind(allo.indica,intersect)
setnames(allo.indica,column.names)

ogno276.indica <- fread('ogno276.indica.fst')
intersect <- fread('intersect.ogno276.indica_intergenic.txt')
ogno276.indica <- cbind(ogno276.indica,intersect)
setnames(ogno276.indica,column.names)

sympno276.indica <- fread('sympno276.indica.fst')
intersect <- fread('intersect.sympno276.indica_intergenic.txt')
sympno276.indica <- cbind(sympno276.indica,intersect)
setnames(sympno276.indica,column.names)

#filter out invariant remnants and keep only Chr, bp, & FST columns
allo.sympno276 <- allo.sympno276 %>% mutate(pop='allo.sympno276') %>% select(Chr, bp, FST, pop) %>% filter(FST>=0 & FST<=1)
allo.indica <- allo.indica %>% mutate(pop='allo.indica') %>% select(Chr, bp, FST, pop) %>% filter(FST>=0 & FST<=1)
ogno276.indica <- ogno276.indica %>% mutate(pop='ogno276.indica') %>% select(Chr, bp, FST, pop) %>% filter(FST>=0 & FST<=1)
sympno276.indica <- sympno276.indica %>% mutate(pop='sympno276.indica') %>% select(Chr, bp, FST, pop) %>% filter(FST>=0 & FST<=1)

#subset out top 15%
allo.sympno276.15 <- subset(allo.sympno276, FST >= quantile(FST, .98))
allo.indica.15 <- subset(allo.indica, FST >= quantile(FST, .95))
ogno276.indica.15 <- subset(ogno276.indica, FST >= quantile(FST, .95))
sympno276.indica.15 <- subset(sympno276.indica, FST >= quantile(FST, .95))

#format as BED and write to file
allo.sympno276.bed <- allo.sympno276.15 %>% mutate(start=bp-1) %>% select(Chr, start, bp)
write.table(allo.sympno276.bed, quote= FALSE, file='allo.sympno276.bed', sep='\t', col.names=FALSE, row.names=FALSE)
allo.indica.bed <- allo.indica.15 %>% mutate(start=bp-1) %>% select(Chr, start, bp)
write.table(allo.indica.bed, quote= FALSE, file='allo.indica.bed', sep='\t', col.names=FALSE, row.names=FALSE)
sympno276.indica.bed <- sympno276.indica.15 %>% mutate(start=bp-1) %>% select(Chr, start, bp)
write.table(sympno276.indica.bed, quote= FALSE, file='sympno276.indica.bed', sep='\t', col.names=FALSE, row.names=FALSE)
ogno276.indica.bed <- ogno276.indica.15 %>% mutate(start=bp-1) %>% select(Chr, start, bp)
write.table(ogno276.indica.bed, quote= FALSE, file='ogno276.indica.bed', sep='\t', col.names=FALSE, row.names=FALSE)

#summary stats of full set and by chr
allo.sympno276.sum <- allo.sympno276 %>% group_by(Chr) %>% summarise(mean=mean(FST))
allo.sympno276.totsum <- summary(allo.sympno276$FST)
allo.indica.totsum <- summary(allo.indica$FST)
allo.indica.sum <- allo.indica %>% group_by(Chr) %>% summarise(mean=mean(FST))
ogno276.indica.totsum <- summary(ogno276.indica$FST)
ogno276.indica.sum <- ogno276.indica %>% group_by(Chr) %>% summarise(mean=mean(FST))
sympno276.indica.totsum <- summary(sympno276.indica$FST)
sympno276.indica.sum <- sympno276.indica %>% group_by(Chr) %>% summarise(mean=mean(FST))

#subset by chr
#Chr 1
allo.indica.1 <- subset(allo.indica, Chr==1)
sympno276.indica.1 <- subset(sympno276.indica, Chr==1)
ogno276.indica.1 <- subset(ogno276.indica, Chr==1)
allo.sympno276.1 <- subset(allo.sympno276, Chr==1)
#Chr 2
allo.indica.2 <- subset(allo.indica, Chr==2)
sympno276.indica.2 <- subset(sympno276.indica, Chr==2)
ogno276.indica.2 <- subset(ogno276.indica, Chr==2)
allo.sympno276.2 <- subset(allo.sympno276, Chr==2)
#Chr 3
allo.indica.3 <- subset(allo.indica, Chr==3)
sympno276.indica.3 <- subset(sympno276.indica, Chr==3)
ogno276.indica.3 <- subset(ogno276.indica, Chr==3)
allo.sympno276.3 <- subset(allo.sympno276, Chr==3)
#Chr 4
allo.indica.4 <- subset(allo.indica, Chr==4)
sympno276.indica.4 <- subset(sympno276.indica, Chr==4)
ogno276.indica.4 <- subset(ogno276.indica, Chr==4)
allo.sympno276.4 <- subset(allo.sympno276, Chr==4)
#Chr 5
allo.indica.5 <- subset(allo.indica, Chr==5)
sympno276.indica.5 <- subset(sympno276.indica, Chr==5)
ogno276.indica.5 <- subset(ogno276.indica, Chr==5)
allo.sympno276.5 <- subset(allo.sympno276, Chr==5)
#Chr 6
allo.indica.6 <- subset(allo.indica, Chr==6)
sympno276.indica.6 <- subset(sympno276.indica, Chr==6)
ogno276.indica.6 <- subset(ogno276.indica, Chr==6)
allo.sympno276.6 <- subset(allo.sympno276, Chr==6)
#Chr 7
allo.indica.7 <- subset(allo.indica, Chr==7)
sympno276.indica.7 <- subset(sympno276.indica, Chr==7)
ogno276.indica.7 <- subset(ogno276.indica, Chr==7)
allo.sympno276.7 <- subset(allo.sympno276, Chr==7)
#Chr 8
allo.indica.8 <- subset(allo.indica, Chr==8)
sympno276.indica.8 <- subset(sympno276.indica, Chr==8)
ogno276.indica.8 <- subset(ogno276.indica, Chr==8)
allo.sympno276.8 <- subset(allo.sympno276, Chr==8)
#Chr 9
allo.indica.9 <- subset(allo.indica, Chr==9)
sympno276.indica.9 <- subset(sympno276.indica, Chr==9)
ogno276.indica.9 <- subset(ogno276.indica, Chr==9)
allo.sympno276.9 <- subset(allo.sympno276, Chr==9)
#Chr 10
allo.indica.10 <- subset(allo.indica, Chr==10)
sympno276.indica.10 <- subset(sympno276.indica, Chr==10)
ogno276.indica.10 <- subset(ogno276.indica, Chr==10)
allo.sympno276.10 <- subset(allo.sympno276, Chr==10)
#Chr 11
allo.indica.11 <- subset(allo.indica, Chr==11)
sympno276.indica.11 <- subset(sympno276.indica, Chr==11)
ogno276.indica.11 <- subset(ogno276.indica, Chr==11)
allo.sympno276.11 <- subset(allo.sympno276, Chr==11)
#Chr 12
allo.indica.12 <- subset(allo.indica, Chr==12)
sympno276.indica.12 <- subset(sympno276.indica, Chr==12)
ogno276.indica.12 <- subset(ogno276.indica, Chr==12)
allo.sympno276.12 <- subset(allo.sympno276, Chr==12)

#merge and get spearman
chr1 <- merge(allo.indica.1,sympno276.indica.1,by='bp')
cor.test(chr1$FST.x, chr1$FST.y, method='spearman')
chr2 <- merge(allo.indica.2,sympno276.indica.2,by='bp')
cor.test(chr2$FST.x, chr2$FST.y, method='spearman')
chr3 <- merge(allo.indica.3,sympno276.indica.3,by='bp')
cor.test(chr3$FST.x, chr3$FST.y, method='spearman')
chr4 <- merge(allo.indica.4,sympno276.indica.4,by='bp')
cor.test(chr4$FST.x, chr4$FST.y, method='spearman')
chr5 <- merge(allo.indica.5,sympno276.indica.5,by='bp')
cor.test(chr5$FST.x, chr5$FST.y, method='spearman')
chr6 <- merge(allo.indica.6,sympno276.indica.6,by='bp')
cor.test(chr6$FST.x, chr6$FST.y, method='spearman')
chr7 <- merge(allo.indica.7,sympno276.indica.7,by='bp')
cor.test(chr7$FST.x, chr7$FST.y, method='spearman')
chr8 <- merge(allo.indica.8,sympno276.indica.8,by='bp')
cor.test(chr8$FST.x, chr8$FST.y, method='spearman')
chr9 <- merge(allo.indica.9,sympno276.indica.9,by='bp')
cor.test(chr9$FST.x, chr9$FST.y, method='spearman')
chr10 <- merge(allo.indica.10,sympno276.indica.10,by='bp')
cor.test(chr10$FST.x, chr10$FST.y, method='spearman')
chr11 <- merge(allo.indica.11,sympno276.indica.11,by='bp')
cor.test(chr11$FST.x, chr11$FST.y, method='spearman')
chr12 <- merge(allo.indica.12,sympno276.indica.12,by='bp')
cor.test(chr12$FST.x, chr12$FST.y, method='spearman')
#Plot full set by chr
#Chr 1
ggplot(allo.indica.1, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.1, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.1, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.1, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 1')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black')) 
#Chr 2
ggplot(allo.indica.2, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.2, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.2, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.2, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 2')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 3
ggplot(allo.indica.3, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.3, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.3, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.3, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 3')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 4
ggplot(allo.indica.4, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.4, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.4, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.4, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 4')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 5
ggplot(allo.indica.5, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.5, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.5, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.5, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 5')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 6
ggplot(allo.indica.6, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.6, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.6, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.6, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 6')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 7
ggplot(allo.indica.7, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.7, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.7, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.7, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 7')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 8
ggplot(allo.indica.8, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.8, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.8, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.8, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 8')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 9
ggplot(allo.indica.9, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.9, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.9, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.9, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 9')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 10
ggplot(allo.indica.10, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.10, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.10, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.10, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 10')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 11
ggplot(allo.indica.11, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.11, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.11, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.11, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 11')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
#Chr 12
ggplot(allo.indica.12, aes(x=bp, y=FST))+
  stat_smooth(data = allo.indica.12, aes(x=bp, y=FST, colour=pop, fill = pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = sympno276.indica.12, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  stat_smooth(data = allo.sympno276.12, aes(x=bp, y=FST, colour=pop, fill=pop),
              method="gam", formula=y~s(x,k=80))+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  ggtitle('FST along chromosome 12')+
  scale_colour_manual(name='Population', values=c('darkred','darkgreen','darkblue'),
                      guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_fill_manual(values = c('red','green','blue'),guide=FALSE)+
  theme(legend.key = element_rect(colour = 'white'),
        legend.background=element_rect(colour='black'))
