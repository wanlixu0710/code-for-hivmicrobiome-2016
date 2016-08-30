setwd("~/R_directory/Wei Jiang/")
library(vegan)
library(ecodist)
library(ggplot2)
library(plyr)
library(dplyr)
library(scatterplot3d)
library(lattice)
library(reshape2)
library(RSvgDevice)
library(reshape)
library(indicspecies)
library(heatmap3)
library(devEMF)
#rm(list=ls())

#### first, lets read in data clinical and bacteria
otu <- read.csv(file="120215JW515F-pr.fasta.otus.fa.OTU.percentages.csv", sep=",", header=T) # bacteria species level
rownames(otu) <- otu$Sample_ID
otu <- otu[,-1]
neg <- otu[(35:36),]
negmean <- colMeans(neg)
negmean <- data.frame(negmean)
negmean <- t(negmean)
otu2 <- otu-negmean
otu <- otu2[-(35:36),]
otu[otu < 0] <- 0  # 2490 otus

# calculate simpson diversity
otu$simpson <- diversity(otu,"simp")
otu$shannon <- diversity(otu, index = "shannon")
otu$Sample_ID <- rownames(otu)
a.div <- otu[,c("Sample_ID","simpson","shannon")]

clinical <- read.csv(file="sample information-microbiome analysis.csv", sep=",", header=T)
clinical$bcell.dif <- clinical$bcell.d7-clinical$bcell.d0
clinical$cd3.dif <- clinical$cd3.d7-clinical$cd3.d0
clinical$cd4.dif <- clinical$cd4.d7-clinical$cd4.d0

a.div <- inner_join(a.div, clinical)

L6 <- read.csv(file="120215JW515F-pr.fasta.otus.fa.genus.percentages.csv", sep=",", header=T)
rownames(L6) <- L6$Sample_ID
L6 <- L6[,-1]
neg <- L6[(35:36),]
negmean <- colMeans(neg)
negmean <- data.frame(negmean)
negmean <- t(negmean)
L62 <- L6-negmean
L6 <- L62[-(35:36),]
L6[L6 < 0] <- 0
L6$Sample_ID <- rownames(L6)
a <- a.div[,c("Sample_ID","group", "gender", "Race")]
L6.group <- inner_join(a, L6)

#############-------This is only focus on pt group------##############
pt.a.div <- filter(a.div, group=="Patient")
pt.otu <- filter(otu, Sample_ID %in% pt.a.div$Sample_ID)
pt.L6 <- filter(L6.group, group=="Patient")

#### 1. will use the titer increase as the outcome varible, and microbiome (a.div in this case) as predictor, firstly focusing on flu induced titer
# change the formating to long format in the excel and read in R
### lets focus on flu.igg, flu.igm, flu.iga first, and also bcell, cd3, cd4, h1n1, b, h3n2, these variables only have day 0 and day 14.
pt.a.div.long <- read.csv(file="pt.a.div.long.csv")
library(tidyr)
pt.a.div.w <- spread(pt.a.div.long, variable, value)

### 1.1 try with flu.igm using repeated measures of anova
summary(with(pt.a.div.w, aov(flu.igm~time+simpson))) # it doesn't seem that you could model the intercept of the vaccine titer.

### 1.2 choose the variables and run the correlations, HLM doesn't work for this type of data since it only has 2 time points.
#cor1 <- pt.a.div %>% select(-1) %>% select(-(3:6)) %>% select(-(27:53)) #only flu vaccine induced antibody
cor1 <- pt.a.div %>% select(-1) %>% select(-(3:6)) %>% select(-(51:53))
cor2 <- rcorr(as.matrix(cor1), type="pearson")
cor3 <- cbind(data.frame(cor2$r)[,(1:2)], data.frame(cor2$n)[,(1:2)], data.frame(cor2$P)[,(1:2)])
write.csv(cor3, file="test.csv")
# from the correlation matrix we found for simpson, h3n2d0(0.07), flu.iga.d0(0.08), flu.igg.fold(0.03), TLR2(0.07), TLR4(0.04); and for shannon, h3n2d0(0.054), h3n2.d14(0.08), flu.igg.fold(0.03), TLR2(0.07), TLR4(0.04)

## then look at the pattern of the above titers and simpson diversity
# h3n2.d0
ggplot(pt.a.div, aes(pt.a.div$h3n2.d0, pt.a.div$simpson))+geom_point()+geom_line()+ stat_smooth()
ggplot(pt.a.div, aes(factor(Sample_ID), h3n2.d0))+geom_point(aes(color=pt.a.div$simpson))#from the figure, it doesn't seems that it have an obvious trend between h3n2.d0 and simpson diversity
# flu.Iga.d0
ggplot(pt.a.div, aes(pt.a.div$flu.Iga.d0, pt.a.div$simpson))+geom_point()+geom_line()+ stat_smooth()
ggplot(pt.a.div, aes(factor(Sample_ID), flu.Iga.d0))+geom_point(aes(color=pt.a.div$simpson))#from the figure, it doesn't seems that it have an obvious trend between h3n2.d0 and simpson diversity
# flu.Igg.fold
ggplot(pt.a.div, aes(pt.a.div$flu.Igg.fold, pt.a.div$simpson))+geom_point()+geom_line()+ stat_smooth()
ggplot(pt.a.div, aes(factor(Sample_ID), flu.Igg.fold))+geom_point(aes(color=pt.a.div$simpson))




### 1.2 try to model the titer with simpson diversity and other factor
## 1.2.1 model h3n2.d0
# select the variable that is correlated to h3n2.d0
# cor4 <- data.frame(cor2$P)%>%select(h3n2.d0, flu.Iga.d0,flu.Igg.fold) %>% mutate(draw_rownames=row.names(cor2$P)) %>%filter(h3n2.d0<0.1|flu.Iga.d0<0.1|flu.Igg.fold<0.1)#this draw everthing
data.frame(cor2$P)%>% mutate(draw_rownames=row.names(cor2$P)) %>%filter(h3n2.d0<0.1)%>%select(draw_rownames)
summary(lm(h3n2.d0~simpson, data=pt.a.div))
summary(lm(h3n2.d0~simpson+age, data=pt.a.div)) # better model?
summary(lm(h3n2.d0~simpson+age+Race, data=pt.a.div))
summary(lm(h3n2.d0~simpson+age+gender, data=pt.a.div))  

#try to model h3n2.d14
summary(lm(h3n2.d14~simpson, data=pt.a.div))

## 1.2.1 model flu.Iga.d0
data.frame(cor2$P)%>% mutate(draw_rownames=row.names(cor2$P)) %>%filter(flu.Iga.d0<0.1)%>%select(draw_rownames)
summary(lm(flu.Iga.d0~simpson, data=pt.a.div))
summary(lm(flu.Iga.d0~simpson+age, data=pt.a.div))
summary(lm(flu.Iga.d0~simpson+gender, data=pt.a.div))# better model?

  