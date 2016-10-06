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
library(devEMF)
library(tidyr)
rm(list=ls())

#### first, lets read in data clinical and bacteria
otu <- read.table(file="2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/120215JW515F-pr.fasta.otus.fa.OTU.percentages.txt", sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)
otu2 <- otu[,-1]%>%select(-37)
colSums (otu2, na.rm = FALSE, dims = 1) # rarified to 100
otu <- as.data.frame(t(otu2)) 
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

clinical <- read.csv(file="updated sample information_8_31_2016.csv", sep=",", header=T)

a.div <- inner_join(a.div, clinical)

L6 <- read.csv(file="2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/120215JW515F-pr.fasta.otus.fa.genus.percentages.txt", sep="\t", header=T, stringsAsFactors = F)
L6.2 <- L6[,-38] %>% select(-1) 
L6.3 <- data.frame(t(L6.2))
colnames(L6.3) <- L6$X
neg <- L6.3[(35:36),] ### stop here!!!
negmean <- colMeans(neg)
negmean <- data.frame(negmean)
negmean <- t(negmean)
L6.2 <- L6.3-negmean
L6 <- L6.2[-(35:36),]
L6[L6 < 0] <- 0
L6$Sample_ID <- rownames(L6)
a <- a.div[,c("Sample_ID","group", "gender", "Race")]
L6.group <- inner_join(a, L6)

### only focus on pt group
pt.a.div <- filter(a.div, group=="Patient")
pt.otu <- filter(otu, Sample_ID %in% pt.a.div$Sample_ID)
pt.L6 <- filter(L6.group, group=="Patient")

#### difference between patients and healthy control############
# a-diversity
hist(a.div$simpson)
summary(lm(simpson~group+gender+age, data=a.div))
t.test(a.div$simpson~ a.div$group)# not significant
# beta-diversity
otu.per <- inner_join(otu,a.div)
br_dist <- vegdist(otu.per[,(1:2490)], method="bray")
str(br_dist)
adonis(br_dist ~otu.per$group, na.rm=TRUE, permutations=999)#sig
betadisper(br_dist, otu.per$group) # not significant
# indicator species
summary( multipatt(L6.group[,(5:245)], L6.group$group,control = how(nperm=350)))

#### research question 1 ###	Anti-CD4 IgG-D0 (AC) and vaccine-induced anti-CD4 IgG (AD and AE)
### 1.1 data analysis on a-diversity vs titer.
## 1.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
a.div.long <-  melt(a.div, id.vars=c("Sample_ID", "simpson","shannon","group","age","gender","Race")) %>% separate(col=variable, into =c("variable", "time.point") , sep="\\.") %>% filter(time.point=="d0"|time.point=="d7"|time.point=="d14")
a.div.long$time.point <- revalue(a.div.long$time.point,c("d0"=0, "d7"=7,"d14"=14))
library(tidyr)
a.div.w <- spread(a.div.long, variable, value)
a.div.w$time.point <- as.numeric(a.div.w$time.point)
#plot
ggplot(a.div.w, aes(x=time.point, y=cd4_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of cd4_igg
hist(a.div.w$cd4_igg, main=paste("Histogram of cd4_igg"))
hist(log(sqrt(a.div.w$cd4_igg)))
hist((log(a.div.w$cd4_igg))^(1/3))

library(lme4)
library(lmerTest)


null.mod <- lmer(((log(a.div.w$cd4_igg))^(1/3))~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod)

null.mod1 <-  lmer(((log(a.div.w$cd4_igg))^(1/3))~time.point+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod1)
anova(null.mod,null.mod1)

null.mod2 <- lmer(((log(a.div.w$cd4_igg))^(1/3))~time.point+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(null.mod2)
anova(null.mod1, null.mod2)

mod1 <- lmer(((log(a.div.w$cd4_igg))^(1/3))~time.point+group+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(mod1)
anova(null.mod2, mod1)

mod2 <- lmer(((log(a.div.w$cd4_igg))^(1/3))~time.point+group+simpson+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(mod2)

anova(mod2, mod1) # not siginificant

## 1.1.2 check the fold and a-diversity with linear regression model # not significant
hist((log(pt.a.div$cd4_igg.fold))^(1/3))
summary(lm(((log(a.div$cd4_igg.d0))^(1/3))~simpson, data=a.div))
summary(lm(((log(a.div$cd4_igg.d0))^(1/3))~simpson+gender+age, data=a.div))

## 1.1.3 latent class model
library(gridExtra)
library(grid)
library(survival)
library(lcmm)
str(a.div.w$time.point)
a.div.w$cd4_igg.trs <- (log(a.div.w$cd4_igg))^(1/3)
m1cd4 <- hlme(cd4_igg.trs~time.point, random=~time.point, subject='Sample_ID',ng=1,data=a.div.w, cor="AR"(time.point), maxiter=5e3)
summary(m1cd4)

m2cd4<-hlme(cd4_igg.trs~time.point, random=~time.point, mixture=~time.point, classmb=~simpson+group+gender, cor="AR"(time.point),subject='Sample_ID',ng=2,data=a.div.w,B=m1cd4)
summary(m2cd4) # not significant

people1 <- as.data.frame(m2cd4$pprob[,1:2]) 
fd1 <- left_join(a.div.w, people1)
p1 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F) + labs(x="time point",y="simpson",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Sample_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T) + labs(x="time point",y="simpson",colour="Latent Class");p2 

### 1.2 cd4_igg and beta diversity
## 1.2.1 nms with bray-curtis
otu.all<- inner_join(otu, a.div)
nms <- metaMDS(otu[,1:2490], dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d.bc <- inner_join(stat, a.div)

# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd4_igg.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d0, shape=group))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d0),shape=group)) #use log to reduce the difference

# cd4_igg.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d14),shape=group))

# cd4_igg.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.fold),shape=group))

## 1.2.2 nms with jaccard
nms <- metaMDS(otu[,1:2490], dist="jaccard", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d.ja <- inner_join(stat, a.div)

# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd4_igg.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d0))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d0),shape=group)) #use log to reduce the difference

# cd4_igg.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d14),shape=group))

# cd4_igg.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.fold),shape=group))

## 1.2.3 permanova bray-curtis
otu.per <- inner_join(otu,pt.a.div[-18,])
br_dist <- vegdist(otu.per[,(1:2490)], method="bray")
str(br_dist)
# cd4_igg.d0
adonis(br_dist ~otu.per$cd4_igg.d0, na.rm=TRUE, permutations=999)#sig
# cd4_igg.d14
adonis(br_dist ~otu.per$cd4_igg.d14, na.rm=TRUE, permutations=999)#sig
# cd4_igg.fold
adonis(br_dist ~otu.per$cd4_igg.fold, na.rm=TRUE, permutations=999)# not sig

## 1.2.4 permanova jaccard
otu.per <- inner_join(otu,pt.a.div[-18,])
ja_dist <- vegdist(otu.per[,(1:2490)], method="jaccard")
# cd4_igg.d0
adonis(ja_dist ~otu.per$cd4_igg.d0, na.rm=TRUE, permutations=999)#sig
# cd4_igg.d14
adonis(ja_dist ~otu.per$cd4_igg.d14, na.rm=TRUE, permutations=999)#sig
# cd4_igg.fold
adonis(ja_dist ~otu.per$cd4_igg.fold, na.rm=TRUE, permutations=999)# not sig

### 1.3 indicator species
## 1.3.1 cd4_igg_d0 at genus level
# plot the cd4_igg to define the cutting point
pt.a.div.w <- a.div.w %>% filter(group=="Patient")
plot(pt.a.div.w$cd4_igg, col=factor(pt.a.div.w$time.point))
summary(pt.a.div$cd4_igg.d0)
plot(pt.a.div$cd4_igg.d0)
# cut cd4_igg.d0 by 50 
pt.L6 <- inner_join(L6.group,pt.a.div)
pt.L6$cd4_igg.cut <- cut(pt.L6$cd4_igg.d0, breaks=c(-Inf, 50, Inf), labels=c("low", "high"))
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd4_igg.cut, control = how(nperm=350)))

## 1.3.2 cd4_igg_d14 at genus level
# plot the cd4_igg to define the cutting point
summary(pt.a.div$cd4_igg.d14)
plot(pt.a.div$cd4_igg.d14)
# cut cd4_igg.d0 by 50 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$cd4_igg.cut <- cut(pt.L6$cd4_igg.d14, breaks=c(-Inf, 50, Inf), labels=c("low", "high"))
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd4_igg.cut,control = how(nperm=350)))

## 1.3.3 cd4_igg_fold at genus level
# plot the cd4_igg to define the cutting point
summary(pt.a.div$cd4_igg.fold)
plot(pt.a.div$cd4_igg.fold)

# cut cd4_igg.d0 by 50 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$cd4_igg.cut <- cut(pt.L6$cd4_igg.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"))

# indicator specis
summary( multipatt(L6.group[,(5:245)], L6.group$group,control = how(nperm=350)))

#### research question 2 ###	ANA IgG-D0 and vaccine-induced ANA IgG 
### 2.1 data analysis on a-diversity vs titer.
## 2.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(a.div.w, aes(x=time.point, y=ana_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of ana_igg
hist(a.div.w$ana_igg)
hist(log((a.div.w$ana_igg)^(1/3)))
a.div.w$ana_igg_trs <- log((a.div.w$ana_igg)^(1/3))

null.mod <- lmer(ana_igg_trs~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w[a.div.w$group=="Patient",], REML=F); summary(null.mod)
null.mod <- lmer(ana_igg_trs~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod)

null.mod1 <-  lmer(ana_igg_trs~time.point+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod1)
anova(null.mod,null.mod1)

null.mod2 <- lmer(ana_igg_trs~time.point+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(null.mod2)
anova(null.mod1, null.mod2)

mod1 <- lmer(ana_igg_trs~time.point+simpson+age+group+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(mod1)
anova(null.mod2, mod1)
# not siginificant with linear effect model.

## 2.1.2 check the fold and a-diversity with linear regression model # not significant
summary(lm((log((ana_igg.d0)^(1/3)))~simpson, data=pt.a.div))
summary(lm((log((ana_igg.fold)^(1/3)))~simpson, data=pt.a.div))
summary(lm((log((ana_igg.fold)^(1/3)))~simpson+gender+age, data=pt.a.div))

## 2.1.3 latent class model
m1cd4 <- hlme(ana_igg_trs~time.point, random=~time.point, subject='Sample_ID',ng=1,data=a.div.w, cor="AR"(time.point), maxiter=5e3)
summary(m1cd4)

m2cd4<-hlme(ana_igg_trs~time.point, random=~time.point, mixture=~time.point, classmb=~group+simpson+gender, cor="AR"(time.point),subject='Sample_ID',ng=2,data=a.div.w,B=m1cd4)
summary(m2cd4) ## didn't converge

people1 <- as.data.frame(m2cd4$pprob[,1:2]) 
fd1 <- left_join(a.div.w, people1)
p1 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F) + labs(x="time point",y="simpson",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Sample_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T) + labs(x="time point",y="simpson",colour="Latent Class");p2 

### 2.2 ana_igg and beta diversity
## 2.2.1 nms with bray-curtis
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# ana_igg.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=ana_igg.d0))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.d0)))+scale_colour_gradient(low="green", high="red") #use log to reduce the difference

# ana_igg.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.d14)))+scale_colour_gradient(low="green", high="red")

# ana_igg.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 2.2.2 nms with jaccard
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# ana_igg.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=ana_igg.d0))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.d0)))+scale_colour_gradient(low="green", high="red") #use log to reduce the difference

# ana_igg.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.d14)))+scale_colour_gradient(low="green", high="red")

# ana_igg.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 2.2.3 permanova bray-curtis
# ana_igg.d0
adonis(br_dist ~otu.per$ana_igg.d0, na.rm=TRUE, permutations=999)# not sig
# ana_igg.d14
adonis(br_dist ~otu.per$ana_igg.d14, na.rm=TRUE, permutations=999)#not sig
# ana_igg.fold
adonis(br_dist ~otu.per$ana_igg.fold, na.rm=TRUE, permutations=999)# not sig

## 2.2.4 permanova jaccard
# ana_igg.d0
adonis(ja_dist ~otu.per$ana_igg.d0, na.rm=TRUE, permutations=999) #not sig
# ana_igg.d14
adonis(ja_dist ~otu.per$ana_igg.d14, na.rm=TRUE, permutations=999) #not sig
# ana_igg.fold
adonis(ja_dist ~otu.per$ana_igg.fold, na.rm=TRUE, permutations=999) # not sig

### 2.3 indicator species
## 2.3.1 ana_igg_d0 at genus level
# plot the ana_igg to define the cutting point
plot(pt.a.div.w$ana_igg, col=factor(pt.a.div.w$time.point))
summary(pt.a.div$ana_igg.d0)
plot(pt.a.div$ana_igg.d0)
# cut ana_igg.d0 by 0.4 
pt.L6 <- inner_join(L6.group,pt.a.div)
pt.L6$ana_igg.cut <- cut(pt.L6$ana_igg.d0, breaks=c(-Inf, 0.4, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igg.cut,control = how(nperm=350))) ## gut bacteria in low group

## 2.3.2 ana_igg_d14 at genus level
# plot the ana_igg to define the cutting point
summary(pt.a.div$ana_igg.d14)
plot(pt.a.div$ana_igg.d14)
# cut ana_igg.d0 by 0.9 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igg.cut <- cut(pt.L6$ana_igg.d14, breaks=c(-Inf, 0.9, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igg.cut,control = how(nperm=350)))

## 2.3.3 ana_igg_fold at genus level
# plot the ana_igg to define the cutting point
summary(pt.a.div$ana_igg.fold)
plot(pt.a.div$ana_igg.fold)
# cut ana_igg.d0 by 2 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igg.cut <- cut(pt.L6$ana_igg.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igg.cut,control = how(nperm=350)))


#### research question 3 ###	ANA IgM-D0 and vaccine-induced ANA IgG 
### 3.1 data analysis on a-diversity vs titer.
## 3.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(a.div.w, aes(x=time.point, y=ana_igm, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of ana_igm
qqnorm(a.div.w$ana_igm)
hist(a.div.w$ana_igm)
hist(10^(a.div.w$ana_igm))


null.mod <- lmer(ana_igm~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w[a.div.w$group=="Patient",], REML=F); summary(null.mod)
null.mod <- lmer(ana_igm~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod)

null.mod1 <-  lmer(ana_igm~time.point+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod1)
anova(null.mod,null.mod1)

null.mod2 <- lmer(ana_igm~time.point+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(null.mod2)
anova(null.mod1, null.mod2)

mod1 <- lmer(ana_igm~time.point+simpson+gender+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(mod1)
anova(null.mod1, mod1)
# not siginificant with linear effect model.

## 3.1.2 check the fold and a-diversity with linear regression model # not significant
summary(lm(ana_igm.d0~simpson, data=pt.a.div))
summary(lm(ana_igm.fold~simpson, data=pt.a.div))
summary(lm(ana_igm.fold~simpson+gender+age, data=pt.a.div))

## 3.1.3 latent class model ****
m1cd4 <- hlme(ana_igm~time.point, random=~time.point, subject='Sample_ID',ng=1,data=a.div.w, cor="AR"(time.point), maxiter=5e3)
summary(m1cd4)

m2cd4<-hlme(ana_igm~time.point, random=~time.point, mixture=~time.point, classmb=~simpson+gender, cor="AR"(time.point),subject='Sample_ID',ng=2,data=a.div.w,B=m1cd4)
summary(m2cd4) 

people1 <- as.data.frame(m2cd4$pprob[,1:2])
fd1 <- left_join(a.div.w, people1)
p1 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F) + labs(x="time point",y="simpson",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(time.point, ana_igm, group=Sample_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Sample_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T) + labs(x="time point",y="simpson",colour="Latent Class");p2 

### 3.2 ana_igm and beta diversity
## 3.2.1 nms with bray-curtis
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# ana_igm.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d0))+scale_colour_gradient(low="green", high="red")

# ana_igm.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d14))+scale_colour_gradient(low="green", high="red")

# ana_igm.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igm.fold)))+scale_colour_gradient(low="green", high="red")

## 3.2.2 nms with jaccard
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# ana_igm.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d0))+scale_colour_gradient(low="green", high="red")

# ana_igm.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d14))+scale_colour_gradient(low="green", high="red")

# ana_igm.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(ana_igm.fold)))+scale_colour_gradient(low="green", high="red")

## 3.2.3 permanova bray-curtis
# ana_igm.d0
adonis(br_dist ~otu.per$ana_igm.d0, na.rm=TRUE, permutations=999)# not sig
# ana_igm.d14
adonis(br_dist ~otu.per$ana_igm.d14, na.rm=TRUE, permutations=999)#not sig
# ana_igm.fold
adonis(br_dist ~otu.per$ana_igm.fold, na.rm=TRUE, permutations=999)# not sig

## 2.2.4 permanova jaccard
# ana_igm.d0
adonis(ja_dist ~otu.per$ana_igm.d0, na.rm=TRUE, permutations=999) #not sig
# ana_igm.d14
adonis(ja_dist ~otu.per$ana_igm.d14, na.rm=TRUE, permutations=999) #not sig
# ana_igm.fold
adonis(ja_dist ~otu.per$ana_igm.fold, na.rm=TRUE, permutations=999) # not sig

### 3.3 indicator species
## 3.3.1 ana_igm_d0 at genus level
# plot the ana_igm to define the cutting point
plot(pt.a.div.w$ana_igm, col=factor(pt.a.div.w$time.point))
summary(pt.a.div$ana_igm.d0)
plot(pt.a.div$ana_igg.d0)
# cut ana_igg.d0 by 0.5 
pt.L6 <- inner_join(L6.group,pt.a.div)
pt.L6$ana_igm.cut <- cut(pt.L6$ana_igm.d0, breaks=c(-Inf, 0.5, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igm.cut,control = how(nperm=350))) ## gut bacteria in low group

## 3.3.2 ana_igm_d14 at genus level
# plot the ana_igm to define the cutting point
summary(pt.a.div$ana_igm.d14)
plot(pt.a.div$ana_igm.d14)
# cut ana_igg.d0 by 1.0 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igm.cut <- cut(pt.L6$ana_igm.d14, breaks=c(-Inf, 1.0, Inf), labels=c("low", "high"));summary(pt.L6$ana_igm.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igm.cut,control = how(nperm=350)))

## 3.3.3 ana_igm_fold at genus level
# plot the ana_igm to define the cutting point
summary(pt.a.div$ana_igm.fold)
plot(pt.a.div$ana_igm.fold)
# cut ana_igm_fold by 2.0
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igm.cut <- cut(pt.L6$ana_igm.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"));summary(pt.L6$ana_igm.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igm.cut,control = how(nperm=350)))


#### research question 4 ###	Anti-dsDNA IgG-D0 and vaccin-induced anti-dsDNA IgG 
### 4.1 data analysis on a-diversity vs titer.
## 4.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(a.div.w, aes(x=time.point, y=dsdna_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of dsdna_igg
hist(a.div.w$dsdna_igg)
hist(log(a.div.w$dsdna_igg)) #better
hist(sqrt(1/a.div.w$dsdna_igg))

null.mod <- lmer(log(dsdna_igg)~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w[a.div.w$group=="Patient",], REML=F); summary(null.mod)
null.mod <- lmer(log(dsdna_igg)~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod)

null.mod1 <-  lmer(log(dsdna_igg)~time.point+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod1)
anova(null.mod,null.mod1)

null.mod2 <- lmer(log(dsdna_igg)~time.point+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(null.mod2)
anova(null.mod1, null.mod2)

mod1 <- lmer(log(dsdna_igg)~time.point+simpson+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(mod1)
anova(null.mod2, mod1)
# not siginificant with linear effect model.

## 4.1.2 check the fold and a-diversity with linear regression model # not significant
summary(lm(log(dsdna_igg.d0)~simpson, data=pt.a.div))
summary(lm(log(dsdna_igg.fold)~simpson, data=pt.a.div))
summary(lm(log(dsdna_igg.fold)~simpson+gender+age, data=pt.a.div))

## 4.1.3 latent class model
a.div.w$log.dsdna_igg <- log(a.div.w$dsdna_igg)
m1cd4 <- hlme(log.dsdna_igg~time.point, random=~time.point, subject='Sample_ID',ng=1,data=a.div.w, cor="AR"(time.point), maxiter=5e3)
summary(m1cd4)

m2cd4<-hlme(log.dsdna_igg~time.point, random=~time.point, mixture=~time.point, classmb=~simpson, cor="AR"(time.point),subject='Sample_ID',ng=2,data=a.div.w,B=m1cd4)
summary(m2cd4) 

people1 <- as.data.frame(m2cd4$pprob[,1:2])
fd1 <- left_join(a.div.w, people1)
p1 <- ggplot(fd1, aes(time.point, log.dsdna_igg, group=Sample_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F) + labs(x="time point",y="simpson",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(time.point, log.dsdna_igg, group=Sample_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Sample_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T) + labs(x="time point",y="simpson",colour="Latent Class");p2 

### 4.2 dsdna_igg and beta diversity
## 4.2.1 nms with bray-curtis
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# dsdna_igg.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.d0)))+scale_colour_gradient(low="green", high="red")

# dsdna_igg.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.d14)))+scale_colour_gradient(low="green", high="red")

# dsdna_igg.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 4.2.2 nms with jaccard
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# dsdna_igg.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.d0)))+scale_colour_gradient(low="green", high="red")

# dsdna_igg.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.d14)))+scale_colour_gradient(low="green", high="red")

# dsdna_igg.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(dsdna_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 4.2.3 permanova bray-curtis
# dsdna_igg.d0
adonis(br_dist ~otu.per$dsdna_igg.d0, na.rm=TRUE, permutations=999) # not sig
# dsdna_igg.d14
adonis(br_dist ~otu.per$dsdna_igg.d14, na.rm=TRUE, permutations=999)#not sig
# dsdna_igg.fold
adonis(br_dist ~otu.per$dsdna_igg.fold, na.rm=TRUE, permutations=999)# not sig

## 4.2.4 permanova jaccard
# dsdna_igg.d0
adonis(ja_dist ~otu.per$dsdna_igg.d0, na.rm=TRUE, permutations=999) #not sig
# dsdna_igg.d14
adonis(ja_dist ~otu.per$dsdna_igg.d14, na.rm=TRUE, permutations=999) #not sig
# dsdna_igg.fold
adonis(ja_dist ~otu.per$dsdna_igg.fold, na.rm=TRUE, permutations=999) # not sig

### 4.3 indicator species
## 4.3.1 dsdna_igg.d0 at genus level
# plot the dsdna_igg to define the cutting point
plot(pt.a.div.w$dsdna_igg, col=factor(pt.a.div.w$time.point))
summary(pt.a.div$dsdna_igg.d0)
plot(pt.a.div$dsdna_igg.d0)
# cut dsdna_igg.d0 by 25 
pt.L6 <- inner_join(L6.group,pt.a.div)
pt.L6$dsdna_igg.cut <- cut(pt.L6$dsdna_igg.d0, breaks=c(-Inf, 25, Inf), labels=c("low", "high"));summary(pt.L6$dsdna_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$dsdna_igg.cut,control = how(nperm=350))) ## gut bacteria in low group

## 4.3.2 dsdna_igg_d14 at genus level
# plot the dsdna_igg to define the cutting point
summary(pt.a.div$dsdna_igg.d14)
plot(pt.a.div$dsdna_igg.d14)
# cut dsdna_igg.d0 by 50 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$dsdna_igg.cut <- cut(pt.L6$dsdna_igg.d14, breaks=c(-Inf, 50, Inf), labels=c("low", "high"));summary(pt.L6$dsdna_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$dsdna_igg.cut,control = how(nperm=350)))

## 4.3.3 dsdna_igg.fold at genus level
# plot the dsdna_igg to define the cutting point
summary(pt.a.div$dsdna_igg.fold)
plot(pt.a.div$dsdna_igg.fold)
# cut dsdna_igg_fold by 2.0
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$dsdna_igg.cut <- cut(pt.L6$dsdna_igg.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"));summary(pt.L6$dsdna_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$dsdna_igg.cut,control = how(nperm=350)))


#### research question 5 ###	Anti-CD8-D0 and vaccine-induced anti-CD8 IgG 
### 5.1 data analysis on a-diversity vs titer.
## 5.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(a.div.w, aes(x=time.point, y=cd8_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of cd8_igg
hist(a.div.w$cd8_igg)
hist(log(a.div.w$cd8_igg)) #better
hist(sqrt(1/a.div.w$cd8_igg))

null.mod <- lmer(log(cd8_igg)~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w[a.div.w$group=="Patient",], REML=F); summary(null.mod)
null.mod <- lmer(log(cd8_igg)~time.point+(1|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod)

null.mod1 <-  lmer(log(cd8_igg)~time.point+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(null.mod1)
anova(null.mod,null.mod1)

null.mod2 <- lmer(log(cd8_igg)~time.point+(1+time.point|Sample_ID), data=a.div.w, REML=F); summary(null.mod2)
anova(null.mod1, null.mod2)

mod1 <- lmer(log(cd8_igg)~time.point+simpson+(1+time.point|Sample_ID)+(1|time.point), data=a.div.w, REML=F); summary(mod1)
anova(null.mod2, mod1)
# not siginificant with linear effect model.

## 5.1.2 check the fold and a-diversity with linear regression model # not significant
summary(lm(cd8_igg.d0~simpson, data=pt.a.div))
summary(lm(log(cd8_igg.fold)~simpson, data=pt.a.div))
summary(lm(log(cd8_igg.fold)~simpson+gender+age, data=pt.a.div))

## 5.1.3 latent class model
a.div.w$log.cd8_igg <- log(a.div.w$cd8_igg)
m1cd4 <- hlme(log.cd8_igg~time.point, random=~time.point, subject='Sample_ID',ng=1,data=a.div.w, cor="AR"(time.point), maxiter=5e3)
summary(m1cd4)

m2cd4<-hlme(log.cd8_igg~time.point, random=~time.point, mixture=~time.point, classmb=~gender, cor="AR"(time.point),subject='Sample_ID',ng=2,data=a.div.w,B=m1cd4)
summary(m2cd4) 

people1 <- as.data.frame(m2cd4$pprob[,1:2])
fd1 <- left_join(a.div.w, people1)
p1 <- ggplot(fd1, aes(time.point, log.cd8_igg, group=Sample_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F) + labs(x="time point",y="log.cd8_igg",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(time.point, log.cd8_igg, group=Sample_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Sample_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T) + labs(x="time point",y="simpson",colour="Latent Class");p2 

### 5.2 cd8_igg and beta diversity
## 5.2.1 nms with bray-curtis
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd8_igg.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.d0)))+scale_colour_gradient(low="green", high="red")

# cd8_igg.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.d14)))+scale_colour_gradient(low="green", high="red")

# cd8_igg.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 5.2.2 nms with jaccard
# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd8_igg.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.d0)))+scale_colour_gradient(low="green", high="red")

# cd8_igg.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.d14)))+scale_colour_gradient(low="green", high="red")

# cd8_igg.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=log(cd8_igg.fold)))+scale_colour_gradient(low="green", high="red")

## 5.2.3 permanova bray-curtis
# cd8_igg.d0
adonis(br_dist ~otu.per$cd8_igg.d0, na.rm=TRUE, permutations=999) # not sig
# cd8_igg.d14
adonis(br_dist ~otu.per$cd8_igg.d14, na.rm=TRUE, permutations=99999)#almost sig
# cd8_igg.fold
adonis(br_dist ~otu.per$cd8_igg.fold, na.rm=TRUE, permutations=9999)# not sig

## 5.2.4 permanova jaccard
# cd8_igg.d0
adonis(ja_dist ~otu.per$cd8_igg.d0, na.rm=TRUE, permutations=999) #not sig
# cd8_igg.d14
adonis(ja_dist ~otu.per$cd8_igg.d14, na.rm=TRUE, permutations=99999) # sig
# cd8_igg.fold
adonis(ja_dist ~otu.per$cd8_igg.fold, na.rm=TRUE, permutations=999) # not sig

### 5.3 indicator species
## 5.3.1 cd8_igg.d0 at genus level
# plot the cd8_igg to define the cutting point
plot(pt.a.div.w$cd8_igg, col=factor(pt.a.div.w$time.point))
summary(pt.a.div$cd8_igg.d0)
plot(pt.a.div$cd8_igg.d0)
# cut cd8_igg.d0 by 0.5 
pt.L6 <- inner_join(L6.group,pt.a.div)
pt.L6$cd8_igg.cut <- cut(pt.L6$cd8_igg.d0, breaks=c(-Inf, 0.5, Inf), labels=c("low", "high"));summary(pt.L6$cd8_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd8_igg.cut,control = how(nperm=350))) ## gut bacteria in low group

## 5.3.2 cd8_igg_d14 at genus level
# plot the cd8_igg to define the cutting point
summary(pt.a.div$cd8_igg.d14)
plot(pt.a.div$cd8_igg.d14)
# cut cd8_igg.d0 by 1.0
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$cd8_igg.cut <- cut(pt.L6$cd8_igg.d14, breaks=c(-Inf, 1.0, Inf), labels=c("low", "high"));summary(pt.L6$cd8_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd8_igg.cut,control = how(nperm=350)))

## 5.3.3 cd8_igg.fold at genus level
# plot the cd8_igg to define the cutting point
summary(pt.a.div$cd8_igg.fold)
plot(pt.a.div$cd8_igg.fold)
# cut cd8_igg_fold by 2.0
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$cd8_igg.cut <- cut(pt.L6$cd8_igg.fold, breaks=c(-Inf, 1.49, Inf), labels=c("low", "high"));summary(pt.L6$cd8_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd8_igg.cut,control = how(nperm=350)))

