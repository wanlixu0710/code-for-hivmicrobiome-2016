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

###I have tried the combination of patients and healthy, but patients group making more sense with data analysis, and the stats is better
#demo
demo <- read.csv(file="sample information-microbiome analysis.csv") %>% filter(Sample_ID %in% pt.L6$Sample_ID) %>% select(1:5)
summary(demo)
sd(demo$age)


#### research question 1 ###	Anti-CD4 IgG-D0 (AC) and vaccine-induced anti-CD4 IgG (AD and AE)
### 1.1 data analysis on a-diversity vs titer.
## 1.1.1 first plot the titer level over 3 data points to see the pattern, and will use the ANCOVA for the model
# make a different format of the data
pt.a.div.long <-  melt(pt.a.div, id.vars=c("Sample_ID", "simpson","shannon","group","age","gender","Race")) %>% separate(col=variable, into =c("variable", "time.point") , sep="\\.") %>% filter(time.point=="d0"|time.point=="d7"|time.point=="d14")
pt.a.div.long$time.point <- revalue(pt.a.div.long$time.point,c("d0"=0, "d7"=7,"d14"=14))
library(tidyr)
pt.a.div.w <- spread(pt.a.div.long, variable, value)
pt.a.div.w$time.point <- as.numeric(pt.a.div.w$time.point)
#plot
ggplot(pt.a.div.w, aes(x=time.point, y=cd4_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()#+facet_grid(~group)

#check the normality of cd4_igg
hist(pt.a.div.w$cd4_igg, main=paste("Histogram of cd4_igg"))
hist(log(sqrt(pt.a.div.w$cd4_igg)))
hist((log(pt.a.div.w$cd4_igg))^(1/3))
## 1.1.1 ANCOVA
fit1 <- aov(((log(pt.a.div.w$cd4_igg))^(1/3)) ~ factor(time.point)+simpson, data=pt.a.div.w); summary(fit1)
fit2<- aov(((log(pt.a.div.w$cd4_igg))^(1/3)) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit2)
# get rid of the outliers and try it again

fit2<- aov(((log(pt.a.div.w[-c(22,24,27),]$cd4_igg))^(1/3)) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w[-c(22,24,27),]); summary(fit2)
fit2$coefficients

### 1.2 cd4_igg and beta diversity
## 1.2.1 nms with bray-curtis
pt.otu.all<- inner_join(otu, pt.a.div)
rownames(pt.otu.all) <- pt.otu.all$Sample_ID
nms <- metaMDS(pt.otu.all[,1:2490], dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d.bc <- inner_join(stat, pt.a.div)

# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd4_igg.d0
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d0, shape=group))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d0),shape=group)) #use log to reduce the difference

# cd4_igg.d14
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d14),shape=group))

# cd4_igg.dif
ggplot(stat.d.bc, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.fold),shape=group))

## 1.2.2 nms with jaccard
nms <- metaMDS(pt.otu.all[,1:2490], dist="jaccard", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d.ja <- inner_join(stat, pt.a.div)

# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# cd4_igg.d0
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d0))+scale_colour_gradient(low="green", high="red")
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d0),shape=group)) #use log to reduce the difference

# cd4_igg.d14
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.d14),shape=group))

# cd4_igg.dif
ggplot(stat.d.ja, aes(MDS1, MDS2))+geom_point(aes(color=factor(group), size=log(cd4_igg.fold),shape=group))

## 1.2.3 permanova bray-curtis
# cd4_igg.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$cd4_igg.d0 , na.rm=TRUE, permutations=99999, method="bray")#sig
# cd4_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd4_igg.d14 , na.rm=TRUE, permutations=99999, method="bray")#sig
# cd4_igg.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd4_igg.fold , na.rm=TRUE, permutations=99999, method="bray") # not sig


## 1.2.4 permanova jaccard
# cd4_igg.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$cd4_igg.d0 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# cd4_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd4_igg.d14 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# cd4_igg.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd4_igg.fold , na.rm=TRUE, permutations=99999, method="jaccard") # not sig

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
# cut cd4_igg.d14 by 50 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$cd4_igg.cut <- cut(pt.L6$cd4_igg.d14, breaks=c(-Inf, 50, Inf), labels=c("low", "high"))
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd4_igg.cut,control = how(nperm=350)))

## 1.3.3 cd4_igg_fold at genus level
# plot the cd4_igg to define the cutting point
summary(pt.a.div$cd4_igg.fold)
plot(pt.a.div$cd4_igg.fold)

# cut cd4_igg.fold by 2
pt.L6$cd4_igg.cut <- cut(pt.L6$cd4_igg.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"))
# indicator specis
summary( multipatt(pt.L6[,(5:245)], pt.L6$cd4_igg.cut,control = how(nperm=350)))

#### research question 2 ###	ANA IgG-D0 and vaccine-induced ANA IgG 
### 2.1 data analysis on a-diversity vs titer.
## 2.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(pt.a.div.w, aes(x=time.point, y=ana_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of ana_igg
hist(pt.a.div.w$ana_igg)
hist(pt.a.div.w$ana_igg^(1/3))
hist(log2 (pt.a.div.w$ana_igg))
hist(log2 ((pt.a.div.w$ana_igg)^(1/3)))

## 2.1.1 ANCOVA
fit1 <- aov((log2 ((pt.a.div.w$ana_igg)^(1/3))) ~ factor(time.point)+simpson, data=pt.a.div.w); summary(fit1)
fit2<- aov((log2 ((pt.a.div.w$ana_igg)^(1/3))) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)

fit1 <- glm(ana_igg ~ factor(time.point)+simpson,family=Gamma, data=pt.a.div.w); summary(fit1)
fit2 <- glm(ana_igg ~ factor(time.point)+simpson+gender+age,family=Gamma, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)


# get rid of the outliers and try it again

fit2<- aov(((log(pt.a.div.w[-c(22,24,27),]$cd4_igg))^(1/3)) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w[-c(22,24,27),]); summary(fit2)
fit2$coefficients

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
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$ana_igg.d0 , na.rm=TRUE, permutations=99999, method="bray")#sig
# ana_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igg.d14 , na.rm=TRUE, permutations=99999, method="bray")#sig
# ana_igg.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igg.fold , na.rm=TRUE, permutations=99999, method="bray") # not sig


## 2.2.4 permanova jaccard
# ana_igg.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$ana_igg.d0 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# ana_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igg.d14 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# ana_igg.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igg.fold , na.rm=TRUE, permutations=99999, method="jaccard") # not sig

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
# cut ana_igg.d14 by 0.9 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igg.cut <- cut(pt.L6$ana_igg.d14, breaks=c(-Inf, 0.9, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igg.cut,control = how(nperm=350)))

## 2.3.3 ana_igg_fold at genus level
# plot the ana_igg to define the cutting point
summary(pt.a.div$ana_igg.fold)
plot(pt.a.div$ana_igg.fold)
# cut ana_igg.fold by 2 
pt.L6 <- inner_join(L6.group,pt.a.div[-18,])
pt.L6$ana_igg.cut <- cut(pt.L6$ana_igg.fold, breaks=c(-Inf, 2, Inf), labels=c("low", "high"));summary(pt.L6$ana_igg.cut)
summary( multipatt(pt.L6[,(5:245)], pt.L6$ana_igg.cut,control = how(nperm=350)))


#### research question 3 ###	ANA IgM-D0 and vaccine-induced ANA IgG 
### 3.1 data analysis on a-diversity vs titer.
## 3.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using linear mixed effect model?
# make a different format of the data
#plot
ggplot(pt.a.div.w, aes(x=time.point, y=ana_igm, color=as.factor(Sample_ID)))+geom_point()+geom_line()

#check the normality of ana_igm
pt.a.div$ana_igm.d0
hist((1/(pt.a.div.w$ana_igm)))
#gamma distribution

## 3.1.1 ANCOVA
fit1 <- aov((1/(pt.a.div.w$ana_igm)) ~ factor(time.point)+simpson, data=pt.a.div.w); summary(fit1)
fit2<- aov((1/(pt.a.div.w$ana_igm)) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)
fit1$coefficients

fit1 <- glm((1/(pt.a.div.w$ana_igm))~ factor(time.point)+simpson,family=Gamma, data=pt.a.div.w); summary(fit1)
fit2 <- glm((1/(pt.a.div.w$ana_igm)) ~ factor(time.point)+simpson+gender+age,family=Gamma, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)


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
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$ana_igm.d0 , na.rm=TRUE, permutations=99999, method="bray")#sig
# ana_igm.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igm.d14 , na.rm=TRUE, permutations=99999, method="bray")
# ana_igm.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igm.fold , na.rm=TRUE, permutations=99999, method="bray")

## 2.2.4 permanova jaccard
# ana_igm.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$ana_igm.d0 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# ana_igm.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igm.d14 , na.rm=TRUE, permutations=99999, method="jaccard")
# ana_igm.fold
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$ana_igm.fold , na.rm=TRUE, permutations=99999, method="jaccard")

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
# cut ana_igg.d14 by 1.0 
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
ggplot(pt.a.div.w, aes(x=time.point, y=dsdna_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()

#check the normality of dsdna_igg
hist(pt.a.div.w$dsdna_igg)
hist(log(pt.a.div.w$dsdna_igg)) #better
hist(sqrt(1/pt.a.div.w$dsdna_igg))

## 4.1.1 ANCOVA
fit1 <- aov(log(pt.a.div.w$dsdna_igg) ~ factor(time.point)+simpson, data=pt.a.div.w); summary(fit1)
fit2<- aov(log(pt.a.div.w$dsdna_igg) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)
fit1$coefficients


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
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$dsdna_igg.d0 , na.rm=TRUE, permutations=99999, method="bray")#sig
# dsdna_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$dsdna_igg.d14 , na.rm=TRUE, permutations=99999, method="bray")

# dsdna_igg.dif
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$dsdna_igg.fold , na.rm=TRUE, permutations=99999, method="bray")

## 4.2.4 permanova jaccard
# dsdna_igg.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$dsdna_igg.d0 , na.rm=TRUE, permutations=99999, method="jaccard")#sig
# dsdna_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$dsdna_igg.d14 , na.rm=TRUE, permutations=99999, method="jaccard")

# dsdna_igg.dif
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$dsdna_igg.fold , na.rm=TRUE, permutations=99999, method="jaccard")

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
## 5.1.1 first plot the titer level over 3 data points to see the pattern, and possibly using ANCOVA
# make a different format of the data
#plot
ggplot(pt.a.div.w, aes(x=time.point, y=cd8_igg, color=as.factor(Sample_ID)))+geom_point()+geom_line()+facet_grid(~group)

#check the normality of cd8_igg
hist(pt.a.div.w$cd8_igg)
hist(log10 (pt.a.div.w$cd8_igg)) #better

## 5.1.1 ANCOVA
fit1 <- aov((log10 (pt.a.div.w$cd8_igg)) ~ factor(time.point)+simpson, data=pt.a.div.w); summary(fit1)
fit2<- aov((log10 (pt.a.div.w$cd8_igg)) ~ factor(time.point)+simpson+gender+age, data=pt.a.div.w); summary(fit2)
BIC(fit1)
BIC(fit2)
plot(fit1)
fit1$coefficients


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
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$cd8_igg.d0, na.rm=TRUE, permutations=99999, method="bray")#sig
# cd8_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd8_igg.d14 , na.rm=TRUE, permutations=99999, method="bray")

# cd8_igg.dif
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd8_igg.fold , na.rm=TRUE, permutations=99999, method="bray")

## 5.2.4 permanova Jaccard
# cd8_igg.d0
adonis(pt.otu.all[,(1:2490)] ~pt.otu.all$age+ pt.otu.all$gender+ pt.otu.all$cd8_igg.d0, na.rm=TRUE, permutations=99999, method="jaccard")#sig
# cd8_igg.d14
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd8_igg.d14 , na.rm=TRUE, permutations=99999, method="jaccard")

# cd8_igg.dif
adonis(pt.otu.all[pt.otu.all$Sample_ID!="Vac34",(1:2490)] ~pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$age+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$gender+ pt.otu.all[pt.otu.all$Sample_ID!="Vac34",]$cd8_igg.fold , na.rm=TRUE, permutations=99999, method="jaccard")

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

