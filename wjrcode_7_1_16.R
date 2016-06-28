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
rm(list=ls())

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
otu$Sample_ID <- rownames(otu)
a.div <- otu[,c("Sample_ID","simpson")]

clinical <- read.csv(file="sample information-microbiome analysis.csv", sep=",", header=T)
clinical <- clinical[,-(29:49)]
clinical$bcell.dif <- clinical$bcell.d7-clinical$bcell.d0
clinical$cd3.dif <- clinical$cd3.d7-clinical$cd3.d0
clinical$cd4.dif <- clinical$cd4.d7-clinical$cd4.d0

a.div <- inner_join(a.div, clinical)


### 1. compare patient and control ###
### 1.1 first, let's see whether there are demographic difference between the two groups
# age vs group
plot(a.div$group, a.div$age, main="age and group")
t.test(a.div$age~a.div$group) # there is no age difference between groups
# race vs. group
table(factor(a.div$group), factor(a.div$Race))
chisq.test(table(factor(a.div$group), factor(a.div$Race))) #p=0.53
# gender vs group
table(factor(a.div$group), factor(a.div$gender))
chisq.test(table(factor(a.div$group), factor(a.div$gender))) # significant, healthy has more female while patient has more male, need to be controled in later model
t.test(a.div$age~a.div$gender) #almost significant; female are 7 years older than male
summary(lm(a.div$age~a.div$gender+a.div$group)) #not significant

## 1.2 the difference in simpson between pt and healthy
plot(a.div$group, a.div$simpson, main="simpson and group") # similar between two groups
# then plot the simpson with gender/group, to see the difference
ggplot(a.div,  aes( y=as.numeric(simpson), x=factor(group)))+ geom_boxplot(aes(colour=group))+ scale_x_discrete()+labs(x="pt vs healthy", y="simpson", title="simpson between pt vs control")+ facet_wrap(~gender, scales="free_y")+theme_bw()
# liner regression model
lm.simpson <- lm(a.div$simpson~a.div$group*a.div$gender)
summary(lm.simpson)# it looks like there is trend, maybe with increased samples size, it will be signigicant, also checked race and age, no significant difference

### 1.3 nms comparision, didn't show specific patterns
otu.nms <-otu[,-(2491:2492)]
nms <- metaMDS(otu.nms, dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d <- inner_join(stat, a.div)
plot(x=stat.d$MDS1,y=stat.d$MDS2, col=stat.d$group) # didn't show specific pattern
plot(x=stat.d$MDS1,y=stat.d$MDS2, col=stat.d$gender) 

nms <- metaMDS(otu.nms, dist="bray", k=3, trymax=250, wascores=TRUE, trymin=50)
stressplot(nms)
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d <- inner_join(stat, a.div)
scatterplot3d(stat.d[,1:3], color=as.integer(stat.d$group), pch=16, grid=TRUE, box=FALSE) #didn't show specific patterns
scatterplot3d(stat.d[,1:3], color=as.integer(stat.d$gender), pch=16, grid=TRUE, box=FALSE) 

### 1.4 permanova
## Calculation of bray curtis dissimilarity
br_dist <- vegdist(otu.nms, method="bray")
adonis(br_dist ~a.div$group, na.rm=TRUE, permutations=99) #not significant


### 1.5 genus level comparision 
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
genus <- melt(L6.group, id.vars=c("Sample_ID","group", "gender", "Race"))
# healthy vs. patients genus level plot
genus <- genus[order(genus$variable, decreasing = TRUE),]
ge1 <- ggplot(genus,  aes( y=as.numeric(value), x=factor(group)))+ geom_boxplot(aes(colour=group))+ scale_x_discrete()+ labs(x="healthy vs. patients", y="percentage", title="genus boxplot between Healthy and patient after controlling for neg")+ facet_wrap(~variable, scales="free_y")+ theme_bw();ge1
# several bacteria are different between groups from the figure, however, the liner regression doesn't show statistical difference.
summary(lm(L6.group$atopostipes~L6.group$group)) # high in healthy
summary(lm(L6.group$bifidobacterium~L6.group$group)) # high in healthy

# male vs. female level plot
ge2 <- ggplot(genus,  aes( y=as.numeric(value), x=factor(gender)))+ geom_boxplot(aes(colour=gender))+ scale_x_discrete()+ labs(x="male vs. female", y="percentage", title="genus boxplot between male and female")+ facet_wrap(~variable, scales="free_y")+ theme_bw();ge2

### 1.6 indicator species 
# indicator species vs. group
L6.is <- L6.group[,-(1:4)]
indisp  <-  multipatt(L6.is, L6.group$group,
                  control = how(nperm=350))
summary(indisp)

# indicator species vs. gender
L6.is <- L6.group[,-(1:4)]
indisp <-  multipatt(L6.is, L6.group$gender,
                  control = how(nperm=350))
summary(indisp)



### 1.3 the difference in the antibodys between different group
clinical.fig <- melt(clinical, id.vars=c("Sample_ID","group", "gender", "Race"))
fig1 <- ggplot(clinical.fig,  aes( y=as.numeric(value), x=factor(group)))+ geom_boxplot(aes(colour=group))+ scale_x_discrete()+ labs(x="male vs. female", y="antibody", title="antibody boxplot between patients and healthy")+ facet_wrap(~variable, scales="free_y")+ theme_bw();fig1
