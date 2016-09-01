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
library(tidyr)
rm(list=ls())

### research question
## could try to see the patients group or combine patients & control together (if group itself doesnt have sig)
# 1.which bateria is related ANA. DsDNA, antiCD4 (D0, D14 or fold) #genus level
# 2. baseline  antiCD4 is related to bacteria?
# 3. which batereia is related to CD4 (most important), bcell (less important), focus on baseline

#### first, lets read in data clinical and bacteria
otu <- read.table(file="2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/120215JW515F-pr.fasta.otus.fa.OTU.percentages.txt", sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)
otu2 <- otu[,-1]%>%select(-37)
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

######## 1.which bateria is related ANA. DsDNA, antiCD4 (D0, D14 or fold), andtiCD8  genus level
### 1.1 a diversity (only focus on patients)
## simple correlation
cor1 <- pt.a.div %>% select(simpson, shannon, ana_igg.d0, ana_igg.d7, ana_igg.d14, ana_igg.fold, ana_igm.d0, ana_igm.d7, ana_igm.d14, ana_igm.fold, dsdna_igg.d0, dsdna_igg.d7, dsdna_igg.d14, dsdna_igg.fold, cd4_igg.d0, cd4_igg.d7,cd4_igg.d14, cd4_igg.fold, cd8_igg.d0,cd8_igg.d7, cd8_igg.d14, cd8_igg.fold)
library(Hmisc)
cor2 <- rcorr(as.matrix(cor1), type="pearson")
cor3 <- cbind(data.frame(cor2$r)[,(1:2)], data.frame(cor2$n)[,(1:2)], data.frame(cor2$P)[,(1:2)])
cor4 <- data.frame(cor2$P)[,(1:2)]%>% mutate(draw.rowname=row.names(cor2$P))%>%filter(simpson<0.1|shannon<0.1) # results shows non of the variables are correlated with simpson or shannon diversity

## try to model the titer as the outcome variable, as the outcome variables have 3 data points, so maybe could try to use mixed effect model.
# make a different format of the data
pt.a.div.long <- melt(pt.a.div, id.vars=c("Sample_ID", "simpson","shannon","group","age","gender","Race")) %>% separate(col=variable, into =c("variable", "time.point") , sep="\\.") %>% filter(time.point=="d0"|time.point=="d7"|time.point=="d14")
pt.a.div.long$time.point <- revalue(pt.a.div.long$time.point,c("d0"=0, "d7"=7,"d14"=14))
library(tidyr)
pt.a.div.w <- spread(pt.a.div.long, variable, value)
pt.a.div.w$time.point <- as.numeric(pt.a.div.w$time.point)

# ana_igg, time (time is not a fixed effect, as it is vaccine induced not time, and vaccine is baseline for everyone) and sample are random effect, and model simpson, gender and age as race as fixed effect
library(lme4)
library(lmerTest) #get the p value from anova test
fm1 <- lmer(ana_igg~1+(time.point|Sample_ID), pt.a.div.w, REML=F);fm1
anova(fm1)
fm2 <- lmer(ana_igg~simpson+gender+age+Race+(time.point|Sample_ID), pt.a.div.w, REML=F);fm2
anova(fm2) # the model AIC is higher than the basic model

# ana_igm
fm1 <- lmer(ana_igm~1+(time.point|Sample_ID), pt.a.div.w, REML=F);fm1
anova(fm1)
fm2 <- lmer(ana_igm~simpson+gender+age+Race+(time.point|Sample_ID), pt.a.div.w, REML=F);fm2
anova(fm2) # not significant

# cd4_igg
fm1 <- lmer(cd4_igg~1+(time.point|Sample_ID), pt.a.div.w, REML=F);fm1
anova(fm1)
fm2 <- lmer(cd4_igg~simpson+gender+Race+(time.point|Sample_ID), pt.a.div.w, REML=F);fm2
anova(fm2) # not significant

# dsDNA_igg
fm1 <- lmer(dsdna_igg~1+(time.point|Sample_ID), pt.a.div.w, REML=F);fm1
anova(fm1)
fm2 <- lmer(dsdna_igg~simpson+gender+Race+(time.point|Sample_ID), pt.a.div.w, REML=F);fm2
anova(fm2) # not significant

# cd8_igg
fm1 <- lmer(cd8_igg~1+(time.point|Sample_ID), pt.a.div.w, REML=F);fm1
anova(fm1)
fm2 <- lmer(cd8_igg~simpson+age+gender+Race+(time.point|Sample_ID), pt.a.div.w, REML=F);fm2
anova(fm2) # not significant

### 1.2 simpson diversity including patients and healthy control
## simpson correlation
cor1 <- a.div %>% select(simpson, shannon, ana_igg.d0, ana_igg.d7, ana_igg.d14, ana_igg.fold, ana_igm.d0, ana_igm.d7, ana_igm.d14, ana_igm.fold, dsdna_igg.d0, dsdna_igg.d7, dsdna_igg.d14, dsdna_igg.fold, cd4_igg.d0, cd4_igg.d7,cd4_igg.d14, cd4_igg.fold, cd8_igg.d0,cd8_igg.d7, cd8_igg.d14, cd8_igg.fold)
cor2 <- rcorr(as.matrix(cor1), type="pearson")
cor3 <- cbind(data.frame(cor2$r)[,(1:2)], data.frame(cor2$n)[,(1:2)], data.frame(cor2$P)[,(1:2)])
cor4 <- data.frame(cor2$P)[,(1:2)]%>% mutate(draw.rowname=row.names(cor2$P))%>%filter(simpson<0.1|shannon<0.1) # not significant

# NO trend has been found with simpson and auto-antibodies. Move to beta-diversity

### 1.2 b-diversity pt only
otu.pt <- otu[otu$Sample_ID %in% pt.a.div$Sample_ID,]
otu.nms <-otu.pt[,-(2491:2493)]
nms <- metaMDS(otu.nms, dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50) # first try to fit in 2 dimentions
stat <- data.frame(nms$points)
stat$Sample_ID <- rownames(stat)
stat.d <- inner_join(stat, pt.a.div)

# Plot the nms plot with different vaccine induced titers (AUTO ANTIBODIES)
# try to wrap all the plots in the same picture (different scales therefore color were not presented well, need to use individual plot instead)
#stat.m <- stat.d %>% select(MDS1, MDS2, age, gender, Race,ana_igg.d0, ana_igg.d7, ana_igg.d14, ana_igg.fold, ana_igm.d0, ana_igm.d7, ana_igm.d14, ana_igm.fold, dsdna_igg.d0, dsdna_igg.d7, dsdna_igg.d14, dsdna_igg.fold, cd4_igg.d0, cd4_igg.d7,cd4_igg.d14, cd4_igg.fold, cd8_igg.d0,cd8_igg.d7, cd8_igg.d14, cd8_igg.fold) %>% melt(id.vars=c("MDS1","MDS2","age","gender","Race"))
#ggplot(stat.m, aes(stat.m$MDS1, stat.m$MDS2))+geom_point(aes(color=stat.m$value))+scale_colour_gradient(low="green", high="red")+facet_wrap(~variable)

# cd4_igg.d0
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d0))+scale_colour_gradient(low="green", high="red")
# cd4_igg.d14
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.d14))+scale_colour_gradient(low="green", high="red")
# cd4_igg.fold
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd4_igg.fold))+scale_colour_gradient(low="green", high="red")

# ana_igg.d0
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igg.d0))+scale_colour_gradient(low="green", high="red")
# ana_igg.d14
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igg.d14))+scale_colour_gradient(low="green", high="red")
# ana_igg.fold
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igg.fold))+scale_colour_gradient(low="green", high="red")

# ana_igm.d0
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d0))+scale_colour_gradient(low="green", high="red")
# ana_igm.d14
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.d14))+scale_colour_gradient(low="green", high="red")
# ana_igm.fold
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=ana_igm.fold))+scale_colour_gradient(low="green", high="red")

# dsdna_igg.d0
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=dsdna_igg.d0))+scale_colour_gradient(low="green", high="red")
# dsdna_igg.d14
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=dsdna_igg.d14))+scale_colour_gradient(low="green", high="red")
# dsdna_igg.fold
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=dsdna_igg.fold))+scale_colour_gradient(low="green", high="red")

# cd8_igg.d0
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd8_igg.d0))+scale_colour_gradient(low="green", high="red")
# cd8_igg.d14
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd8_igg.d14))+scale_colour_gradient(low="green", high="red")
# cd8_igg.fold
ggplot(stat.d, aes(MDS1, MDS2))+geom_point(aes(color=cd8_igg.fold))+scale_colour_gradient(low="green", high="red")

#????? needs to look at the pattens of each pictures and have the envfit work!

