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
otu$Sample_ID <- rownames(otu)
a.div <- otu[,c("Sample_ID","simpson")]

clinical <- read.csv(file="sample information-microbiome analysis.csv", sep=",", header=T)
clinical$bcell.dif <- clinical$bcell.d7-clinical$bcell.d0
clinical$cd3.dif <- clinical$cd3.d7-clinical$cd3.d0
clinical$cd4.dif <- clinical$cd4.d7-clinical$cd4.d0

a.div <- inner_join(a.div, clinical)

### 1. compare patient and control ###
## demographic
summary(a.div[,(2:5)])
sd(a.div$age)
summary(factor(a.div$female))
17/(15+17)
summary(factor(a.div$group))
21/(11+21)
pa <- a.div %>% filter (group=="patient")
summary(pa[,(2:5)])
sd(pa$age)
summary(factor(pa$female))
14/21
12/21
he <- a.div %>% filter (group=="healthy")
summary(he[,(2:5)])
sd(he$age)
summary(factor(he$female))
3/11
4/11

## 1.1 first, let's see whether there are demographic difference between the two groups
# age
fig1 <- plot(a.div$group, a.div$age)
plot(a.div$group, a.div$age)
t.test(a.div$age~a.div$group) # there is no age difference between groups
# race
table(factor(a.div$group), factor(a.div$Race))
chisq.test(table(factor(a.div$group), factor(a.div$Race))) #p=0.53
# gender
table(factor(a.div$group), factor(a.div$gender))
chisq.test(table(factor(a.div$group), factor(a.div$gender))) # significant, healthy has more female while patient has more male, need to be controled in later model
t.test(a.div$age~a.div$gender) #almost significant; female are 7 years older than male
summary(lm(a.div$age~a.div$gender+a.div$group)) #not significant

## 1.2 the difference in simpson between pt and healthy
plot(a.div$group, a.div$simpson) # similar between two groups
# then plot the simpson with gender/group, to see the difference
ggplot(a.div,  aes( y=as.numeric(simpson), x=factor(group)))+ geom_boxplot(aes(colour=group))+ scale_x_discrete()+labs(x="pt vs healthy", y="simpson", title="simpson between pt vs control")+ facet_wrap(~gender, scales="free_y")+theme_bw()
# liner regression model
lm.simpson <- lm(a.div$simpson~a.div$group*a.div$gender)
summary(lm.simpson)# it looks like there is trend, maybe with increased samples size, it will be signigicant, also checked race and age, no significant difference

## 1.3 nms comparision, didn't show specific patterns
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

## 1.4 permanova
## Calculation of bray curtis dissimilarity
br_dist <- vegdist(otu.nms, method="bray")
adonis(br_dist ~a.div$group, na.rm=TRUE, permutations=99) #not significant


## 1.5 genus level comparision 
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

## indicator species
# indicator species
L6.is <- L6.group[,-(1:4)]
FT  <-  multipatt(L6.is, L6.group$group,
                  control = how(nperm=350))
summary(FT)

## 1.6 now lets check the differences in antibody titer between the pt. and healthy, the data use clinical (full data with samples dont have microbiome data)
clinical.fig <- melt(clinical, id.vars=c("Sample_ID","group", "gender", "Race"))
fig1 <- ggplot(clinical.fig,  aes( y=as.numeric(value), x=factor(group)))+ geom_boxplot(aes(colour=group))+ scale_x_discrete()+ labs(x="male vs. female", y="antibody", title="antibody boxplot between patients and healthy")+ facet_wrap(~variable, scales="free_y")+ theme_bw();fig1
# choose the one looks different on the figure and run linear regression
summary(lm(clinical$h1n1.d7~clinical$group)) #pt is higher
summary(lm(clinical$flu.Igg.d7~clinical$group)) #pt is higher
summary(lm(clinical$flu.Igm.d7~clinical$group)) #pt is higher
summary(lm(clinical$flu.Igm.fold~clinical$group))#pt is higher
summary(lm(clinical$cd4.igg.d0~clinical$group)) #pt is higher
summary(lm(clinical$cd4.igg.d14~clinical$group)) #pt is higher
summary(lm(clinical$cd4.igg.dif~clinical$group)) #pt is higher
summary(lm(clinical$ana.igg.d14~clinical$group)) #pt is higher
summary(lm(clinical$ana.igg.dif~clinical$group)) #pt is higher
summary(lm(clinical$ana.igm.d0~clinical$group)) #pt is almost significant higher
summary(lm(clinical$ana.igm.d14~clinical$group)) #pt is higher
summary(lm(clinical$ana.igg.dif~clinical$group)) #pt is higher
summary(lm(clinical$cd8.igg.d14~clinical$group)) #pt is higher
summary(lm(clinical$cd8.igg.dif~clinical$group)) #pt is almost significant higher
summary(lm(clinical$TLR2~clinical$group)) #pt is higher
summary(lm(clinical$TLR4~clinical$group)) #pt is higher
summary(lm(clinical$LPS.plasma~clinical$group)) #pt is higher

# fig2 <- ggplot(clinical.fig,  aes( y=as.numeric(value), x=factor(gender)))+ geom_boxplot(aes(colour=gender))+ scale_x_discrete()+ labs(x="male vs. female", y="antibody", title="antibody boxplot between male and female")+ facet_wrap(~variable, scales="free_y")+ theme_bw();fig2


### 2. Now we are moving to the patients only
pt.a.div <- filter(a.div, group=="Patient")
pt.otu <- filter(otu, Sample_ID %in% pt.a.div$Sample_ID)
pt.L6 <- filter(L6.group, group=="Patient")

## 2.1 simpson score
# check the simpson and differnt titers by plotting
pt.a.div.plot <- melt(pt.a.div, id.vars=c("Sample_ID","group","simpson","Race","age","gender"))
fig2 <- ggplot(pt.a.div.plot,  aes( y=simpson, x=value)) + geom_point()+  geom_smooth(method="lm", se=FALSE)+  labs(x="value", y="simpson", title="titers and predictors vs. simpson")+ facet_wrap(~variable, scales="free")+ theme(axis.text.x = element_text(angle = 45)) 
fig2

## 2.1 cut antibody titer to two levels for the graphing purpose
table(pt.a.div$cd4.igg.dif)
pt.a.div$cd4.igg <- cut(pt.a.div$cd4.igg.dif,c(-Inf,37,Inf), c("low","high"))

table(pt.a.div$dsdna.igg.dif)
pt.a.div$dsdna.igg <- cut(pt.a.div$dsdna.igg.dif,c(-Inf,18,Inf), c("low","high"))

table(pt.a.div$ana.igg.dif)
pt.a.div$ana.igg <- cut(pt.a.div$ana.igg.dif,c(-Inf,0.5,Inf), c("low","high"))

table(pt.a.div$ana.igm.dif)
pt.a.div$ana.igm <- cut(pt.a.div$ana.igm.dif,c(-Inf,0.5,Inf), c("low","high"))

table(pt.a.div$cd8.igg.dif)
pt.a.div$cd8.igg <- cut(pt.a.div$cd8.igg.dif,c(-Inf,0.25,Inf), c("low","high"))

table(pt.a.div$bcell.dif)
pt.a.div$bcell <- cut(pt.a.div$bcell.dif,c(-Inf,50,Inf), c("low","high"))

table(pt.a.div$cd3.dif)
pt.a.div$bcell <- cut(pt.a.div$bcell.dif,c(-Inf,0,Inf), c("negative","positive"))

table(pt.a.div$cd4.dif)
pt.a.div$bcell <- cut(pt.a.div$bcell.dif,c(-Inf,0,Inf), c("negative","positive"))

# for neutralization and ELISA specific antibody, 4 folders is the cut-off 
pt.a.div$h1n1 <- cut(pt.a.div$h1n1.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$h1n1, pt.a.div$h1n1.fold)

pt.a.div$b <- cut(pt.a.div$b.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$b, pt.a.div$b.fold)

pt.a.div$h3n2 <- cut(pt.a.div$h3n2.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$h3n2, pt.a.div$h3n2.fold)

pt.a.div$flu.iga <- cut(pt.a.div$flu.Iga.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$flu.iga, pt.a.div$flu.Iga.fold)

pt.a.div$flu.igg <- cut(pt.a.div$flu.Igg.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$flu.igg, pt.a.div$flu.Igg.fold)

pt.a.div$flu.igm <- cut(pt.a.div$flu.Igm.fold,c(-Inf,4,Inf), c("low","high"), right=FALSE)
table(pt.a.div$flu.igm, pt.a.div$flu.Igm.fold)

#TLR baseline
table(pt.a.div$TLR2)
pt.a.div$TLR2.base <- cut(pt.a.div$TLR2,c(-Inf,3,Inf), c("low","high"))

table(pt.a.div$TLR4)
pt.a.div$TLR4.base <- cut(pt.a.div$TLR4,c(-Inf,4,Inf), c("low","high"))

table(pt.a.div$TLR7)
pt.a.div$TLR7.base <- cut(pt.a.div$TLR7,c(-Inf,0.9,Inf), c("low","high"))

table(pt.a.div$TLR9)
pt.a.div$TLR9.base <- cut(pt.a.div$TLR9,c(-Inf,1,Inf), c("low","high"))

# plot the simpson boxplot
pt.a.div.group <- pt.a.div[, c(1:2, 5:6, 56, 61:75)]
pt.a.div.plot <- melt(pt.a.div.group, id.vars=c("Sample_ID", "simpson","LPS.plasma"))
ggplot(pt.a.div.plot,  aes(y=simpson, x=value)) + geom_boxplot() +  labs(x="value", y="simpson", title="titers and predictors vs. simpson")+ facet_wrap(~variable, scales="free")+ theme_bw()

# now try with some of the possibilities that has been indicated, throw the data in the spss, and found simpson is correlated to flu.igg.fold and TL4 (p<0.05), there are three are close to significant: h3n2.d0(p=0.072), flu.iga.d0 (p=0.078) and TLR2(p=0.068), then I tryed the partial correlation control for variables and try to build the model, now I am echoing the model using linear regresssion in R
summary(lm(simpson~ flu.Igg.fold, data=pt.a.div))
summary(lm(simpson~ flu.Igg.fold+ h3n2.d0, data=pt.a.div))
AIC(lm(simpson~ flu.Igg.fold+ h3n2.d0, data=pt.a.div))
summary(lm(simpson~ flu.Igg.fold+ h3n2.d0 + b.fold, data=pt.a.div)) #model 1
AIC(lm(simpson~ flu.Igg.fold+ h3n2.d0 + b.fold, data=pt.a.div))

summary(lm(simpson~ flu.Igg.fold + flu.Iga.d0, data=pt.a.div))
summary(lm(simpson~ flu.Igg.fold + flu.Iga.d0+ h3n2.d0, data=pt.a.div))
summary(lm(simpson~ flu.Igg.fold + flu.Iga.d0+ h3n2.d0 +b.fold, data=pt.a.div))
AIC(lm(simpson~ flu.Igg.fold + flu.Iga.d0+ h3n2.d0 +b.fold, data=pt.a.div)) # model 2, not as good as the first model

## 2.2 LPS-Plasma
ggplot(pt.a.div.plot,  aes(y=LPS.plasma, x=value)) + geom_boxplot() +  labs(x="value", y="LPS.plasma", title="titers and predictors vs. LPS.plasma")+ facet_wrap(~variable, scales="free")+ theme_bw()

# now try with some of the possibilities that has been indicated, throw the data in the spss, and there are several variables are significant, therefore I used forward stepwise selection from spss (the forward stepwise selection doesn't work in R), since TLRs have a lot of missing data and the correlation is not significant, they are excluded in the forward stepwise selection). Now echo the selection in linear regression. also use the bivariate correlation (pearson) in the spss , there are several variables that are sig: gender, bd0, cd4.igg.d0, cd4.igg.d7, cd4.igg.d14, cd4.igg.dif, dsdna.igg.d14, ana.igg.d0 (p=0.051), ana.igg.d7, ana.igg.d14, ana.igg.dif, ana.igm.d14, cd8.igg.dif)
pt.a.div$LPS.plasma <- as.numeric(as.character(pt.a.div$LPS.plasma))
summary(lm(LPS.plasma ~ cd8.igg.dif, data=pt.a.div))
summary(lm(LPS.plasma ~ cd8.igg.dif + gender, data=pt.a.div))
summary(lm(LPS.plasma ~ cd8.igg.dif + gender+ cd4.igg.dif, data=pt.a.div))
summary(lm(LPS.plasma ~gender+ cd4.igg.dif, data=pt.a.div))
summary(lm(LPS.plasma ~gender+ cd4.igg.dif + h3n2.d0, data=pt.a.div))
AIC(lm(LPS.plasma ~gender+ cd4.igg.dif + h3n2.d0, data=pt.a.div))

## 2.2 genus level vs titer
# firstly put the data in the spss and try to get the correlation from spss, the result is saved as excel in the folder


# 2.3 indicator species @ genus level
pt.L6 <- pt.L6[,-(2:4)]
pt.L6.full <- inner_join(pt.a.div[,c(1, 4:5, 61:76)], pt.L6)
# indicator species
pt.L6.clinical <- pt.L6.full[,(1:19)]
pt.L6.indi <- pt.L6.full[,(20:260)]
rownames(pt.L6.indi) <- pt.L6.clinical$Sample_ID

# indicator species for cd4.igg difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, cd4.igg) %>% filter(!is.na(cd4.igg))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$cd4.igg, control = how(nperm=350))
summary(indi.sp)

# indicator species for dsdna.igg difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, dsdna.igg) %>% filter(!is.na(dsdna.igg))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$dsdna.igg, control = how(nperm=350))
summary(indi.sp)

# indicator species for ana.igg difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, ana.igg) %>% filter(!is.na(ana.igg))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$ana.igg, control = how(nperm=350))
summary(indi.sp)

# indicator species for ana.igm difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, ana.igm) %>% filter(!is.na(ana.igm))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$ana.igm, control = how(nperm=350))
summary(indi.sp)

# indicator species for cd8.igg difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, cd8.igg) %>% filter(!is.na(cd8.igg))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$cd8.igg, control = how(nperm=350))
summary(indi.sp)

# indicator species for bcell difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, bcell) %>% filter(!is.na(bcell))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$bcell, control = how(nperm=350))
summary(indi.sp)

# indicator species for h1n1 difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, h1n1) %>% filter(!is.na(h1n1))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$h1n1, control = how(nperm=350))
summary(indi.sp)

# indicator species for b difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, b) %>% filter(!is.na(b))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$b, control = how(nperm=350))
summary(indi.sp)

# indicator species for h3n2 difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, h3n2) %>% filter(!is.na(h3n2))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$h3n2, control = how(nperm=350))
summary(indi.sp)

# indicator species for flu.iga difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, flu.iga) %>% filter(!is.na(flu.iga))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$flu.iga, control = how(nperm=350))
summary(indi.sp)

# indicator species for flu.igg difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, flu.igg) %>% filter(!is.na(flu.igg))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$flu.igg, control = how(nperm=350))
summary(indi.sp)

# indicator species for flu.igm difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, flu.igm) %>% filter(!is.na(flu.igm))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$flu.igm, control = how(nperm=350))
summary(indi.sp)

# indicator species for TLR2.base difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, TLR2.base) %>% filter(!is.na(TLR2.base))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$TLR2.base, control = how(nperm=350))
summary(indi.sp)

# indicator species for TLR4.base difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, TLR4.base) %>% filter(!is.na(TLR4.base))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$TLR4.base, control = how(nperm=350))
summary(indi.sp)

# indicator species for TLR7.base difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, TLR7.base) %>% filter(!is.na(TLR7.base))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$TLR7.base, control = how(nperm=350))
summary(indi.sp)

# indicator species for TLR9.base difference 
indi.data <- pt.L6.clinical %>% select(Sample_ID, TLR9.base) %>% filter(!is.na(TLR9.base))
indi.L6 <- pt.L6.indi[rownames(pt.L6.indi) %in% indi.data$Sample_ID,]
indi.sp <- multipatt(indi.L6, indi.data$TLR9.base, control = how(nperm=350))
summary(indi.sp)

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

## then look at the pattern of the vaccine titers and simpson diversity
pt.a.div.plot <- melt(pt.a.div, id.vars=c("Sample_ID","group","simpson","Race","age","gender","shannon"))
fig1 <- ggplot(pt.a.div.plot,  aes( y=simpson, x=value)) + geom_point()+  geom_smooth(method="lm", se=FALSE)+  labs(x="value", y="simpson", title="titers and predictors vs. simpson")+ facet_wrap(~variable, scales="free")+ theme(axis.text.x = element_text(angle = 45));fig1

fig2 <- ggplot(pt.a.div.plot,  aes( y=simpson, x=value)) + geom_point()+  geom_smooth()+  labs(x="value", y="simpson", title="titers and predictors vs. simpson")+ facet_wrap(~variable, scales="free")+ theme(axis.text.x = element_text(angle = 45));fig2

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


####notes
####we could group the patients and healthy together
cor4 <- a.div %>% select(-1) %>% select(-(3:6)) %>% select(-(51:53))
cor6 <- rcorr(as.matrix(cor4), type="pearson")
cor7 <- cbind(data.frame(cor6$r)[,(1:2)], data.frame(cor6$n)[,(1:2)], data.frame(cor6$P)[,(1:2)])
write.csv(cor7, file="test.csv")

