##Calculating phylogenetic diversity in each ecoregion before and after extinction of VU and above

#Loading essential package
library(ape)
library(apTreeshape)
library(geiger)
library(caper)
library(picante)

#Reading nexus tree file of reef species
cat("Reading trees...","\n")
supertree<-read.nexus(file="Supertree.tre")

#Reading species distribution file (Veron et al. 2009)
distribution<-read.csv(file="Distribution.csv",header=TRUE,row.names=1)

#Removing non-reef species from data
reef<-array(0,dim=length(distribution[,1]))
for (i in 1:length(distribution[,1])){
	reef[i]<-sum(distribution[i,1:141])
	}
distribution<-cbind(distribution,reef)
distribution<-subset(distribution,reef!=0)
distribution<-distribution[,1:141]

#Removing non-reef species from trees
cat("Pruning trees to reef...","\n")
overlap<-geiger:::.treedata(supertree[[1]],distribution)
reeftree<-rtreeshape(n=length(supertree),tip.number=length(distribution[,1]),model="yule")
for (i in 1:1000){
	reeftree[[i]]<-drop.tip(supertree[[i]],overlap$tree_not_data)
	}
print(geiger:::.treedata(reeftree[[1]],distribution))

#Calculating phylogenetic diversity in each ecoregion before extinction of VU and above
pd_original<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
colnames(pd_original)<-colnames(distribution)
for (i in 1:1000){
	cat("Observed PD : Tree",i,"\n")
	pd_original[i,]<-pd(t(distribution),reeftree[[i]])[,1]
	}
write.csv(pd_original,file="Ecoregion_PD_original.csv")

#Reading species data file (Carpenter et al. 2008)
species<-read.csv(file="Species.csv",header=TRUE,row.names=1)
species[,1][is.na(species[,1])]<-"DD"
species[,1]<-factor(species[,1],levels=c("LC","DD","NT","VU","EN","CR"))
print(geiger:::.treedata(reeftree[[1]],species))

#Calculating phylogenetic diversity in each ecoregion after extinction of VU and above
pd_observed<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
pd_null<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
colnames(pd_observed)<-colnames(distribution)
colnames(pd_null)<-colnames(distribution)
for (i in 1:length(distribution[1,])){
	cat("Ecoregion",i,"\n")
	cat("Pruning trees to region...","\n")
	overlap<-geiger:::.treedata(reeftree[[1]],distribution[distribution[,i]==1,])
	regiontree<-rtreeshape(n=length(reeftree),tip.number=sum(distribution[,i]),model="yule")
	for (j in 1:1000){
		regiontree[[j]]<-drop.tip(reeftree[[j]],overlap$tree_not_data)
		}
	print(geiger:::.treedata(regiontree[[1]],distribution[distribution[,i]==1,]))
	sparedtips<-rownames(distribution[(distribution[,i]==1&as.numeric(species[,1])<4),])
	for (j in 1:1000){
		cat("Ecoregion",i,": Observed extinction : Tree",j,"\n")
		pd_observed[j,i]<-pd.calc(cm=regiontree[[j]],tip.subset=sparedtips,method="TBL")[1]
		}
	for (j in 1:1000){
		cat("Ecoregion",i,": Random extinction : Tree",j,"\n")
		pd_null[j,i]<-mean(pd.bootstrap(cm=regiontree[[j]],ntips=length(sparedtips),reps=1000,method="TBL"))
		}
	write.csv(pd_observed,file="Ecoregion_PD_observed.csv")
	write.csv(pd_null,file="Ecoregion_PD_null.csv")
	}

#Initialising array to store excess PD loss data
excessPD<-array(0,dim=c(1000,length(distribution[1,])))
colnames(excessPD)<-colnames(distribution)
for (i in 1:length(distribution[1,])){
	excessPD[,i]<-(pd_null[,i]-pd_observed[,i])/pd_null[,i]
	}
write.csv(excessPD,file="Ecoregion_PD_excessPD.csv")

#Initialising array to store results
output<-array(0,dim=c(length(distribution[1,]),14))
colnames(output)<-c(
	"Species richness","Extinction rate",
	"Mean original PD","SD original PD",
	"Mean observed PD","SD observed PD",
	"Mean null PD","SD null PD",
	"Mean Excess PD","SD Excess PD","Lower CI","Upper CI",
	"T-statistic","P"
	)
rownames(output)<-colnames(distribution)

#Performing T-test to compare observed PD with null PD
for (i in 1:length(distribution[1,])){
	output[i,1]<-sum(distribution[,i])
	output[i,2]<-sum(distribution[(distribution[,i]==1&as.numeric(species[,1])>=4),i])/sum(distribution[,i])
	output[i,3]<-mean(pd_original[,i])
	output[i,4]<-sd(pd_original[,i])
	output[i,5]<-mean(pd_observed[,i])
	output[i,6]<-sd(pd_observed[,i])
	output[i,7]<-mean(pd_null[,i])
	output[i,8]<-sd(pd_null[,i])
	output[i,9]<-mean(excessPD[,i])
	output[i,10]<-sd(excessPD[,i])
	output[i,11]<-t.test(excessPD[,i],mu=0,conf.level=0.95)$conf.int[1]
	output[i,12]<-t.test(excessPD[,i],mu=0,conf.level=0.95)$conf.int[2]
	output[i,13]<-t.test(excessPD[,i],mu=0,conf.level=0.95)$statistic
	output[i,14]<-t.test(excessPD[,i],mu=0,conf.level=0.95)$p.value
	}

#Writing output to file
write.csv(output,file="Ecoregion_PD_output.csv")

##################################################

##Calculating phylogenetic species variability in each ecoregion before and after extinction of VU and above

#Loading essential package
library(ape)
library(apTreeshape)
library(geiger)
library(picante)

#Reading nexus tree file of reef species
cat("Reading trees...","\n")
supertree<-read.nexus(file="Supertree.tre")

#Reading species distribution file (Veron et al. 2009)
distribution<-read.csv(file="Distribution.csv",header=TRUE,row.names=1)

#Removing non-reef species from data
reef<-array(0,dim=length(distribution[,1]))
for (i in 1:length(distribution[,1])){
	reef[i]<-sum(distribution[i,1:141])
	}
distribution<-cbind(distribution,reef)
distribution<-subset(distribution,reef!=0)
distribution<-distribution[,1:141]

#Removing non-reef species from trees
cat("Pruning trees to reef...","\n")
overlap<-geiger:::.treedata(supertree[[1]],distribution)
reeftree<-rtreeshape(n=length(supertree),tip.number=length(distribution[,1]),model="yule")
for (i in 1:length(supertree)){
	reeftree[[i]]<-drop.tip(supertree[[i]],overlap$tree_not_data)
	}
print(geiger:::.treedata(reeftree[[1]],distribution))

#Calculating phylogenetic species variability in each ecoregion before extinction of VU and above
psv_original<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
colnames(psv_original)<-colnames(distribution)
for (i in 1:1000){
	cat("Observed PSV : Tree",i,"\n")
	psv_original[i,]<-psv(t(distribution),reeftree[[i]])[,1]
	}
write.csv(psv_original,file="Ecoregion_PSV_original.csv")

#Reading species data file (Carpenter et al. 2008)
species<-read.csv(file="Species.csv",header=TRUE,row.names=1)
species[,1][is.na(species[,1])]<-"DD"
species[,1]<-factor(species[,1],levels=c("LC","DD","NT","VU","EN","CR"))
print(geiger:::.treedata(reeftree[[1]],species))

#Generating species distribution after extinction of VU and above
distribution_spared<-distribution
distribution_spared[as.numeric(species[,1])>=4,]<-0

#Calculating phylogenetic species variability in each ecoregion after extinction of VU and above
psv_observed<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
psv_null<-array(NA,dim=c(length(reeftree),length(distribution[1,])))
colnames(psv_observed)<-colnames(distribution)
colnames(psv_null)<-colnames(distribution)
for (i in 1:length(reeftree)){
	cat("Tree",i,": Observed extinction","\n")
	psv_observed[i,]<-psv(t(distribution_spared),reeftree[[i]])[,1]
	nreps=1000
	psv_random<-array(NA,dim=c(nreps,length(distribution[1,])))
	for (j in 1:nreps){
		cat("Tree",i,": Random extinction : Rep",j,"\n")
		distribution_random<-distribution
		for (k in 1:length(distribution[1,])){
			randomtips<-sample(which(distribution[,k]==1),sum(distribution[,k])-sum(distribution_spared[,k]))
			if (length(randomtips)>0) {distribution_random[randomtips,k]<-0}
			}
		psv_random[j,]<-psv(t(distribution_random),reeftree[[i]])[,1]
		}
	psv_null[i,]<-colMeans(psv_random)
	write.csv(psv_observed,file="Ecoregion_PSV_observed.csv")
	write.csv(psv_null,file="Ecoregion_PSV_null.csv")
	}

#Initialising array to store excess PSV loss data
excessPSV<-array(0,dim=c(length(reeftree),length(distribution[1,])))
colnames(excessPSV)<-colnames(distribution)
for (i in 1:length(distribution[1,])){
	excessPSV[,i]<-(psv_null[,i]-psv_observed[,i])/psv_null[,i]
	}
write.csv(excessPSV,file="Ecoregion_PSV_excessPSV.csv")

#Initialising array to store results
output<-array(0,dim=c(length(distribution[1,]),14))
colnames(output)<-c(
	"Species richness","Extinction rate",
	"Mean original PSV","SD original PSV",	
	"Mean observed PSV","SD observed PSV",
	"Mean null PSV","SD null PSV",
	"Mean Excess PSV","SD Excess PSV","Lower CI","Upper CI",
	"T-statistic","P"
	)
rownames(output)<-colnames(distribution)

#Performing T-test to compare observed PSV with null PSV
for (i in 1:length(distribution[1,])){
	output[i,1]<-sum(distribution[,i])
	output[i,2]<-sum(distribution[(distribution[,i]==1&as.numeric(species[,1])>=4),i])/sum(distribution[,i])
	output[i,3]<-mean(psv_original[,i])
	output[i,4]<-sd(psv_original[,i])
	output[i,5]<-mean(psv_observed[,i])
	output[i,6]<-sd(psv_observed[,i])
	output[i,7]<-mean(psv_null[,i])
	output[i,8]<-sd(psv_null[,i])
	output[i,9]<-mean(excessPSV[,i])
	output[i,10]<-sd(excessPSV[,i])
	output[i,11]<-t.test(excessPSV[,i],mu=0,conf.level=0.95)$conf.int[1]
	output[i,12]<-t.test(excessPSV[,i],mu=0,conf.level=0.95)$conf.int[2]
	output[i,13]<-t.test(excessPSV[,i],mu=0,conf.level=0.95)$statistic
	output[i,14]<-t.test(excessPSV[,i],mu=0,conf.level=0.95)$p.value
	}

#Writing output to file
write.csv(output,file="Ecoregion_PSV_output.csv")

##################################################

##Fitting linear models for excess PD loss vs. richness and extinction rate with trees as random effect

#Loading essential package
library(lme4)

#Reading excess PD loss data (%) and output
excessPD_data<-read.csv("Ecoregion_PD_excessPD.csv",header=TRUE,row.names=1)*100
excessPD_output<-read.csv("Ecoregion_PD_output.csv",header=TRUE,row.names=1)

#Converting data array and subsetting dataset for at least 20 species
excessPD_array<-array(NA,dim=c(0,4))
colnames(excessPD_array)<-c("tree","richness","extinction","excessPD")
for (i in 1:length(excessPD_data[1,])){
	excessPD_array<-rbind(excessPD_array,cbind(1:1000,excessPD_output[i,1],excessPD_output[i,2]*100,excessPD_data[,i]))
	}
excessPD_all<-data.frame(excessPD_array)
excessPD_20<-excessPD_all[excessPD_all$richness>20,]

#Fitting linear models with mixed effects
excessPD_richness_20<-lmer(excessPD~richness+(1|tree),data=excessPD_20,REML=FALSE)
excessPD_extinction_20<-lmer(excessPD~extinction+(1|tree),data=excessPD_20,REML=FALSE)

#Fitting linear null model
excessPD_null_20<-lmer(excessPD~1+(1|tree),data=excessPD_20,REML=FALSE)

#Performing likelihood ratio test
anova(excessPD_null_20,excessPD_richness_20)
anova(excessPD_null_20,excessPD_extinction_20)

#Computing confidence intervals by bootstrapping
excessPD_richness_20_bs<-confint(excessPD_richness_20,level=0.95,method="boot",nsim=1000)
excessPD_extinction_20_bs<-confint(excessPD_extinction_20,level=0.95,method="boot",nsim=1000)

##################################################

##Fitting linear models for excess PSV loss vs. richness and extinction rate with trees as random effect

#Loading essential package
library(lme4)

#Reading excess PSV loss data (%) and output
excessPSV_data<-read.csv("Ecoregion_PSV_excessPSV.csv",header=TRUE,row.names=1)*100
excessPSV_output<-read.csv("Ecoregion_PSV_output.csv",header=TRUE,row.names=1)

#Converting data array and subsetting dataset for at least 20 species
excessPSV_array<-array(NA,dim=c(0,4))
colnames(excessPSV_array)<-c("tree","richness","extinction","excessPSV")
for (i in 1:length(excessPSV_data[1,])){
	excessPSV_array<-rbind(excessPSV_array,cbind(1:1000,excessPSV_output[i,1],excessPSV_output[i,2]*100,excessPSV_data[,i]))
	}
excessPSV_all<-data.frame(excessPSV_array)
excessPSV_20<-excessPSV_all[excessPSV_all$richness>20,]

#Fitting linear models with mixed effects
excessPSV_richness_20<-lmer(excessPSV~richness+(1|tree),data=excessPSV_20,REML=FALSE)
excessPSV_extinction_20<-lmer(excessPSV~extinction+(1|tree),data=excessPSV_20,REML=FALSE)

#Fitting linear null model
excessPSV_null_20<-lmer(excessPSV~1+(1|tree),data=excessPSV_20,REML=FALSE)

#Performing likelihood ratio test
anova(excessPSV_null_20,excessPSV_richness_20)
anova(excessPSV_null_20,excessPSV_extinction_20)

#Computing confidence intervals by bootstrapping
excessPSV_richness_20_bs<-confint(excessPSV_richness_20,level=0.95,method="boot",nsim=1000)
excessPSV_extinction_20_bs<-confint(excessPSV_extinction_20,level=0.95,method="boot",nsim=1000)
