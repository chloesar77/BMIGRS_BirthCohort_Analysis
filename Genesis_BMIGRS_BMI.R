#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(GWASTools)
library(GENESIS)
library(SNPRelate)
library(gdsfmt)
library(Rcpp)
library(data.table)


#########
#load the phenotype file
#########

pheno<-fread(args[1])

pheno2=pheno[,c("framid","idtype.x","cohort","sex.x","age50","BMI_age50","mean_birth","birthcohort3","birthcohort4","wGRS_941","famid","id","fa","mo","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","SMK")]
pheno2$wGRS=pheno2$wGRS_941
pheno2$logBMI=log(pheno2$BMI_age50)
#winsorization approach
pheno2$logBMI[pheno2$logBMI>=4]=4
#standardization of the GRS
pheno2$wGRS=pheno2$wGRS/sd(pheno2$wGRS)

pheno2$mean_birth=pheno2$mean_birth/1000
pheno2$age2_50=pheno2$age50*pheno2$age50
pheno2$sex.x[pheno2$sex.x==1]="M"
pheno2$sex.x[pheno2$sex.x==2]="F"

scanAnnot<-ScanAnnotationDataFrame(data.frame(scanID=pheno2$framid,idtype=pheno2$idtype.x,cohort=as.factor(pheno2$cohort),sex=pheno2$sex.x,age50=as.numeric(pheno2$age50),age2_50=as.numeric(pheno2$age2_50),BMI_age50=as.numeric(pheno2$BMI_age50),mean_birth=as.numeric(pheno2$mean_birth),birthcohort3=as.factor(pheno2$birthcohort3),birthcohort4=as.factor(pheno2$birthcohort4),wGRS=as.numeric(pheno2$wGRS),famid=pheno2$famid,id=pheno2$id,fa=pheno2$fa,mo=pheno2$mo,PC1=as.numeric(pheno2$PC1),PC2=as.numeric(pheno2$PC2),PC3=as.numeric(pheno2$PC3),PC4=as.numeric(pheno2$PC4),PC5=as.numeric(pheno2$PC5),PC6=as.numeric(pheno2$PC6),PC7=as.numeric(pheno2$PC7),PC8=as.numeric(pheno2$PC8),PC9=as.numeric(pheno2$PC9),PC10=as.numeric(pheno2$PC10),logBMI=as.numeric(pheno2$logBMI),SMK=as.factor(pheno2$SMK)))

library('MASS')
library('kinship2')
library('coxme')

kmat=makekinship(famid=pheno2$famid, id=pheno2$framid, mother.id=pheno2$mo, father.id=pheno2$fa)

#########
#load grm 
#########

usemypcrel<-kmat
usemypcrel<-as.matrix(usemypcrel)

#########
#NULL model
#########
covMatList <- list("Kin" = usemypcrel)
mat<-as.matrix(usemypcrel)

#the model includes the BMI GRS, birth year (continuous) / birth cohort in 4 categories / or FHS cohorts (Gen1/2/3) and the interaction

nullmod <- fitNullModel(scanAnnot, family = "gaussian", outcome = "logBMI", covars=c("sex","age50","SMK","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS","cohort","wGRS*cohort"), group.var="cohort", cov.mat =mat)

nullmod1 <- fitNullModel(scanAnnot, family = "gaussian", outcome = "logBMI", covars=c("sex","age50","SMK","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS","mean_birth","wGRS*mean_birth"), group.var="cohort", cov.mat =mat)

nullmod3 <- fitNullModel(scanAnnot, family = "gaussian", outcome = "logBMI", covars=c("sex","age50","SMK","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS","birthcohort4","wGRS*birthcohort4"), group.var="birthcohort4", cov.mat =mat)

outnullmod<-as.data.frame(nullmod$fixef)
outnullmod1<-as.data.frame(nullmod1$fixef)
outnullmod3<-as.data.frame(nullmod3$fixef)

write.csv(outnullmod,args[2])
write.csv(outnullmod1,args[3])
write.csv(outnullmod3,args[4])

