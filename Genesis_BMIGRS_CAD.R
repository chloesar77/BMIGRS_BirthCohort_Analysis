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
colnames(pheno)[2]="idtype.x"

pheno2=pheno[,c("framid","idtype.x","cohort","sex.x","age65","SMK","CAD_event","mean_birth","birthcohort3","birthcohort4","birthcohort2groups","wGRS_941","famid","id","fa","mo","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")]
pheno2$wGRS.x=pheno2$wGRS_941
#standardization of the GRS
pheno2$wGRS.x=pheno2$wGRS.x/sd(pheno2$wGRS.x)

pheno2$mean_birth=pheno2$mean_birth/1000

pheno2$sex.x[pheno2$sex.x==1]="M"
pheno2$sex.x[pheno2$sex.x==2]="F"

scanAnnot<-ScanAnnotationDataFrame(data.frame(scanID=pheno2$framid,idtype=pheno2$idtype.x,SMK=as.factor(pheno2$SMK),cohort=as.factor(pheno2$cohort),sex=pheno2$sex.x,age65=as.numeric(pheno2$age65),CAD_event=as.numeric(pheno2$CAD_event),mean_birth=as.numeric(pheno2$mean_birth),birthcohort3=as.factor(pheno2$birthcohort3),birthcohort4=as.factor(pheno2$birthcohort4),birthcohort2groups=as.numeric(pheno2$birthcohort2groups),wGRS.x=as.numeric(pheno2$wGRS.x),famid=pheno2$famid,id=pheno2$id,fa=pheno2$fa,mo=pheno2$mo,PC1=as.numeric(pheno2$PC1),PC2=as.numeric(pheno2$PC2),PC3=as.numeric(pheno2$PC3),PC4=as.numeric(pheno2$PC4),PC5=as.numeric(pheno2$PC5),PC6=as.numeric(pheno2$PC6),PC7=as.numeric(pheno2$PC7),PC8=as.numeric(pheno2$PC8),PC9=as.numeric(pheno2$PC9),PC10=as.numeric(pheno2$PC10)))

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

#the model includes the BMI GRS, birth year (continuous) / birth cohort in 2 categories / or FHS cohorts (Gen1/2/3) and the interaction

nullmod <- fitNullModel(scanAnnot, family = "binomial", outcome = "CAD_event", covars=c("sex","SMK","age65","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS.x","cohort","wGRS.x*cohort"), cov.mat =mat)

nullmod1 <- fitNullModel(scanAnnot, family = "binomial", outcome = "CAD_event", covars=c("sex","SMK","age65","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS.x","mean_birth","wGRS.x*mean_birth"), cov.mat =mat)

nullmod4 <- fitNullModel(scanAnnot, family = "binomial", outcome = "CAD_event", covars=c("sex","SMK","age65","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","wGRS.x","birthcohort2groups","wGRS.x*birthcohort2groups"), cov.mat =mat)

outnullmod<-as.data.frame(nullmod$fixef)
outnullmod1<-as.data.frame(nullmod1$fixef)
outnullmod4<-as.data.frame(nullmod4$fixef)

write.csv(outnullmod,args[2])
write.csv(outnullmod1,args[3])
write.csv(outnullmod4,args[4])


