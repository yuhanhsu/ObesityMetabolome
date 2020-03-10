##########################################################################################
### Calculate association between GB (BMI GRS) and BMI or BMI-associated metabolites in 
### two datasets using linear regression and inverse variance-weighted meta-analysis
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(meta) # for meta-analysis


# read in input arguments (input and output file paths)
args <- commandArgs(trailingOnly=T)

pairFile <- args[1] # mapping of BMI-associated metabolites across two datasets
grsFile1 <- args[2] # BMI GRS values for Dataset 1
metFile1 <- args[3] # metabolite z-score data for Dataset 1
ptypeFile1 <- args[4] # BMI z-score data for Dataset 1
grsFile2 <- args[5] # BMI GRS values for Dataset 2
metFile2 <- args[6] # metabolite z-score data for Dataset 2
ptypeFile2 <- args[7] # BMI z-score data for Dataset 2

out_assocFile <- args[8] # BMI-GRS and metabolite-GRS association statistics in Dataset 1, Dataset 2, and their meta-analysis


# -----------------------------------------------------------
# read in BMI-associated metabolite pairs
pairTable <- read.table(pairFile,header=T,sep="\t",stringsAsFactors=F) # "Metabolite1","Metabolite2"
print(paste("# of BMI-associated metabolites:",dim(pairTable)[1]))


# -----------------------------------------------------------
# calculate GRS association stats in Dataset 1

grs1 <- read.table(grsFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
print(paste("Dataset 1 samples with GRS:",dim(grs1)[1]))

metData1 <- read.table(metFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metData1 <- metData1[match(rownames(grs1),rownames(metData1)),match(pairTable$Metabolite1,colnames(metData1))]
print(paste("Dataset 1 metabolite data:",dim(metData1)[1],"samples x",dim(metData1)[2],"metabolites"))

tempTable <- read.table(ptypeFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmi1 <- tempTable[match(rownames(metData1),rownames(tempTable)),"invBMI"]
print(paste("# of Dataset 1 BMI scores:",sum(!is.na(bmi1))))

metData1 <- cbind(bmi1,metData1)
colnames(metData1)[1] <- "invBMI"

statsTable1 <- NULL
for (m in colnames(metData1)) {
	reg <- summary(lm(metData1[,m] ~ grs1$GRS))
	df <- data.frame(Variable1=m,Beta1=coef(reg)[2,"Estimate"],
		SE1=coef(reg)[2,"Std. Error"],P1=coef(reg)[2,"Pr(>|t|)"])
	statsTable1 <- rbind(statsTable1,df)
}


# -----------------------------------------------------------
# calculate GRS association stats in Dataset 2

grs2 <- read.table(grsFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
print(paste("Dataset 2 samples with GRS:",dim(grs2)[1]))

metData2 <- read.table(metFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metData2 <- metData2[match(rownames(grs2),rownames(metData2)),match(pairTable$Metabolite2,colnames(metData2))]
print(paste("Dataset 2 metabolite data:",dim(metData2)[1],"samples x",dim(metData2)[2],"metabolites"))

tempTable <- read.table(ptypeFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmi2 <- tempTable[match(rownames(metData2),rownames(tempTable)),"invBMI"]
print(paste("# of Dataset 2 BMI scores:",sum(!is.na(bmi2))))

metData2 <- cbind(bmi2,metData2)
colnames(metData2)[1] <- "invBMI"

statsTable2 <- NULL
for (m in colnames(metData2)) {
	reg <- summary(lm(metData2[,m] ~ grs2$GRS))
	df <- data.frame(Variable2=m,Beta2=coef(reg)[2,"Estimate"],
		SE2=coef(reg)[2,"Std. Error"],P2=coef(reg)[2,"Pr(>|t|)"])#,MCDS_Rsquared=reg$r.squared)
	statsTable2 <- rbind(statsTable2,df)
}


# -----------------------------------------------------------
# perform inverse variance-weighted meta-analysis of GRS association stats

allStatsTable <- NULL
for (i in 1:dim(statsTable1)[1]) {
	res <- metagen(c(statsTable1[i,"Beta1"],statsTable2[i,"Beta2"]),
		c(statsTable1[i,"SE1"],statsTable2[i,"SE2"]))
	df <- data.frame(META_Beta=res$TE.fixed,META_SE=res$seTE.fixed,META_P=res$pval.fixed)
	allStatsTable <- rbind(allStatsTable,df)
}

# output all assoication stats
allStatsTable <- cbind(statsTable1,statsTable2,allStatsTable)
write.table(allStatsTable,out_assocFile,sep="\t",row.names=F,quote=F)


print("SCRIPT COMPLETED")
