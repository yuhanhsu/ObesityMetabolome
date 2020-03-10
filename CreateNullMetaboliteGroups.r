##########################################################################################
### Perform permutations across two datasets to generate null "cause", "effect", and 
### "bidirectional" metabolite groups for use in PAIRUP-MS pathway analysis
### (B=BMI, Gb=BMI GRS, M=metabolite, Gm=metabolite SNP)
### (1) permute (B, Gb) pairs in Dataset 1 and 2 -> (B', Gb') x nPerm
### (2) calculate B' ~ Gm and M ~ Gb' stats
### (3) calculate Wald ratios, SEs, and p-values using
###	(a) M ~ Gm and B' ~ Gm stats
###	(b) B' ~ Gb' (same as B ~ Gb) and M ~ Gb' stats
### (4) rank metabolites using Wald p-values from (3) to get:
###	(a) "cause": top quartile (most significant) of B' ~ Gm, bottom quartile of M ~ Gb'
###	(b) "effect": top quartile of M ~ Gb', bottom quartile of B' ~ Gm
###	(c) "bidirectional": top quartile of both
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(meta) # for meta-analysis


# read in input arguments (input and output file paths)
args <- commandArgs(trailingOnly=T)

realAssocFile <- args[1] # file containing observed association stats
bmiFile1 <- args[2] # BMI z-score data for Dataset 1
grsFile1 <- args[3] # BMI GRS values for Dataset 1
metFile1 <- args[4] # metabolite z-score data for Dataset 1
snpFile1 <- args[5] # metabolite SNP genotype dosage for Dataset 1
bmiFile2 <- args[6] # BMI z-score data for Dataset 2
grsFile2 <- args[7] # BMI GRS values for Dataset 2
metFile2 <- args[8] # metabolite z-score data for Dataset 2
snpFile2 <- args[9] # metabolite SNP genotype dosage for Dataset 2
nPerm <- as.numeric(args[10]) # number of permutations
bmiGrsBeta <- as.numeric(args[11]) # observed BMI ~ Gb association stats
bmiGrsSe <- as.numeric(args[12])

out_causeFile <- args[13] # file containing lists (rows) of cause metabolites
out_effectFile <- args[14] # file containing lists of effect metabolites
out_bidirFile <- args[15] # file containing lists of bidirectional metabolites


# -----------------------------------------------------------
### function for calculating Wald stats
calc.wald <- function(Bx,Bxse,By,Byse) {
	ratio <- By/Bx
	se <- sqrt(Byse^2 / Bx^2 + By^2 * Bxse^2 / Bx^4)
	p <- 2*pnorm(-abs(ratio/se)) # estimate p-value from normal distribution
	res <- c(ratio,se,p)
	return(res)
}


# -----------------------------------------------------------
# read in BMI-associated metabolites and correspodning Gm SNP ID (CHR:POS_A1_A2)
tempTable <- read.table(realAssocFile,header=T,sep="\t",stringsAsFactors=F)
tempTable$SNP_ID <- paste(tempTable$SNP,tempTable$A1,tempTable$A2,sep="_")
pairTable <- tempTable[,c("Metabolite1","Metabolite2","SNP_ID","METvSNP_Beta","METvSNP_SE")]
print(paste("# of BMI-associated metabolites",dim(pairTable)[1]))
print("")


# -----------------------------------------------------------
# read in Dataset 1 and Dataset 1 (BMI,GRS) values and
# perform paired permuations

# Dataset1
bmiGrsTable1 <- read.table(grsFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
print(paste("Dataset 1 samples with GRS:",dim(bmiGrsTable1)[1]))

tempTable <- read.table(bmiFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmiGrsTable1$BMI <- tempTable[match(rownames(bmiGrsTable1),rownames(tempTable)),"invBMI"]
print(paste("# of Dataset 1 BMI scores:",sum(!is.na(bmiGrsTable1$BMI))))

indices <- 1:dim(bmiGrsTable1)[1]
permTable1 <- bmiGrsTable1
for (i in 1:nPerm) {
	permTable <- bmiGrsTable1[sample(indices),]
	colnames(permTable) <- c(paste("GRS_",i,sep=""),paste("BMI_",i,sep=""))
	permTable1 <- cbind(permTable1,permTable)
}

# Dataset 2
bmiGrsTable2 <- read.table(grsFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
print(paste("Dataset 2 samples with GRS:",dim(bmiGrsTable2)[1]))

tempTable <- read.table(bmiFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmiGrsTable2$BMI <- tempTable[match(rownames(bmiGrsTable2),rownames(tempTable)),"invBMI"]
print(paste("# of OE BMI scores:",sum(!is.na(bmiGrsTable2$BMI))))

indices <- 1:dim(bmiGrsTable2)[1]
permTable2 <- bmiGrsTable2
for (i in 1:nPerm) {
	permTable <- bmiGrsTable2[sample(indices),]
	colnames(permTable) <- c(paste("GRS_",i,sep=""),paste("BMI_",i,sep=""))
	permTable2 <- cbind(permTable2,permTable)
}


# -----------------------------------------------------------
# calculate B' ~ Gm Wald p-values

tempTable <- read.table(snpFile1,sep="\t",header=T,check.names=F,row.names=1,stringsAsFactors=F)
snpTable1 <- tempTable[,match(pairTable$SNP_ID,colnames(tempTable))]
colnames(snpTable1) <- pairTable$Metabolite1
print(paste("Dataset 1 MET SNP dosage data:",dim(snpTable1)[1],"samples x",dim(snpTable1)[2],"mets"))

tempTable <- read.table(snpFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
snpTable2 <- tempTable[,match(pairTable$SNP_ID,colnames(tempTable))]
colnames(snpTable2) <- pairTable$Metabolite1
print(paste("Dataset 2 MET SNP dosage data:",dim(snpTable2)[1],"samples x",dim(snpTable2)[2],"mets"))

permBmiSnpPTable <- data.frame(Metabolite1=pairTable$Metabolite1)
for (i in 1:nPerm) {
	pvalues <- NULL
	for (j in 1:dim(pairTable)[1]) {	
		# B' ~ Gm
		reg1 <- summary(lm(permTable1[,paste("BMI_",i,sep="")] ~ snpTable1[,pairTable$Metabolite1[j]]))
		reg2 <- summary(lm(permTable2[,paste("BMI_",i,sep="")] ~ snpTable2[,pairTable$Metabolite1[j]]))
		res <- metagen(c(coef(reg1)[2,"Estimate"],coef(reg2)[2,"Estimate"]),
			c(coef(reg1)[2,"Std. Error"],coef(reg2)[2,"Std. Error"]))
	
		# Wald stats
		waldRes <- calc.wald(pairTable$METvSNP_Beta[j],pairTable$METvSNP_SE[j],
				res$TE.fixed,res$seTE.fixed) # Bx, Bxse, By, Byse

		# store Wald P
		pvalues <- c(pvalues,waldRes[3])
	}
	permBmiSnpPTable[,paste("Perm_",i,sep="")] <- pvalues

	if (i %% 100 == 0) { print(i) }
}


# -----------------------------------------------------------
# calculate M ~ Gb' Wald p-values

tempTable <- read.table(metFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metTable1 <- tempTable[match(rownames(permTable1),rownames(tempTable)),match(pairTable$Metabolite1,colnames(tempTable))]
print(paste("Dataset 1 metabolite data:",dim(metTable1)[1],"samples x",dim(metTable1)[2],"mets"))

tempTable <- read.table(metFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metTable2 <- tempTable[match(rownames(permTable2),rownames(tempTable)),match(pairTable$Metabolite2,colnames(tempTable))]
print(paste("Dataset 2 metabolite data:",dim(metTable2)[1],"samples x",dim(metTable2)[2],"mets"))

permMetGrsPTable <- data.frame(Metabolite1=pairTable$Metabolite1)
for (i in 1:nPerm) {
	pvalues <- NULL
	for (j in 1:dim(pairTable)[1]) {
		# M ~ Gb'
		reg1 <- summary(lm(metTable1[,pairTable$Metabolite1[j]] ~ permTable1[,paste("GRS_",i,sep="")]))
		reg2 <- summary(lm(metTable2[,pairTable$Metabolite2[j]] ~ permTable2[,paste("GRS_",i,sep="")]))
		res <- metagen(c(coef(reg1)[2,"Estimate"],coef(reg2)[2,"Estimate"]),
			c(coef(reg1)[2,"Std. Error"],coef(reg2)[2,"Std. Error"]))
			
		# Wald stats
		waldRes <- calc.wald(bmiGrsBeta,bmiGrsSe,res$TE.fixed,res$seTE.fixed) # Bx, Bxse, By, Byse
		
		# store Wald P
		pvalues <- c(pvalues,waldRes[3])
	}
	permMetGrsPTable[,paste("Perm_",i,sep="")] <- pvalues

	if (i %% 100 == 0) { print(i) }
}


# -----------------------------------------------------------
# define null cause, effect, and bidirectional mets using quartile cutoffs
# output met lists (one row per permutation)

causeLists <- NULL
effectLists <- NULL
bidirLists <- NULL
for (i in 1:nPerm) {
	tempX <- -log10(permBmiSnpPTable[,paste("Perm_",i,sep="")])
	tempY <- -log10(permMetGrsPTable[,paste("Perm_",i,sep="")])

	causeLists <- rbind(causeLists,
		paste(pairTable$Metabolite1[tempX > quantile(tempX)[4] & tempY < quantile(tempY)[2]],collapse=";"))
	effectLists <- rbind(effectLists,
		paste(pairTable$Metabolite1[tempX < quantile(tempX)[2] & tempY > quantile(tempY)[4]],collapse=";"))
	bidirLists <- rbind(bidirLists,
		paste(pairTable$Metabolite1[tempX > quantile(tempX)[4] & tempY > quantile(tempY)[4]],collapse=";"))		
}

write.table(causeLists,out_causeFile,col.names=F,row.names=F,quote=F)
write.table(effectLists,out_effectFile,col.names=F,row.names=F,quote=F)
write.table(bidirLists,out_bidirFile,col.names=F,row.names=F,quote=F)

