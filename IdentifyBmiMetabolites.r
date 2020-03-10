##########################################################################################
### (1) Calculate association between BMI and shared known and matched metabolites 
###     across two datasets
### (2) Generate clustered correlation heat map of the BMI-associated metabolites
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(meta) # for meta-analysis
library(RColorBrewer) # for heat map
library(gplots)


# read in input arguments (input and output file paths)
args <- commandArgs(trailingOnly=T)

knownFile <- args[1] # mapping of shared known metabolites across the two datasets
matchFile <- args[2] # mapping of matched metabolites across the two datasets
metFile1 <- args[3] # metabolite z-score data for Dataset 1
ptypeFile1 <- args[4] # BMI z-score data for Dataset 1
metFile2 <- args[5] # metabolite z-score data for Dataset 2
ptypeFile2 <- args[6] # BMI z-score data for Dataset 2

out_assocFile <- args[7] # metabolite-BMI association statistics in Dataset 1, Dataset 2, and their meta-analysis
out_plotFile <- args[8] # clutstered correlation heat map of metabolites associated with BMI at Bonferroni significance in the meta-analysis


# -----------------------------------------------------------
# shared known and matched metabolites to include in analysis

# read in names of shared known metabolites in the two datasets
knownTable <- read.table(knownFile,header=T,sep="\t",stringsAsFactors=F) # "Metabolite1","Metabolite2"
knownTable <- knownTable[complete.cases(knownTable),] # skip entries with "NA"
print(paste("# of shared knowns:",dim(knownTable)[1]))

# read in names of matched metabolites in the two datasets
matchTable <- read.table(matchFile,header=T,sep="\t",stringsAsFactors=F)[,c("Metabolite1","Metabolite2")]
# if a metaoblite (in Dataset 1) appears in both shared known and matched tables, remove the matched entry
matchTable <- matchTable[!(matchTable$Metabolite1 %in% knownTable$Metabolite1),]
print(paste("# of matched metabolites:",dim(matchTable)[1]))

# create table of shared knowns + matched metabolites
allTable <- cbind(data.frame(Type=c(rep("Shared_Known",dim(knownTable)[1]),rep("Matched",dim(matchTable)[1]))),
	rbind(knownTable,matchTable)) 
print(paste("total # of metabolites:",dim(allTable)[1]))
print("")


# -----------------------------------------------------------
# calculate BMI association stats for Dataset 1 metabolites

metData1 <- read.table(metFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metData1 <- metData1[,match(allTable$Metabolite1,colnames(metData1))]
print(paste("Dataset 1 metabolite data:",dim(metData1)[1],"samples x",dim(metData1)[2],"metabolites")) 

tempTable <- read.table(ptypeFile1,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmi1 <- tempTable[match(rownames(metData1),rownames(tempTable)),"invBMI"]
print(paste("# of Dataset 1 BMI scores:",sum(!is.na(bmi1))))

statsTable1 <- NULL
for (m in allTable$Metabolite1) {
	reg <- summary(lm(bmi1 ~ metData1[,m]))
	df <- data.frame(Metabolite1=m,Beta1=coef(reg)[2,"Estimate"],
		SE1=coef(reg)[2,"Std. Error"],P1=coef(reg)[2,"Pr(>|t|)"])
	statsTable1 <- rbind(statsTable1,df)
}


# -----------------------------------------------------------
# calculate BMI association stats for Dataset 2 metabolites

metData2 <- read.table(metFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
metData2 <- metData2[,match(allTable$Metabolite2,colnames(metData2))]
print(paste("Dataset 2 metabolite data:",dim(metData2)[1],"samples x",dim(metData2)[2],"metabolites"))

tempTable <- read.table(ptypeFile2,header=T,sep="\t",check.names=F,row.names=1,stringsAsFactors=F)
bmi2 <- tempTable[match(rownames(metData2),rownames(tempTable)),"invBMI"]
print(paste("# of Dataset 2 BMI scores:",sum(!is.na(bmi2))))

statsTable2 <- NULL
for (m in allTable$Metabolite2) {
	reg <- summary(lm(bmi2 ~ metData2[,m]))
	df <- data.frame(Metabolite2=m,Beta2=coef(reg)[2,"Estimate"],
		SE2=coef(reg)[2,"Std. Error"],P2=coef(reg)[2,"Pr(>|t|)"])
	statsTable2 <- rbind(statsTable2,df)
}


# -----------------------------------------------------------
# perform inverse variance-weighted meta-analysis of BMI association stats

allStatsTable <- NULL
for (i in 1:dim(allTable)[1]) {
	res <- metagen(c(statsTable1[i,"Beta1"],statsTable2[i,"Beta2"]),
		c(statsTable1[i,"SE1"],statsTable2[i,"SE2"]))
	df <- data.frame(META_Beta=res$TE.fixed,META_SE=res$seTE.fixed,META_P=res$pval.fixed)
	allStatsTable <- rbind(allStatsTable,df)
}


# output all assoication stats
allStatsTable <- cbind(allTable$Type,statsTable1,statsTable2,allStatsTable)
colnames(allStatsTable)[1] <- "Type"

write.table(allStatsTable,out_assocFile,sep="\t",row.names=F,quote=F)

print("")


# -----------------------------------------------------------
# identify BMI-associated metabolites using Bonferroni-corrected p-value threshold

bonfThres <- 0.05/dim(allTable)[1] 
print(paste("Bonferroni p-value:",bonfThres))

sigMetaSignals <- allStatsTable$Metabolite1[allStatsTable$META_P < bonfThres]
print(paste("# of signficant metabolites in meta-analysis:",length(sigMetaSignals)))
print(paste("   # of shared knowns:",sum(sigMetaSignals %in% knownTable$Metabolite1)))
print(paste("   # of matched metabolites:",sum(sigMetaSignals %in% matchTable$Metabolite1)))
print("")

sigMetaKnown <- sigMetaSignals[sigMetaSignals %in% knownTable$Metabolite1]
sigMetaMatched <- sigMetaSignals[sigMetaSignals %in% matchTable$Metabolite1]


# -----------------------------------------------------------
# generate clustered correlation heat map of BMI-associated metabolites

tempData <- metData2[,match(allTable$Metabolite2[match(sigMetaSignals,allTable$Metabolite1)],colnames(metData2))]
colnames(tempData) <- allTable$Metabolite1[match(sigMetaSignals,allTable$Metabolite1)]

allBmiMetData <- rbind(metData1[,sigMetaSignals],tempData)
allBmiMetCorMatrix <- cor(allBmiMetData) # correlation matrix

# hierarchical clustering
hr <- hclust(as.dist(1-allBmiMetCorMatrix),method="complete")

# for side bar indicating signal type
temp <- colnames(allBmiMetData) %in% knownTable$Metabolite1
typeColors <- rep("grey",length(sigMetaSignals))
typeColors[temp] <- "black" # black = shared known, grey = matched

# generate heat map PDF
pdf(out_plotFile,width=25,height=25)
heatmap.2(allBmiMetCorMatrix,Rowv=as.dendrogram(hr),Colv=as.dendrogram(hr),col=rev(brewer.pal(9,"RdBu")),
	trace="none",dendrogram="both",labCol=NA,margins=c(0.1,5),cexRow=0.4,
	RowSideColors=typeColors, # color side bar
	density.info="none", key.title=NA, key.xlab=NA, key.ylab=NA,symkey=F,symbreaks=F,
	key.par=list(mar=c(3,4,2,10)), #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
	lmat=rbind(c(0,0,4),c(3,1,2),c(0,0,5)), lhei=c(0.3,4,0.2), lwid=c(0.3,0.1,4))
	# lmat: 1=row side bar, 2=heatmap, 3=row dendrogram, 4=column dendrogram, 5=key
	
	#lmat=rbind(c(0,3),c(2,1),c(0,4)), lhei=c(0.3,4,0.2), lwid=c(0.5,4))
	# lmat: 1=heatmap, 2=row dendrogram, 3=column dendrogram, 4=key
dev.off()


print("SCRIPT COMPLETED")
