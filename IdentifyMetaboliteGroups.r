##########################################################################################
### (1) Generate plot of -log10 GB IV p-value vs. -log10 GM IV p-value for 
###     BMI-associated metabolites
### (2) Identify "cause", "effect", and "bidirectional" metabolite groups using 
###     top and bottom quartile cutoffs of GB and GM IV p-values
###
### author: Yu-Han Hsu
##########################################################################################

rm(list=ls())

library(ggplot2)
library(ggrepel)

# read in input arguments (input and output file paths)
args <- commandArgs(trailingOnly=T)

ivStatsFile <- args[1] # file containing GM and GB Wald IV p-values

out_CauseFile <- args[2] # list of cause metabolites (use names in Metabolite1 column)
out_EffectFile <- args[3] # list of effect metabolites
out_BidirFile <- args[4] # list of bidirectional metabolites
out_PlotFile <- args[5] # plot of -log10(Wald IV p-value) of GB vs. GM


# -----------------------------------------------------------
# read in Wald IV stats
ivStatsTable <- read.table(ivStatsFile,header=T,sep="\t",stringsAsFactors=F)
dim(ivStatsTable)


# -----------------------------------------------------------
# plot of -log10 GB IV p-value vs. -log10 GM IV p-value

# plot shared knowns last (clearer view in plot)
plotTable <- rbind(ivStatsTable[ivStatsTable$Type=="Matched",],ivStatsTable[ivStatsTable$Type=="Shared_Known",])
plotTable$Type <- factor(plotTable$Type,c("Shared_Known","Matched")) # specify order

# define shared known mets to label in the 3 metabolite groups 
tempX <- -log10(plotTable$BMIvSNP_WaldP)
tempY <- -log10(plotTable$METvGRS_WaldP)
isCause <- tempX > quantile(tempX)[4] & tempY < quantile(tempY)[2]
isEffect <- tempX < quantile(tempX)[2] & tempY > quantile(tempY)[4]
isBidir <- tempX > quantile(tempX)[4] & tempY > quantile(tempY)[4]
labelTable <- subset(plotTable,(isCause | isEffect | isBidir) & Type=="Shared_Known")

pdf(out_PlotFile,width=5,height=5.2)

# plot with no metaoblite labels
ggplot(plotTable,aes(x=-log10(BMIvSNP_WaldP),y=-log10(METvGRS_WaldP),color=Type,shape=Type)) +
	geom_rect(aes(xmin=quantile(-log10(BMIvSNP_WaldP))[4],xmax=Inf,ymin=-Inf,ymax=quantile(-log10(METvGRS_WaldP))[2]),
		fill="#FFCC99",alpha=0.1,inherit.aes=F) + # orange shade = cause region
	geom_rect(aes(xmin=-Inf,xmax=quantile(-log10(BMIvSNP_WaldP))[2],ymin=quantile(-log10(METvGRS_WaldP))[4],ymax=Inf),
		fill="#99CCFF",alpha=0.1,inherit.aes=F) + # blue shade = effect region
	geom_rect(aes(xmin=quantile(-log10(BMIvSNP_WaldP))[4],xmax=Inf,ymin=quantile(-log10(METvGRS_WaldP))[4],ymax=Inf),
		fill="#FFCCFF",alpha=0.1,inherit.aes=F) + # pink shade = bidirectional region
	
	geom_point(alpha=0.9) + xlab("-log10 BMI-SNP P") + ylab("-log10 MET-GRS P") +
	scale_colour_manual(values=c("red","grey50"),name="Metabolite Type",labels=c("Shared Known","Matched")) +
	scale_shape_manual(values=c(17,16),name="Metabolite Type",labels=c("Shared Known","Matched")) +
	
	geom_vline(aes(xintercept=quantile(-log10(BMIvSNP_WaldP))[4]),color="grey40",lty=2) +
	geom_hline(aes(yintercept=quantile(-log10(METvGRS_WaldP))[2]),color="grey40",lty=2) +
	geom_vline(aes(xintercept=quantile(-log10(BMIvSNP_WaldP))[2]),color="grey40",lty=2) +
	geom_hline(aes(yintercept=quantile(-log10(METvGRS_WaldP))[4]),color="grey40",lty=2) + 
	
	theme_bw() + theme(legend.position="top")

# plot with shared known labels in the three metaoblite groups
ggplot(plotTable,aes(x=-log10(BMIvSNP_WaldP),y=-log10(METvGRS_WaldP),color=Type,shape=Type)) +
	geom_rect(aes(xmin=quantile(-log10(BMIvSNP_WaldP))[4],xmax=Inf,ymin=-Inf,ymax=quantile(-log10(METvGRS_WaldP))[2]),
		fill="#FFCC99",alpha=0.1,inherit.aes=F) + # orange shade = cause region
	geom_rect(aes(xmin=-Inf,xmax=quantile(-log10(BMIvSNP_WaldP))[2],ymin=quantile(-log10(METvGRS_WaldP))[4],ymax=Inf),
		fill="#99CCFF",alpha=0.1,inherit.aes=F) + # blue shade = effect region
	geom_rect(aes(xmin=quantile(-log10(BMIvSNP_WaldP))[4],xmax=Inf,ymin=quantile(-log10(METvGRS_WaldP))[4],ymax=Inf),
		fill="#FFCCFF",alpha=0.1,inherit.aes=F) + # pink shade = bidirectional region
	
	geom_point(alpha=0.9) + xlab("-log10 BMI-SNP P") + ylab("-log10 MET-GRS P") +
	scale_colour_manual(values=c("red","grey50"),name="Metabolite Type",labels=c("Shared Known","Matched")) +
	scale_shape_manual(values=c(17,16),name="Metabolite Type",labels=c("Shared Known","Matched")) +
	
	geom_vline(aes(xintercept=quantile(-log10(BMIvSNP_WaldP))[4]),color="grey40",lty=2) +
	geom_hline(aes(yintercept=quantile(-log10(METvGRS_WaldP))[2]),color="grey40",lty=2) +
	geom_vline(aes(xintercept=quantile(-log10(BMIvSNP_WaldP))[2]),color="grey40",lty=2) +
	geom_hline(aes(yintercept=quantile(-log10(METvGRS_WaldP))[4]),color="grey40",lty=2) + 
	
	geom_text_repel(data=labelTable,aes(x=-log10(BMIvSNP_WaldP),y=-log10(METvGRS_WaldP),label=Metabolite1),
		size=2.8,fontface="bold",inherit.aes=F) +
	
	theme_bw() + theme(legend.position="top")

dev.off()


# -----------------------------------------------------------
# output list of cause, effect, and bidirectional metabolites

tempX <- -log10(ivStatsTable$BMIvSNP_WaldP)
tempY <- -log10(ivStatsTable$METvGRS_WaldP)

isCause <- tempX > quantile(tempX)[4] & tempY < quantile(tempY)[2]
isEffect <- tempX < quantile(tempX)[2] & tempY > quantile(tempY)[4]
isBidir <- tempX > quantile(tempX)[4] & tempY > quantile(tempY)[4]

print("CAUSE metabolites:")
print(paste("total:",sum(isCause)))
print(paste("shared known:",sum(ivStatsTable$Type[isCause]=="Shared_Known")))
print(paste("matched:",sum(ivStatsTable$Type[isCause]=="Matched")))

write(paste(ivStatsTable$Metabolite1[isCause],collapse=";"),out_CauseFile)


print("EFFECT metabolites:")
print(paste("total:",sum(isEffect)))
print(paste("shared known:",sum(ivStatsTable$Type[isEffect]=="Shared_Known")))
print(paste("matched:",sum(ivStatsTable$Type[isEffect]=="Matched")))

write(paste(ivStatsTable$Metabolite1[isEffect],collapse=";"),out_EffectFile)


print("BIDIRECTIONAL metabolites:")
print(paste("total:",sum(isBidir)))
print(paste("shared known:",sum(ivStatsTable$Type[isBidir]=="Shared_Known")))
print(paste("matched:",sum(ivStatsTable$Type[isBidir]=="Matched")))

write(paste(ivStatsTable$Metabolite1[isBidir],collapse=";"),out_BidirFile)


print("SCRIPT COMPLETED")
