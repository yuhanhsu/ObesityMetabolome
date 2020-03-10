##########################################################################################
### Calculate BMI effect size-weighted genetic risk score (GRS) for samples 
### with imputed genotype data
###
### author: Yu-Han Hsu
##########################################################################################

import sys


# read in input arguments (input and output file paths)
bmiSnpFile = sys.argv[1] # file containing BMI SNPs with allele and effect size info
sampleFile = sys.argv[2] # file containing sample IDs
genotypeFile = sys.argv[3] # file containing genotype data for each BMI SNP and each sample

out_GrsFile = sys.argv[4] # output file containing GRS values for each sample


# read in BMI SNPs and corresonding effect size estimates (BETAs)
snpDict = {} # CHR:POS as key, [BETA, BMI Allele, Other Allele] as value
with open(bmiSnpFile,'r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.split()
		snpDict[fields[1]+':'+fields[2]] = [fields[5],fields[3],fields[4]]

print len(snpDict)


# read in list of sample IDs (in order they appear in genotype data file below)
sampleList = []
with open(sampleFile,'r') as inFile:
	for line in inFile:
		sampleList.append(line.split()[0])

print len(sampleList)


# calculate dosage for BMI allele for each sample
# calculate effect size-weighted genetic risk score
with open(genotypeFile,'r') as inFile:
	GRSs = [0] * len(sampleList)
	
	for line in inFile:
		fields = line.strip().split()
		chrPos = fields[0] + ':' + fields[2]
		beta,effect,other = snpDict[chrPos]
		
		dosages = []
		for i in range(len(sampleList)):
			# missing genotype
			if (fields[5+i*3] == '0' and fields[6+i*3] == '0'
				and fields[7+i*3] == '0'):
				dosages.append('NA')
			
			# second allele = BMI allele
			elif fields[4] == effect and fields[3] == other:
				dosage = '%g' % (float(fields[6+i*3]) * 1
					+ float(fields[7+i*3]) * 2)
				dosages.append(dosage)
			
			# first allele = BMI allele
			elif fields[3] == effect and fields[4] == other:
				dosage = '%g' % (float(fields[5+i*3]) * 2
					+ float(fields[6+i*3]) * 1)
				dosages.append(dosage)
			
			# the two alleles don't correspond to the stored BMI SNP alleles
			else:
				print 'ERROR: mismatched alleles'

		if 'NA' not in dosages:
			GRSs = [y+z for y,z in zip(GRSs,[float(beta) * float(x) for x in dosages])]


# output GRS values for each sample
with open(out_GrsFile,'w') as outFile:
	outFile.write('Sample\tGRS\n')

	for i in range(len(sampleList)):
		outFile.write(sampleList[i] + '\t' + str(GRSs[i]) + '\n')


print 'SCRIPT COMPLETED'
