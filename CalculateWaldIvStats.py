##########################################################################################
### Calculate Wald ratio IV effect estimate and corresponding SE and p-value
### for metabolite and BMI instruments
###
### author: Yu-Han Hsu
##########################################################################################

import sys, math
import scipy.stats


# read in input arguments (input and output file paths)
assocFile = sys.argv[1] # file containing obsevational association statistics betweeen MET, MET SNP, BMI, and BMI GRS
bmiGrsBeta = float(sys.argv[2]) # BMI and BMI instrument (GRS) association
bmiGrsSe = float(sys.argv[3])

out_WaldFile = sys.argv[4] # output file containing observational association + Wald IV statistics


### Wald IV statistic function
def calc_wald(Bx,Bxse,By,Byse):
	ratio = By/Bx # Wald ratio
	se = math.sqrt(Byse**2 / Bx**2 + By**2 * Bxse**2 / Bx**4) # standard error of Wald ratio
	p = 2*scipy.stats.norm.cdf(-abs(ratio/se)) # estimate p-value from normal distribution
	return ratio, se, p


# read in observational association statistics, output with Wald IV effect estimate statistics
with open(assocFile,'r') as inFile, open(out_WaldFile,'w') as outFile:

	headers = next(inFile).strip().split()

	# indices for getting relevant BETA and SE values
	ms_b = headers.index('METvSNP_Beta') # metabolite and metabolite instrument association from GWAS
	ms_se = headers.index('METvSNP_SE')
	bs_b = headers.index('BMIvSNP_Beta') # BMI and metabolite instrument association from GWAS
	bs_se = headers.index('BMIvSNP_SE')
	
	#bmiGrsBeta = 1.179 # BMI and BMI instrument (GRS) association
	#bmiGrsSe = 0.234
	mg_b = headers.index('METvGRS_Beta') # metabolite and BMI instrument (GRS) association
	mg_se = headers.index('METvGRS_SE')

	# output file headers
	outHeaders = headers + ['BMIvSNP_Wald','BMIvSNP_WaldSE','BMIvSNP_WaldP','METvGRS_Wald','METvGRS_WaldSE','METvGRS_WaldP']
	outFile.write('\t'.join(outHeaders) + '\n')

	for line in inFile:
		fields = line.strip().split()

		# BMIvSNP Wald stats (Gm -> M -> B)
		metSnpBeta = float(fields[ms_b])
		metSnpSe = float(fields[ms_se])
		bmiSnpBeta = float(fields[bs_b])
		bmiSnpSe = float(fields[bs_se])

		bmiSnpWald,bmiSnpWaldSe,bmiSnpWaldP = calc_wald(metSnpBeta,metSnpSe,bmiSnpBeta,bmiSnpSe)
		
		# METvGRS Wald stats (Gb -> B -> M)
		metGrsBeta = float(fields[mg_b])
		metGrsSe = float(fields[mg_se])

		metGrsWald,metGrsWaldSe,metGrsWaldP = calc_wald(bmiGrsBeta,bmiGrsSe,metGrsBeta,metGrsSe)

		# output results
		outFields = fields + map(str,[bmiSnpWald,bmiSnpWaldSe,bmiSnpWaldP,metGrsWald,metGrsWaldSe,metGrsWaldP])
		outFile.write('\t'.join(outFields) + '\n')

print 'SCRIPT COMPLETED'
