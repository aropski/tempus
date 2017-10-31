#!/usr/bin/env Rscript
# Run this program in the same directory as the run_tempus_func.r and input vcf files

#Rscript master.r ExAC.vcf hg19 
#Rscript master.r <VEP_annotated.vcf> <genome build;e.g. hg19> 




##################
# Initialization #
##################

	# Load required methods
	require('VariantAnnotation')
	require('GenomicRanges')
	require('vcfR') 

	# Load user-defined functions
	source("run_tempus_func.r")

	# Convert command line arguments 
	args = commandArgs(trailingOnly=TRUE)
	if (length(args) < 2) {
		stop("Missing arguments!!", call.=FALSE)
	} else if (length(args) == 2) {
		vep_filename <- args[1]	#"ExAC.vcf" 				#annotated vcf
		genome_build <- args[2] #"hg19" 					#genome build
		outfile_name <- 'final_Challenge_data.vcf'
	} else {
		stop("Too many arguments!!!!", call.=FALSE)
	}


###########################################################
# Part 1: Read in the data w/ VEP Annotations (ANN & CSQ) #
###########################################################

	my_vep <- readVcf(vep_filename, genome_build)

		
################################################
#  Part 2: Extract which variant was annotated #
################################################

	my_ALLELE_NUM <- get_field(vcf=my_vep, field_name="ALLELE_NUM",  info_type="ANN")
	my_ALLELE_CON <- get_field(vcf=my_vep, field_name="Consequence", info_type="CSQ")
		
		
######################################################################
# Part 3: Compute coverage / read count / %ALT-REF for picked allele #
######################################################################

	my_ALLELE_AO   <- get_ALT_INFO_value(vcf=my_vep, field_name="AO", my_ALLELE_NUM)
	my_ALLELE_RO   <- get_ALT_INFO_value(vcf=my_vep, field_name="RO", my_ALLELE_NUM)
	my_ALLELE_PERC <- my_ALLELE_AO / my_ALLELE_RO

	# Add to vcf file
	 info(my_vep)$my_ALLELE_DP   <- my_vep@info$DP
	 info(my_vep)$my_ALLELE_AO   <- my_ALLELE_AO
	 info(my_vep)$my_ALLELE_PERC <- my_ALLELE_PERC
	 
	
############################################################
# Part 4: Write final annotated vcf & create jpg summaries #
############################################################

	writeVcf(my_vep, outfile_name)

	
	#vcf_final <- read.vcfR("final_Challenge_data.vcf", verbose=FALSE)
	#chrom < create.chromR(name="CHROM", vcf=vcf_final)
	#chrom < proc.chromR(chrom, verbose=FALSE)
	
	#jpeg('smry_2.jpg')
	#plot(chrom)
	#dev.off()

	#jpeg('smry_3.jpg')
	#chromoqc(chrom, dp.alpha = 22)
	#dev.off()


	

	
