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
	source("master_func.r")

	# Convert command line arguments 
	args = commandArgs(trailingOnly=TRUE)
	if (length(args) < 2) {
		stop("Missing arguments!!", call.=FALSE)
	} else if (length(args) == 2) {
		vep_filename <- args[1]	#"ExAC.vcf" 	#annotated vcf from 2 Ensembl runs
		genome_build <- args[2] #"hg19" 		#genome build
		outfile_name <- 'final_Challenge_data.vcf'
		outtbl_name  <- 'final_Challenge_data.txt'
	} else {
		stop("Too many arguments!!!!", call.=FALSE)
	}


###########################################################
# Part 1: Read in the data w/ VEP Annotations (ANN & CSQ) #
###########################################################

	vep_vcf <- readVcf(vep_filename, genome_build)

		
################################################
#  Part 2: Extract which variant was annotated #
################################################

	ALT_NUM <- get_field(vcf=vep_vcf, field_name="ALLELE_NUM",  info_type="ANN")
	ALT_CON <- get_field(vcf=vep_vcf, field_name="Consequence", info_type="CSQ")
		
		
######################################################################
# Part 3: Compute coverage / read count / %ALT-REF for picked allele #
######################################################################

	ALT_AO   <- get_ALT_INFO_value(vcf=vep_vcf, field_name="AO", ALT_NUM)
	ALT_RO   <- get_ALT_INFO_value(vcf=vep_vcf, field_name="RO", ALT_NUM)
	ALT_PERC <- ALT_AO / ALT_RO

	# Add to vcf object 
	info(vep_vcf)$ALT_DP   <- vep_vcf@info$DP
	info(vep_vcf)$ALT_AO   <- ALT_AO
	info(vep_vcf)$ALT_PERC <- ALT_PERC
	 
	 	
############################################################
# Part 4: Write final annotated vcf & create minimal table #
############################################################

	writeVcf(vep_vcf, outfile_name)
	
	# Collect variables of interest for tab-delimited table and write file
	min_vcf_details <- data.frame(
		KEY = rownames(info(vep_vcf)),	#0. Key identifiers (CHROM POS REF ALT)
		ALT_CON  = ALT_CON, 			#1. Type of variation 
		ALT_DP   = vep_vcf@info$DP,		#2. Depth of sequence coverage at variation
		ALT_AO 	 = ALT_AO,				#3. Number of variant reads  
		ALT_PERC = ALT_PERC,			#4. Percentage of variant reads to reference reads
		ExAC_AF  = get_field(vcf=vep_vcf, field_name="ExAC_AF", info_type="CSQ"), #5. Allele frequency of variant ExAC
		Gene 	 = get_field(vcf=vep_vcf, field_name="Gene", info_type="CSQ"),
		Feature  = get_field(vcf=vep_vcf, field_name="Feature", info_type="CSQ"),
		Feature_type = get_field(vcf=vep_vcf, field_name="Feature_type", info_type="CSQ")
	)
	write.table(min_vcf_details, outtbl_name, sep="\t", row.names=F)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



	

	
