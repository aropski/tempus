#!/usr/bin/env Rscript
# Run this program in the same directory as the run_tempus_func.r and input vcf files

#Rscript master.r variant_output_effect_pick.vcf hg19 
#Rscript master.r <VEP_annotated.vcf> <genome build;e.g. hg19> 

##################
# Initialization #
##################

	# Load required methods
	require('VariantAnnotation')
	require('GenomicRanges')
	require('httr')
	require('jsonlite')
	
	# Load user-defined functions
	source("master_func.r")

	# Convert command line arguments 
	args = commandArgs(trailingOnly=TRUE)
	if (length(args) < 2) {
		stop("Missing arguments!!", call.=FALSE)
	} else if (length(args) == 2) {
		vep_filename <- args[1]	#"variant_output_effect_pick.vcf" 	#annotated vcf from Ensembl VEP
		genome_build <- args[2] #"hg19" 							#genome build
		outfile_name <- 'Challenge_data_final.vcf'
		outtbl_name  <- 'Challenge_data_final.txt'
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
	ALT_CON <- get_field(vcf=vep_vcf, field_name="Consequence", info_type="ANN")
		
		
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

	
###########################################
# Part 4: Bulk query ExAC API; add to vcf #
###########################################

	vep_ExAC_anno <- get_ExAC(vep_vcf)

	# Add allele_freq & rsid from ExAC
	allele_freq <- sapply(vep_ExAC_anno, function(x) get_ExAC_field(x$variant$allele_freq))
	rsid 		<- sapply(vep_ExAC_anno, function(x) get_ExAC_field(x$variant$rsid))
	info(vep_vcf)$ExAC_AF <- as.vector(allele_freq)
	info(vep_vcf)$rsid	  <- as.vector(rsid)
	

############################################################
# Part 5: Write final annotated vcf & create minimal table #
############################################################

	writeVcf(vep_vcf, outfile_name)
	
	# Collect variables of interest for tab-delimited table and write file
	min_vcf_details <- data.frame(
		KEY = rownames(info(vep_vcf)),		 #0. Key identifiers (CHROM POS REF ALT)
		ALT_CON      = ALT_CON, 			 #1. Type of variation 
		ALT_DP       = vep_vcf@info$DP,		 #2. Depth of sequence coverage at variation
		ALT_AO 	     = ALT_AO,				 #3. Number of variant reads  
		ALT_PERC     = ALT_PERC,			 #4. Percentage of variant reads to reference reads
		ExAC_AF      = vep_vcf@info$ExAC_AF, #5. Allele frequency of variant ExAC
		GENE_PHENO   = get_field(vcf=vep_vcf, field_name="GENE_PHENO",  info_type="ANN"), # May link to Cancer Gene Census
		RSID 	     = vep_vcf@info$rsid,
		Gene 	     = get_field(vcf=vep_vcf, field_name="Gene", info_type="ANN"),
		Feature      = get_field(vcf=vep_vcf, field_name="Feature", info_type="ANN"),
		Feature_type = get_field(vcf=vep_vcf, field_name="Feature_type", info_type="ANN")
	)
	write.table(min_vcf_details, outtbl_name, sep="\t", row.names=F)
	
	
	
