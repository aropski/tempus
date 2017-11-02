

# Function for extracting allele-specific details (as picked by Ensembl VEP 'ALLELE_NUM')
get_ALT_INFO_value <- function(vcf, field_name, ALLELE_NUM) {

	# Extract all ALT allele counts for field_name
	field_values <- vcf@info[field_name]
	
	# Extract value for picked allele; indexed by 'ALLELE_NUM'
	my_ALT <- NULL
	for (idx in seq(nrow(field_values))) {
		ALLELE_idx <- as.numeric(ALLELE_NUM[idx])
		my_ALT <- c(my_ALT, as.numeric(field_values[[1]][[idx]][ALLELE_idx]))
	}
	return(my_ALT)
}


# Function for extracting annotation descriptions, matching by field_name to INFO header
get_field <- function(vcf, field_name, info_type) {

	# Determine which info row contains info_type of interest (e.g. ANN or CSQ)
	info_idx <- which(rownames(info(header(vcf))) == info_type)
	# Extract Description
	info_Desc <- info(header(vcf))[info_idx,'Description']
	# Drop pre-text and convert to vector
	info_fields <- unlist(strsplit(strsplit(info_Desc, "Format: ")[[1]][2], "|", fixed=T))
	# Determine which Description matches field_name argument
	field_name_idx <- which( info_fields == field_name )
	# Extract field and return
	return(sapply( vcf@info[[info_type]], function(x) strsplit(x ,"|", fixed=T)[[1]][field_name_idx]  ))
	
}


# Convert vcf format to ExAC Variant Id as: CHROMOSOME-POSITION-REFERENCE-VARIANT
convert_vcf_to_exac <-  function(vcf_key) {
	return(gsub(":|/|_", "-", vcf_key))
}

	
# Function for sending bulk request to ExAC API
get_ExAC <- function(vep_vcf) {

	# Convert vcf format to ExAC Variant Id as: CHROMOSOME-POSITION-REFERENCE-VARIANT
	JSON_KEY <- as.vector(sapply(rownames(info(vep_vcf)), function (x) convert_vcf_to_exac(x)))

	# Query in bulk; does not return same order
	req <- POST("http://exac.hms.harvard.edu/rest/bulk/variant", 
				body = toJSON(JSON_KEY), encode='json')
	stop_for_status(req)
	json <- fromJSON(content(req, "text"))

	# Match returned json to vep_vcf order
	json_adj <- json[match(JSON_KEY, names(json))]

	# Audit 
	if(!all(names(json_adj) == JSON_KEY)) {
		stop("ExAC mapping does not match vep_vcf names!!!")
	} 
	
	return(json_adj)
	
}


# Function for preventing NULL values for missing ExAC annotations
get_ExAC_field <- function(x) {
	if(is.null(x)) {
		return(NA)
	} else {
		return(x[[1]])
	}
}
	

# calculated variant consequences rank of SO terms used in --most_severe; 
# this is primary pick_order used in VEP, as shown in: 
# https://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
SO_term_rank <- c('Transcript ablation', 'Splice acceptor variant', 'Splice donor variant', 'Stop gained', 
				'Frameshift variant', 'Stop lost', 'Start lost', 'Transcript amplification', 'Inframe insertion', 
				'Inframe deletion', 'Missense variant', 'Protein altering variant', 'Splice region variant', 
				'Incomplete terminal codon variant', 'Stop retained variant', 'Synonymous variant', 'Coding sequence variant', 
				'Mature miRNA variant', '5 prime UTR variant', '3 prime UTR variant', 'Non coding transcript exon variant', 
				'Intron variant', 'NMD transcript variant', 'Non coding transcript variant', 'Upstream gene variant', 
				'Downstream gene variant', 'TFBS ablation', 'TFBS amplification', 'TF binding site variant', 
				'Regulatory region ablation', 'Regulatory region amplification', 'Feature elongation', 
				'Regulatory region variant', 'Feature truncation', 'Intergenic variant')


