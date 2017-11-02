# make -f Makefile variant_output_effect_pick.vcf Challenge_data_final.vcf

challenge_run: variant_output_effect_pick.vcf Challenge_data_final.vcf

.PHONY: challenge_run


# Annotate most deleterious allele with Ensemble VEP
variant_output_effect_pick.vcf:
	vep -species homo_sapiens --cache --offline --assembly GRCh37 -i Challenge_data.vcf --pick --flag_pick --variant_class --terms SO --allele_number --gene_phenotype --dont_skip --vcf_info_field ANN --fields ALLELE_NUM,CONSEQUENCE,FEATURE,FEATURE_TYPE,GENE_PHENO,GENE,GENE_PHENO,VARIANT_CLASS,HGNC_ID --vcf --output_file variant_output_effect_pick.vcf

# Add read count/freq to vcf
Challenge_data_final.vcf: 
	Rscript master.r variant_output_effect_pick.vcf hg19 

	

