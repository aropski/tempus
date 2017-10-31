#make -f Makefile variant_output_effect_pick.vcf Challenge_data_ExAC.vcf Challenge_data_final.vcf


challenge_run: variant_output_effect_pick.vcf Challenge_data_ExAC.vcf Challenge_data_final.vcf

.PHONY: challenge_run


# Annotate most deleterious allele with Ensemble VEP
variant_output_effect_pick.vcf:
	vep -species homo_sapiens --cache --offline --assembly GRCh37 -i Challenge_data.vcf --pick --flag_pick --pick_order rank,canonical,appris,tsl,biotype,ccds,length --variant_class --terms SO --allele_number --gene_phenotype --dont_skip --vcf_info_field ANN --fields ALLELE_NUM,GENE_PHENO,VARIANT_CLASS,HGNC_ID --vcf --output_file variant_effect_output_pick.vcf

# Annotate ExAC population allele frequency with Ensemble VEP
Challenge_data_ExAC.vcf: 
	vep -species homo_sapiens --cache --offline --assembly GRCh37 -i variant_effect_output_pick.vcf --pick --dont_skip --plugin ExAC,/home/allan/.vep/ExAC.0.3.GRCh37.vcf.gz,AF --vcf --output_file ExAC.vcf

# Add read count/freq to vcf
Challenge_data_final.vcf: 
	Rscript master.r ExAC.vcf hg19 

