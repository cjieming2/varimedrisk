#!/bin/bash

####################################################################
## README
## varimed.sh <flag> <...> >& my.log
## This script tries to build a suite to semi-automate the use of VariMed
## Type in varimed.sh <flag> for more documentation on syntax.
## It assumes the VariMed resource is on a mySQL server
## It requires the following installed and in PATH:
## 1) vcftools
## 2) bcftools
## 3) loadgenotypefile2sql.py (note that the SQL hostname, username and password are hardcoded in here for convenience -- TODO)
## 4) bgzip

####################################################################
### FLAG CODE   ####################################################
### 1) annotate
### >> Use bcftools to annotate ID column of VCF file using VCF.gz file; output: vcf.gz file prefixed with user-defined sample/prefix name
## E.g. varimed.sh annotate dbSNP138 genotypes.vcf.gz dbSNP138.vcf.gz

### 2) upload
### >> Filter only VariMed IDs <rsID file> and upload genotypes to mySQL server; output: vcf.gz file prefixed with 'varimed.snps'. Currently, only SNP analyses, indels are explicitly removed.
### >> NOTE also that currently X, Y and M chrs are included, so if your var_ld tables do not contain these non-autosomes, your results will be wrong (because LD accounting will spew errors).
## E.g. varimed.sh upload P0 genotypes.vcf.gz rubbish <(awk '{OFS="\t"}{FS="\t"}{print "rs"$0}' rsID.txt) P0_genotype proj_dz_risk_mias

### 3) run
### >> Run the risk prediction scripts
## E.g. varimed.sh run host chenj dbname password HG00096_genotype user_chenj.var_varimed_staging_LR_final_ext_sex_eth_spop_EUR var_ld_data_1000genomes_phase3.chr_EUR_spop HG00096_spop_EUR_LR HG00096_spop_EUR_dz_risk HG00096_spop_EUR_snps_all 2 200 0.000001 >&dz_risk_HG00096_1kgp3.EUR.spop.log

### 4) compile
### >> Downloads the results from the mySQL database and then compile
### E.g. varimed.sh compile <(cut -f2,16 samples-whitelist-2565.list | grep -v ADM | sed 1d) proj_dz_risk_ICGC_TCGA >& varimed_compile_ICGC_TCGA.log
### <SAMPLEFILE> = a 2-col file, with col1=sample-id, col2=population; no header


## requires flag code
if [[ -z $1 ]] ; then
	head -32 $0 | sed 1,2d
else
	FLAG=$1
fi

####################################################################
## 1 - use bcftools to annotate ID column of VCF file using VCF.gz file
####################################################################
if [[ ${FLAG} == "annotate" ]] ; then
	
	if [[ "$#" -ne 4 && (-z $2 || -z $3 || -z $4)]]; then
		echo "<Syntax>"
		echo "varimed.sh annotate <PREFIXNAME> <GENOTYPEFILE vcf.gz> <DBSNPFILE vcf.gz>"
		echo ""
		echo "<Example>"
		echo "varimed.sh annotate dbSNP138 genotypes.vcf.gz dbSNP138.vcf.gz"
	
		exit
	fi
	
	## parameters
	NAME=$2
	VCF=$3
	DBSNP=$4
	
	## print to log
	echo "#############################"
	echo "## print parameters to log:"
	echo "#############################"
	echo "FLAG          	      =${FLAG}"
	echo "prefix                =${NAME}"  
	echo "VCF.GZ genotype file  =${VCF}"
	echo "VCF.GZ dbsnp file     =${DBSNP}"
	echo "#############################"
	echo "OUTPUTs:"
	echo "file : ${NAME}_${VCF}"

	## use bcftools to annotate
	bcftools annotate -c ID -a ${DBSNP} ${VCF} | bgzip -c > ${NAME}_${VCF}

fi



####################################################################
## 2 - filter only VariMed IDs and upload genotypes to mySQL server
####################################################################
if [[ ${FLAG} == "upload" ]] ; then
	
	if [[ "$#" -ne 4 && (-z $2 || -z $3 || -z $4)]]; then
		echo "<Syntax>"
		echo "varimed.sh upload <SAMPLENAME> <GENOTYPEFILE vcf.gz> <RSID txt> <DBNAME> <TABLENAME>"
		echo "<RSID txt> = text file with one column of rsID, prepended with 'rs'"
		echo "NOTE also that currently X, Y and M chrs are included in the upload, so if your var_ld tables do not contain these non-autosomes, your results will be wrong (because LD accounting will spew errors)."
		echo ""
		echo "<Example>"
		echo "varimed.sh upload P0 genotypes.vcf.gz rubbish <(awk '{OFS="\t"}{FS="\t"}{print "rs"$0}' rsID.txt) P0_genotype proj_dz_risk_mias"
	
		exit
	fi
	
	
  ## parameters
	NAME=$2
	VCF=$3
	SNP=$4
	DB=$5
	TABLE=$6
	
	
	## print to log
	echo "#############################################"
	echo "## print parameters to log:"
	echo "#############################################"
	echo "FLAG          	        =${FLAG}"
	echo "samplename for vcf2geno =${NAME}"  
	echo "VCF input               =${VCF}"  
	echo "rsID file (to filter)   =${SNP}"  
	echo "upload to database      =${DB}"  
	echo "output table name       =${TABLE}"  
	echo "#############################################"
	echo "OUTPUTs:"
	echo "1) file : varimed.snps.${VCF}"
	echo "2) table: ${TABLE} on db: ${DB}"
	echo "3) tmp files: vcf2geno-${SAMPLE}.txt"
	
	## filter
	echo -e "\n>>> Filtering by VariMed SNPs (rsIDs) using VCFtools..."
	vcftools --gzvcf ${VCF} --out ${VCF} --snps ${SNP} --remove-indels --recode --recode-INFO-all --stdout | gzip -c > varimed.snps.${VCF}
	
	## convert to 3-col file for upload
	bcftools view varimed.snps.${VCF} | bcftools query -f '[%SAMPLE\t%ID\t%TGT]\n' | grep -v "\." | sed -e 's/\///g' -e 's/\|//g' -e 's/rs//g' > vcf2geno-${NAME}.txt
	
	## upload to mySQL
	echo -e "\n>>> Uploading vcf2geno-${NAME}.txt to ${DB}.${TABLE}..."
	loadgenotypefile2sql.py vcf2geno-${NAME}.txt ${TABLE} ${DB}
	
	## delete vcf2geno file
	#rm vcf2geno-${NAME}.txt
	
fi


####################################################################
## 3 - run VariMed codes adapted from Rong Chen's original
####################################################################
if [[ ${FLAG} == "run" ]] ; then
	
	if [[ "$#" -ne 11 && (-z $2 || -z $3 || -z $4 || -z $5 || -z $6 || -z $7 || -z $8 || -z $9 || -z ${10} || -z ${11})]]; then
		echo "<Syntax>"
		echo "varimed.sh run <HOST> <USER> <DB> <password> <INPUT_GENO_TABLE> <INPUT_VARIMED_TABLE> <INPUT_LD_TABLE> <OUTPUT_TABLE_LR> <OUTPUT_TABLE_DZ_RISK> <OUTPUT_TABLE_SNP_ALL> <NUMGWAS> <SIZGWAS> <PMIN>"
		echo ""
		echo "<Example>"
		echo "varimed.sh run buttelab-aws-prod-aurora-cluster.cluster-cd8zgucpvgtu.us-west-2.rds.amazonaws.com chenj proj_dz_risk_1000G_phase3 3VrTh60IlfiHjLATiVkKn8orM HG00096_genotype user_chenj.var_varimed_staging_LR_final_ext_sex_eth_spop_EUR var_ld_data_1000genomes_phase3.chr_EUR_spop HG00096_spop_EUR_LR HG00096_spop_EUR_dz_risk HG00096_spop_EUR_snps_all 2 200 0.000001 >&dz_risk_HG00096_1kgp3.EUR.spop.log"
	
		exit
	fi
	
	
	## parameters
	HOST=$2
	USER=$3
	DB=$4
	PASS=$5
	ITABLEGENOTYPE=$6
	ITABLEVARIMED=$7
	ITABLELD=$8
	OTABLELR=$9
	OTABLEDZRISK=${10}
	OTABLESNPALL=${11}
	NUMGWAS=${12} 
	SIZGWAS=${13}
	PMIN=${14}
	
	## print to log
	echo "#############################"
	echo "## print parameters to log:"
	echo "#############################"
	echo "FLAG          	     =${FLAG}"
	echo "mySQL hostname       =${HOST}"  
	echo "mySQL username       =${USER}"  
	echo "mySQL database       =${DB}"  
	echo "mySQL password       =*"  
	echo "input genotype table =${ITABLEGENOTYPE}"  
	echo "input VariMed table  =${ITABLEVARIMED}"  
	echo "input LD table       =${ITABLELD}"
	echo "output LR table      =${OTABLELR}"         ## LR_R2
	echo "output dz_risk       =${OTABLEDZRISK}"     ## combine_LR_single
	echo "output snps_all      =${OTABLESNPALL}"
	echo "cutoff for num gwas  =${NUMGWAS}"
	echo "cutoff for size gwas =${SIZGWAS}"
	echo "cutoff for min P val =${PMIN}"
	echo "#############################"
	
	## run VariMed
	echo -e "\n>>> Running VariMed/risk scripts..."
	dz_risk_jc.pl -h ${HOST} \
	-u ${USER} \
	-d ${DB} \
	-p ${PASS} \
	-i ${ITABLEGENOTYPE} \
	-j ${ITABLEVARIMED} \
	-k ${ITABLELD} \
	-o ${OTABLELR} \
	-q ${OTABLEDZRISK} \
	-r ${OTABLESNPALL} \
	-g ${NUMGWAS} -s ${SIZGWAS} -t ${PMIN}

	
fi


####################################################################
## 4 - compile results from mySQL server and ready for analyses
####################################################################
if [[ ${FLAG} == "compile" ]] ; then
	
	if [[ "$#" -ne 3 && (-z $2 || -z $3)]]; then
		echo "<Syntax>"
		echo "varimed.sh compile <SAMPLEFILE> <DB>"
		echo "<SAMPLEFILE> = a 2-col file, with col1=sample-id, col2=population; no header"
		echo ""
		echo "<Example>"
		echo "varimed.sh compile <(cut -f2,16 samples-whitelist-2565.list | grep -v ADM | sed 1d) proj_dz_risk_ICGC_TCGA >& varimed_compile_ICGC_TCGA.log"
	
		exit
	fi
	
	## parameters
	SAMPLEFILE=$2
	DB=$3
	
	## print to log
	echo "#############################"
	echo "## print parameters to log:"
	echo "#############################"
	echo "FLAG          	     =${FLAG}"
	echo "sample-file          =${SAMPLEFILE}"  
	echo "mySQL database       =${DB}"
	echo "#############################"
	echo "OUTPUTs:"
	echo "1) combined_dz_risk file"
	echo "2) LR files (in LR folder)"
	echo "3) snps_all files (in snps_all folder)"
	
	## download results from mySQL database
	## combine dz_risk file
	echo -e "\n>>> Downloading and compiling results from mySQL database..."
	( while read name pop; do \
	echo "$name" \ 
	
	## generic
	mysql --database=${DB} -e "select * from "$name"_dz_risk;" | awk '{OFS="\t"}{FS="\t"}{print $0, "'"$name"'", "'"$pop"'"}' | sed 1d >> combined_dz_risk.txt; \
	mysql --database=${DB} -e "select * from "$name"_LR;" > "$name"_LR.txt ; \
	mysql --database=${DB} -e "select * from "$name"_snps_all;" > "$name"_snps_all.txt; \
	done ) < ${SAMPLEFILE} 
	
	sed '1i\broad_phenotype\tLR\tLR_max\tSNP_risk\tSNP_protective\tsample-id\tpopulation' combined_dz_risk.txt | sed 's/ /_/g' > combined_dz_risk.txt_
	mv combined_dz_risk.txt_ combined_dz_risk.txt

	
	## organising
	mkdir LR snps_all
	mv *_LR.txt LR
	mv *_snps_all.txt snps_all
	
	
fi
