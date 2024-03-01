#!/bin/bash


### Germline variant calling pipeline
### GENE 5150
### Maintainer Paulo Joshua Tanicala
### ptanicala@kgi.edu

#####NOTE:
### son sequencing data was used for this pipeline. mother and father data were processed by other groups.
### file paths needs to be modified for mother and father sequencing data.


##### Step 0. Download and create index file for alignment.

wget --no-check-certificate  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa
gatk-4.3.0.0/gatk CreateSequenceDictionary R=hg38.fa O=hg38.dict

##### Step 1. Use bwa mem to align reads against reference genome.
bwa mem -t 32 -R "@RG\tID:SON\tPL:ILLUMINA\tSM:SON" ../cesar/TEST/hg38.fa \
	 ../cesar/DATASETS/SON/SON_R1.fastq.gz ../cesar/DATASETS/SON/SON_R2.fastq.gz > son.paired.sam

#### Step 2. Sort reads using samtools.
samtools sort -o son.sorted.paired.sam son.paired.sam

##### Step 3. Remove duplicates
../cesar/TEST/gatk-4.3.0.0/gatk MarkDuplicates -I son.sorted.paired.sam \
	-O son_sorted_dedup_reads.bam -M metrics

##### Step 4. Recalibrate. Maybe lab error

../cesar/TEST/gatk-4.3.0.0/gatk BaseRecalibrator -I son_sorted_dedup_reads.bam \
	-R ../cesar/TEST/hg38.fa --known-sites ../cesar/TEST/Homo_sapiens_assembly38.dbsnp138.vcf \
	-O recal_data.table


##### Step 5. Use recalibration data to adjust quality scores
../cesar/TEST/gatk-4.3.0.0/gatk ApplyBQSR -I son_sorted_dedup_reads.bam \
	-R ../cesar/TEST/hg38.fa --bqsr-recal-file recal_data.table \
	-O son_sorted_bqsr_dedup_reads.bam

##### Step 6. Collect alignment metrics 

../cesar/TEST/gatk-4.3.0.0/gatk CollectAlignmentSummaryMetrics R=hg38.fa \
	I=son_sorted_bqsr_dedup_reads.bam O=alignment_metrics.txt

../cesar/TEST//gatk-4.3.0.0/gatk CollectInsertSizeMetrics \
INPUT=son_sorted_bqsr_dedup_reads.bam \
OUTPUT=insert_size_metrics.txt HISTOGRAM_FILE=histogram.pdf

##### Step 7. Call variants

../cesar/TEST/gatk-4.3.0.0/gatk HaplotypeCaller -R hg38.fa \
	-I son_sorted_bqsr_dedup_reads.bam -O raw_variants.vcf

##### Step 8. Split variants into SNPS and INDELS

../cesar/TEST/gatk-4.3.0.0/gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type SNP -O raw_snps.vcf
../cesar/TEST/gatk-4.3.0.0/gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type INDEL -O raw_indels.vcf


#######################2/1/24 CLASS ENDED HERE
#######################2/8/24 & 2/14/24 CLASS CONTINUED HERE
##### Step 9. Filters both INDEL and SNPs

####SNPs filter
../cesar/TEST/gatk-4.3.0.0/gatk VariantFiltration \
-R ../cesar/TEST/hg38.fa \
-V raw_snps.vcf \
-O filtered_snps.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 10" \
-genotype-filter-name "GQ_filter"


####INDEL filter
../cesar/TEST/gatk-4.3.0.0/gatk VariantFiltration \
-R ../cesar/TEST/hg38.fa \
-V raw_indels.vcf \
-O filtered_indels.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 200.0" \
-filter-name "SOR_filter" -filter "SOR > 10.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 10" \
-genotype-filter-name "GQ_filter"

##### Step 10. After the variants are filtered, they will have "PASS". These need to be selected

####SNPs selection
../cesar/TEST/gatk-4.3.0.0/gatk SelectVariants \
--exclude-filtered \
-V filtered_snps.vcf \
-O analysis-ready-snps.vcf

####INDEL selection
../cesar/TEST/gatk-4.3.0.0/gatk SelectVariants \
--exclude-filtered \
-V filtered_indels.vcf \
-O analysis-ready-indels.vcf


##### Step 11. Filtering for Genotype

#finds all "DP_filter" and ignores them
grep -v “DP_filter” analysis-ready-snps.vcf
#output into new file (high quality)
grep -v “DP_filter” analysis-ready-snps.vcf > son_final_set_snps.vcf

DP = quality over genotype

cat analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf

##### Step 12. Functional annotation using GATK Funcotator

#File input was output from grep "DP_filter" step.
#For SNPs
../cesar/TEST/gatk-4.3.0.0/gatk Funcotator --variant son_final_set_snps.vcf --reference ../cesar/TEST/hg38.fa --ref-version hg38 --data-sources-path ../cesar/TEST/funcotator_dataSources.v1.7.20200521g --output analysis-ready-snps-filteredGT-functotated.vcf --output-file-format VCF

#For INDELs
../cesar/TEST/gatk-4.3.0.0/gatk Funcotator --variant analysis-ready-indels-filteredGT.vcf --reference hg38.fa --ref-version hg38 --data-sources-path funcotator_dataSources.v1.7.20200521g --output analysis-ready-indels-filteredGT-functotated.vcf --output-file-format VCF

##############Output files from this pipline (.vcf) was manipulated using Google Colab