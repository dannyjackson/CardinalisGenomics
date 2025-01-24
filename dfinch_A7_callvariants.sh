#!/bin/bash
module load bcftools/1.19

ref="/xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
bamdir="/xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment/"
ID="cardinalis"

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/genotype_calls

bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


sbatch --account=mcnew \
--job-name=callvariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.callvariants.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=240:00:00 \
callvariants.sh

# Submitted batch job 12032184