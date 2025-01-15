# A6_snpID

cd /xdisk/mcnew/dannyjackson/cardinals/referencelists

ls /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/*bam > allsamplebams.txt


#!/bin/bash

#SBATCH --job-name=snps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.snps.%j

cd /xdisk/mcnew/dannyjackson/cardinals/referencelists

~/programs/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt -out /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsnps -nThreads 10 

sbatch snps.sh 
Submitted batch job 10711013

zcat allsnps.mafs.gz | awk '{print $1, $2, $3, $4}' > sites.mafs

tail -n +2 sites.mafs > sites_headless.mafs

~/programs/angsd/angsd sites index sites_headless.mafs



#!/bin/bash

#SBATCH --job-name=vcf_likelihoods
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.vcf_likelihoods.%j

cd /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/cardinals/referencelists/sites_headless.mafs -doBcf 1 -doGlf 2 -nThreads 12 -out genolike

sbatch vcf_all.sh 
Submitted batch job 3353611

