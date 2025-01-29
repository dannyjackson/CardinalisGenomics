# B2_fst.sh

FST scans 
bayescan
nucleotide diversity θπ 
Tajimas D 
SweeD 
OmegaPlus

FST
dxy
RAiSD
SweeD
Tajima’s D
# OmegaPlis
# bayescan

awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info.txt
# noca urban
echo -e "NOCA003\nNOCA004\nNOCA006\nNOCA008\nNOCA012\nNOCA013" > /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_urban.txt

# pyrr urban
echo -e "PYRR003\nPYRR004\nPYRR006\nPYRR007\nPYRR009\nPYRR011" > /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_urban.txt

# pyrr rural
echo -e "MSB25201\nUWBM103346\nUWBM77548\nUWBM77718\nUWBM77780\nUWBM77781" > /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_rural.txt

# nocal rural
echo -e "UWBM100619\nUWBM100620\nUWBM100621\nUWBM103345\nUWBM77856\nUWBM77978" > /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_rural.txt

awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info.txt | tail -n +2 > /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_names.txt



#!/bin/bash

#SBATCH --job-name=angsdsaf_pops
#SBATCH --ntasks=10
#SBATCH --nodes=10            
#SBATCH --time=25:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_pops.%j

module load bcftools/1.19

bcftools reheader -s /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_names.txt -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.rehead.bcf /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.bcf

bcftools view -Ou -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/noca_urban.bcf -S /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_urban.txt /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.rehead.bcf --threads 10

bcftools view -Ou -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/pyrr_urban.bcf -S /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_urban.txt /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.rehead.bcf --threads 10

bcftools view -Ou -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/pyrr_rural.bcf -S /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_rural.txt /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.rehead.bcf --threads 10

bcftools view -Ou -o /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/subsets/noca_rural.bcf -S /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_rural.txt /xdisk/mcnew/dannyjackson/cardinals/vcf_likelihoods/all/genolike.rehead.bcf --threads 10

sbatch angsdsaf.sh 
Submitted batch job 3353868

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaurbanbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_urban.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrrurbanbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_urban.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaruralbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/noca_rural.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrrruralbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrr_rural.txt

#!/bin/bash

#SBATCH --job-name=angsdsaf
#SBATCH --ntasks=10
#SBATCH --nodes=1          
#SBATCH --time=50:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf.%j

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaurbanbams.txt -out /xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals/referencelists/sites_headless.mafs -anc /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -nThreads 10

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaruralbams.txt -out /xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals/referencelists/sites_headless.mafs -anc /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -nThreads 10

sbatch angsdsaf.sh 
Submitted batch job 11101345

#!/bin/bash

#SBATCH --job-name=angsdsaf
#SBATCH --ntasks=10
#SBATCH --nodes=1          
#SBATCH --time=50:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf.%j

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrrurbanbams.txt -out /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_urban -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals/referencelists/sites_headless.mafs -anc /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -nThreads 10

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals/referencelists/pyrrruralbams.txt -out /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_rural -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals/referencelists/sites_headless.mafs -anc /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -nThreads 10

sbatch angsdsaf2_pyrr.sh 
Submitted batch job 3407533

cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaruralbams.txt > /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaall_pyrrurbanbams.txt
cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaurbanbams.txt >> /xdisk/mcnew/dannyjackson/cardinals/referencelists/nocaall_pyrrurbanbams.txt



#this is with 2pops
#first calculate per pop saf for each population (done above)
#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS noca_urban.saf.idx noca_rural.saf.idx > noca.pre.post.ml
~/programs/angsd/misc/realSFS pyrr_urban.saf.idx pyrr_rural.saf.idx> pyrr.pre.post.ml



#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index noca_urban.saf.idx noca_rural.saf.idx -sfs noca.pre.post.ml -fstout noca
~/programs/angsd/misc/realSFS fst index pyrr_urban.saf.idx pyrr_rural.saf.idx -sfs pyrr.pre.post.ml -fstout pyrr




#get the global estimate
~/programs/angsd/misc/realSFS fst stats noca.fst.idx 
FST.Unweight[nObs:2617562]:0.044727 Fst.Weight:0.028677

~/programs/angsd/misc/realSFS fst stats pyrr.fst.idx 
FST.Unweight[nObs:2617605]:0.047155 Fst.Weight:0.030829

# NOCA
~/programs/angsd/misc/realSFS fst stats2 noca.fst.idx  -win 50000 -step 50000 >slidingwindow_noca
~/programs/angsd/misc/realSFS fst stats2 noca.fst.idx -win 1 -step 1 >slidingwindow_singlesnps_noca

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst_noca.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'VYXE' slidingwindow_singlesnps_noca >> slidingwindow_singlesnps_fst_noca.txt

sed -i 's/VYXE//g' slidingwindow_singlesnps_fst_noca.txt 
sed -i 's/\.1\t/\t/g' slidingwindow_singlesnps_fst_noca.txt


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst_noca.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'VYXE' slidingwindow_noca >> slidingwindow_fst_noca.txt


sed -i 's/VYXE//g' slidingwindow_fst_noca.txt 
sed -i 's/\.1\t/\t/g' slidingwindow_fst_noca.txt

# First, create a plot of windowed output
module load R/4.4.0
R
library(qqman)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

fst <- read.csv('slidingwindow_fst_noca.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA
max(fst$fst)
# 0.629884

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)


ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nrow(ordered_fst)
[1] 86349
# 86349 snps, so top 0.1% would be 86

outlier_fst_disorder <- ordered_fst[1:86,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.386049
max(outlier_fst_disorder$fst)
# 0.629884

max(fst$neg_log_pvalues_one_tailed)
mean(fst$neg_log_pvalues_one_tailed)
min(fst$neg_log_pvalues_one_tailed)

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, "noca.outlierzfst_window.csv")


pdf(file = "manhattan_windowed_noca_logp.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="neg_log_pvalues_one_tailed", logp=FALSE, ylab = "-log(p-value)", cex = 0.5)
dev.off()

pdf(file = "manhattan_windowed_noca_fst.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()



# second, create list of outlier snps

fst <- read.csv('slidingwindow_singlesnps_fst_noca.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)



ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nrow(ordered_fst)
[1] 2591546
# 2591546 snps, so top 0.1% would be 2591


outlier_fst_disorder <- ordered_fst[1:2592,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.462061
max(outlier_fst_disorder$fst)
# 0.849009

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, "noca.outlierfst.csv")






# create a list of outlier windows

# noca
fst <- read.csv('slidingwindow_fst_noca.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)



ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nrow(ordered_fst)
[1] 86349
# 86349 snps, so top 0.1% would be 86


outlier_fst_disorder <- ordered_fst[1:86,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.386049
max(outlier_fst_disorder$fst)
# 0.629884

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, "noca.outlierfst.windowed.csv")

# pyrr

fst <- read.csv('slidingwindow_fst_pyrr.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)



ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nrow(ordered_fst)
[1] 14612
# 86349 snps, so top 0.1% would be 86


outlier_fst_disorder <- ordered_fst[1:86,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.217238
max(outlier_fst_disorder$fst)
# 0.51291

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, "pyrr.outlierfst.windowed.csv")

# finally, plot snps


#!/bin/bash

#SBATCH --job-name=snpplot
#SBATCH --ntasks=10
#SBATCH --nodes=1          
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.snpplot.%j

module load R/2.2.0

Rscript /xdisk/mcnew/dannyjackson/cardinals/fst/snpplot_noca.R

#!/usr/bin/env Rscript

setwd('/xdisk/mcnew/dannyjackson/cardinals/fst')
library(qqman)

fst <- read.csv('slidingwindow_singlesnps_fst_noca.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA

# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)

# png(file = "manhattan_snps_noca.png", width = 1425, height =975)
# manhattan(fst, chr="chr", bp="midPos", snp="Nsites", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
# dev.off()


pdf(file = "manhattan_snps_noca_logp.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="neg_log_pvalues_one_tailed", logp=FALSE, ylab = "-log(p-value)", cex = 0.5)
dev.off()

sbatch snpplot_noca.sh 
Submitted batch job 3408647
# with logp
3411013







# PYRR


~/programs/angsd/misc/realSFS fst stats2 pyrr.fst.idx  -win 50000 -step 50000 >slidingwindow_pyrr
~/programs/angsd/misc/realSFS fst stats2 pyrr.fst.idx -win 1 -step 1 >slidingwindow_singlesnps_pyrr

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst_pyrr.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'VYXE' slidingwindow_singlesnps_pyrr >> slidingwindow_singlesnps_fst_pyrr.txt

sed -i 's/VYXE//g' slidingwindow_singlesnps_fst_pyrr.txt 
sed -i 's/\.1\t/\t/g' slidingwindow_singlesnps_fst_pyrr.txt


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst_pyrr.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'VYXE' slidingwindow_pyrr >> slidingwindow_fst_pyrr.txt


sed -i 's/VYXE//g' slidingwindow_fst_pyrr.txt 
sed -i 's/\.1\t/\t/g' slidingwindow_fst_pyrr.txt



# first plot the windowed output
module load R/4.4.0
R
library(qqman)
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

fst <- read.csv('slidingwindow_singlesnps_fst_pyrr.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA
max(fst$fst)
0.988323


# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)



ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(neg_log_pvalues_one_tailed)) 

nrow(ordered_fst)
[1] 2591590
# 2591590 snps, so top 0.1% would be 2591

outlier_fst_disorder <- ordered_fst[1:2592,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.471459
max(outlier_fst_disorder$fst)
# 0.988323

outlier_fst_disorder2 <- subset(outlier_fst_disorder, select = -c(region))

write.csv(outlier_fst_disorder2, "pyrr.outlierfst.csv")



fst <- read.csv('slidingwindow_fst_pyrr.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA


# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)




pdf(file = "manhattan_windowed_pyrr_logp.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="neg_log_pvalues_one_tailed", logp=FALSE, ylab = "-log(p-value)", cex = 0.5)
dev.off()

pdf(file = "manhattan_windowed_pyrr.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()


#!/bin/bash

#SBATCH --job-name=snpplot
#SBATCH --ntasks=10
#SBATCH --nodes=1          
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.snpplot.%j

module load R/2.2.0

Rscript /xdisk/mcnew/dannyjackson/cardinals/fst/snpplot_pyrr.R

#!/usr/bin/env Rscript

setwd('/xdisk/mcnew/dannyjackson/cardinals/fst')
library(qqman)

fst <- read.csv('slidingwindow_singlesnps_fst_pyrr.txt', sep ='\t')
fst_noNA <- na.omit(fst)
nrow(fst) - nrow(fst_noNA)
fst = fst_noNA


# z transform fst values

fst_xbar = mean(fst$fst)
fst_sd = sd(fst$fst)

fst$z <- (fst$fst - fst_xbar)/fst_sd

p_values_one_tailed <- pnorm(q=fst$z, lower.tail=FALSE)

# Calculate -log10(p-value)
fst$neg_log_pvalues_one_tailed <- -log10(p_values_one_tailed)


pdf(file = "manhattan_snps_pyrr_logp.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="neg_log_pvalues_one_tailed", logp=FALSE, ylab = "-log(p-value)", cex = 0.5)
dev.off()


png(file = "manhattan_snps_pyrr.png", width = 1425, height =975)
manhattan(fst, chr="chr", bp="midPos", snp="Nsites", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()

sbatch snpplot_pyrr.sh 
Submitted batch job 3411043




# analyze significant genes in snps
# moved all output to /xdisk/mcnew/dannyjackson/cardinals/fst/snp_relevantgenes


# noca

#!/bin/bash

#SBATCH --job-name=relevantgenesnoca
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.relevantgenesnoca.%j

cd /xdisk/mcnew/dannyjackson/cardinals/fst

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$2".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')' >> relevantgenes_noca_fst.txt

done < noca.outlierfst.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_fst.txt

sbatch relevantgenes_noca.sh 
sSubmitted batch job 3410502

# pyrr

#!/bin/bash

#SBATCH --job-name=relevantgenespyrr
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.relevantgenespyrr.%j

cd /xdisk/mcnew/dannyjackson/cardinals/fst

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$2".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')' >> relevantgenes_pyrr_fst.txt

done < pyrr.outlierfst.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_fst.txt


sbatch relevantgenes_pyrr.sh 
Submitted batch job 3410501



while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$2".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`
echo $chr
echo $midpos 

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')'

done < pyrr.outlierfst.csv



grep "VYXE1001099.1" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff 



# windowed relevant genes


awk 'BEGIN {FS = ","} {$1=""}1' noca.outlierfst.windowed.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > noca.outlierfst.windowed.tmp.csv

mv noca.outlierfst.windowed.tmp.csv noca.outlierfst.windowed.csv

awk 'BEGIN {FS = ","} {$1=""}1' pyrr.outlierfst.windowed.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > pyrr.outlierfst.windowed.tmp.csv

mv pyrr.outlierfst.windowed.tmp.csv pyrr.outlierfst.windowed.csv

sed -i 's/\"//g' noca.outlierfst.windowed.csv
sed -i 's/\"//g' pyrr.outlierfst.windowed.csv

# noca

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $2}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_noca_windowed_fst.txt

done < noca.outlierfst.windowed.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_windowed_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_windowed_fst.txt

# pyrr
while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $2}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_pyrr_windowed_fst.txt

done < pyrr.outlierfst.windowed.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_windowed_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_windowed_fst.txt