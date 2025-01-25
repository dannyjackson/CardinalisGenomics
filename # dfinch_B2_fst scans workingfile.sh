# dfinch_B2_fst scans


while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/subsetsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaurbanbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/noca_urban.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrrurbanbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrr_urban.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaruralbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/noca_rural.txt

while read -r file;
do
grep $file /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrrruralbams.txt
done < /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrr_rural.txt

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

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaurbanbams.txt -out /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/noca_urban -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaruralbams.txt -out /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/noca_rural -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

sbatch angsdsaf_noca.sh 
Submitted batch job 3682520


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

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrrurbanbams.txt -out /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/pyrr_urban -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrrruralbams.txt -out /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/pyrr_rural -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

sbatch angsdsaf_pyrr.sh 
Submitted batch job 3682750

# cat /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaruralbams.txt > /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaall_pyrrurbanbams.txt
# cat /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaurbanbams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaall_pyrrurbanbams.txt
# cat /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/pyrrurbanbams.txt >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/nocaall_pyrrurbanbams.txt



mkdir noca pyrr
mv noca* noca
mv pyrr* pyrr

#this is with 2pops

# NOCA 
#first calculate per pop saf for each population (done above)
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/noca
mkdir depth ml saf
mv *depth* depth/
mv noca_* saf/
#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS saf/noca_urban.saf.idx saf/noca_rural.saf.idx > ml/noca.pre.post.ml

#prepare the fst for easy window analysis etc
mkdir fst_out
~/programs/angsd/misc/realSFS fst index saf/noca_urban.saf.idx saf/noca_rural.saf.idx -sfs ml/noca.pre.post.ml -fstout fst_out/noca

#get the global estimate
~/programs/angsd/misc/realSFS fst stats fst_out/noca.fst.idx 
FST.Unweight[nObs:1144981]:0.026822 Fst.Weight:0.000656

~/programs/angsd/misc/realSFS fst stats2 fst_out/noca.fst.idx  -win 50000 -step 50000 > fst_out/windowed/slidingwindow_noca
~/programs/angsd/misc/realSFS fst stats2 fst_out/noca.fst.idx -win 1 -step 1 >fst_out/snps/singlesnps_noca

# snps 
echo -e 'region\tchr\tmidPos\tNsites\tfst' > fst_out/snps/singlesnps_fst_noca.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' fst_out/snps/singlesnps_noca >> fst_out/snps/singlesnps_fst_noca.txt

sed -i 's/NC_//g' fst_out/snps/singlesnps_fst_noca.txt
sed -i 's/\.1\t/\t/g' fst_out/snps/singlesnps_fst_noca.txt


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

Rscript /xdisk/mcnew/dannyjackson/cardinals_dfinch/fst/snpplot_noca.R

#!/usr/bin/env Rscript

setwd('/xdisk/mcnew/dannyjackson/cardinals_dfinch/fst')
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

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst/pyrr
mkdir depth ml saf
mv *depth* depth/
mv pyrr_* saf/

#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS saf/pyrr_urban.saf.idx saf/pyrr_rural.saf.idx> ml/pyrr.pre.post.ml

#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index pyrr_urban.saf.idx pyrr_rural.saf.idx -sfs pyrr.pre.post.ml -fstout pyrr

#get the global estimate
~/programs/angsd/misc/realSFS fst stats pyrr.fst.idx 
FST.Unweight[nObs:2617605]:0.047155 Fst.Weight:0.030829




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

Rscript /xdisk/mcnew/dannyjackson/cardinals_dfinch/fst/snpplot_pyrr.R

#!/usr/bin/env Rscript

setwd('/xdisk/mcnew/dannyjackson/cardinals_dfinch/fst')
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
# moved all output to /xdisk/mcnew/dannyjackson/cardinals_dfinch/fst/snp_relevantgenes


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

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/fst

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$2".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')' >> relevantgenes_noca_fst.txt

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

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/fst

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$2".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $3}' <<<"${line}"`

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')' >> relevantgenes_pyrr_fst.txt

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

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')'

done < pyrr.outlierfst.csv



grep "VYXE1001099.1" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff 



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


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_noca_windowed_fst.txt

done < noca.outlierfst.windowed.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_windowed_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_windowed_fst.txt

# pyrr
while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print "VYXE0"$1".1"}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $2}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_pyrr_windowed_fst.txt

done < pyrr.outlierfst.windowed.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_windowed_fst.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_windowed_fst.txt