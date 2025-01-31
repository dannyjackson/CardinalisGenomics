# B1_PCA.sh no trans, all individuals

# for pcangsd, I need to subset by unlinked snps
# and to prune linked SNPs using a window size of 50kb, a step size of 5 SNPs, and an r2 threshold of 0.5
~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsamplebams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs -doBcf 1 -doGlf 2 -nThreads 12 -out genolike_notrans -noTrans 1

#!/bin/bash

module load python/3.11/3.11.4
module load R
module load plink/1.9
module load samtools

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans

~/programs/angsd/angsd -vcf-gl /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_notrans.bcf -doPlink 2 -out genolike_plink -doGeno -1 -dopost 1 -domaf 1 -doMajorMinor 1 -sites /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sites_headless.mafs


plink --tped /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans/genolike_plink.tped --tfam /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered

plink --tped /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans/genolike_plink.tped --tfam /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 

# bgzip genolike_pruned.beagle.bed

sbatch --account=mcnew \
--job-name=pruneforlinkage \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/pruneforlinkage%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=24:00:00 \
pruneforlinkage.sh

# Submitted batch job 12043557

#!/bin/bash
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/all_notrans
# for selection / this is only able to look at differences between species not urban/rural
module load python

pcangsd -b /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_notrans.beagle.gz -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/all_notrans/all_notrans -t 12 -e 2 --selection --pcadapt --sites_save --snp_weights

# for pop gen
pcangsd -p /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all_notrans/genolike_pruned -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/all_notrans/pruned -t 12 -e 2 --selection --admix


sbatch --account=mcnew \
--job-name=selection_pca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/selection_pca%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=24:00:00 \
selection_pca_subset.sh


# it's possible that UWBM77548 and UWBM77718 are related...

# save the following as: /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info.txt

Sample  Species Treatment  
MSB25201    PYRR    Rural
NOCA003 NOCA    Urban
NOCA004 NOCA    Urban
NOCA006 NOCA    Urban
NOCA008 NOCA    Urban
NOCA012 NOCA    Urban
NOCA013 NOCA    Urban
PYRR003 PYRR    Urban
PYRR004 PYRR    Urban
PYRR006 PYRR    Urban
PYRR007 PYRR    Urban
PYRR009 PYRR    Urban
PYRR011 PYRR    Urban
UWBM100619  NOCA    Rural
UWBM100620  NOCA    Rural
UWBM100621  NOCA    Rural
UWBM103345  NOCA    Rural
UWBM103346  PYRR    Rural
UWBM77548   PYRR    Rural
UWBM77718   PYRR    Rural
UWBM77780   PYRR    Rural
UWBM77781   PYRR    Rural
UWBM77856   NOCA    Rural
UWBM77978   NOCA    Rural

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/all_notrans/

C <- as.matrix(read.table("pruned.cov")) # Reads estimated covariance matrix
D <- as.matrix(read.table("output.selection")) # Reads PC based selection statistics
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# Obtain p-values from PC-based selection scan
p <- pchisq(D, 1, lower.tail=FALSE)
write.table(p, "pvalues.txt")

# plot admixture
tbl=read.table("pruned.admix.3.Q")
pdf(file = "admix_pruned.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()





# can we drop one of UWBM77718, UWBM77548?
# drop 718 from this and all future analyses. It has a slightly lower average depth and a higher standard deviation.

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/
echo 'UWBM77718 1' > subsetindv.txt
awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info.txt | tail -n +2 > /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/samplenames.txt

awk 'BEGIN{FS=OFS=" "} NR==FNR{a[NR]=$0; next} {$1=a[FNR]} 1' /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/samplenames.txt genolike_pruned.fam > genolike_pruned.renamed.fam

mv genolike_pruned.renamed.fam genolike_pruned.fam

plink --bed /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.bim --fam /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.fam --allow-extra-chr --snps-only 'just-acgt' --remove subsetindv.txt --out genolike_pruned_subset --make-bed 


pcangsd -p /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned_subset -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/pruned_subset -t 12 -e 2 --selection --admix


Sample  Species Treatment  
MSB25201    PYRR    Rural
NOCA003 NOCA    Urban
NOCA004 NOCA    Urban
NOCA006 NOCA    Urban
NOCA008 NOCA    Urban
NOCA012 NOCA    Urban
NOCA013 NOCA    Urban
PYRR003 PYRR    Urban
PYRR004 PYRR    Urban
PYRR006 PYRR    Urban
PYRR007 PYRR    Urban
PYRR009 PYRR    Urban
PYRR011 PYRR    Urban
UWBM100619  NOCA    Rural
UWBM100620  NOCA    Rural
UWBM100621  NOCA    Rural
UWBM103345  NOCA    Rural
UWBM103346  PYRR    Rural
UWBM77548   PYRR    Rural
UWBM77780   PYRR    Rural
UWBM77781   PYRR    Rural
UWBM77856   NOCA    Rural
UWBM77978   NOCA    Rural

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/

C <- as.matrix(read.table("pruned_subset.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info_subset.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_subset_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_subset_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_subset_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("pruned_subset.admix.3.Q")
pdf(file = "admix_pruned_subset.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()



# rerun selection using subset individuals


sbatch --account=mcnew \
--job-name=selection_pca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.selection_pca.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=2:00:00 \
selection_pca_subset.sh

#!/bin/bash

# for selection / this is only able to look at differences between species not urban/rural
module load python

pcangsd -b /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike.beagle.gz -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/subset -t 12 -e 2 --selection --pcadapt --sites_save --snp_weights

Submitted batch job 12043363


# https://github.com/Rosemeis/pcangsd/blob/master/scripts/pcadapt.R

####################
# Run pcadapt scan #
####################
#
#	args[1]: zscores matrix
#	args[2]: prefix for output files
#
module load R

R

args = commandArgs(trailingOnly=TRUE)
library(bigutilsr)

zscores <- read.table('all.pcadapt.zscores')
K <- ncol(zscores)
z_mat <- as.matrix(zscores)

d2 <- dist_ogk(z_mat)


write.table(d2, file="all.pcadapt.test.txt", quote=F, row.names=F, col.names=F)
write.table(pchisq(d2, df=K, lower.tail=F), file="all.pcadapt.pval.txt", quote=F, row.names=F, col.names=F)

pval <- read.table('all.pcadapt.pval.txt')
p.tmp <- unlist(pval)
p <- as.numeric(p.tmp)
p_adj <- p.adjust(p, "fdr", n = length(p))
write.table(p_adj, file="all.pcadapt.padj.txt", quote=F, row.names=F, col.names=F)

# make table of chrom, pos, test, pval, padj
zcat /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike.beagle.gz | awk '{print $1}' > all.snps.txt

tail -n +2 all.snps.txt > all.snps_headless.txt

sed -i 's/\.1//g' all.snps_headless.txt
# sed -i 's/\_/\,/g' all.snps_headless.chroms.txt

awk '$1>0 {print $1}' all.sites | wc -l
paste -d ' ' all.snps_headless.txt all.sites| awk '$2>0 {print $1}' > all.snps_filtered.txt



paste -d ',' all.snps_filtered.txt all.pcadapt.test.txt all.pcadapt.pval.txt all.pcadapt.padj.txt > all.pcadapt_stats_matrix.csv

grep 'NC_' all.pcadapt_stats_matrix.csv > all.pcadapt_stats_matrix.chroms.csv
sed -i 's/NC\_//g' all.pcadapt_stats_matrix.chroms.csv
sed -i 's/_/,/g' all.pcadapt_stats_matrix.chroms.csv


# plot manhattan

#!/bin/bash

#SBATCH --job-name=pcadapt_manhattan
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pcadapt_manhattan.%j

module load R

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca

Rscript all.manhattan.R

#!/usr/bin/env Rscript
#manhattan.R

library(qqman)
stats<-read.csv('all.pcadapt_stats_matrix.chroms.csv', header=FALSE)
SNP<-c(1: (nrow(stats)))
colnames(stats) <- c("CHROM","POS","TEST","PVAL", "PADJ")

mydf<-data.frame(SNP,stats)


pdf(file = "all.pcadapt_pval.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="PVAL",snp="SNP",logp=FALSE,ylab="pcadapt selection scan"), cex = 0.5)
dev.off()


pdf(file = "all.pcadapt_padj.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="PADJ",snp="SNP",logp=FALSE,ylab="pcadapt selection scan"), cex = 0.5)
dev.off()


# plot properly

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
stats<-read.csv('all.pcadapt_stats_matrix.chroms.csv', header=FALSE)
SNP<-c(1: (nrow(stats)))
colnames(stats) <- c("CHROM","POS","TEST","PVAL", "PADJ")

mydf<-data.frame(SNP,stats)

# draw it with cutoff line 


blues <- c("#4EAFAF", "#082B64")
df.tmp <- stats %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(stats, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
max(df.tmp$TEST)

png("all.pcadapt_pval.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(TEST))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 6766 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,13300)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Test statistic") +
  
  # add genome-wide sig and sugg lines
    # geom_hline(yintercept = 0.3808062) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(POS), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()







sbatch all.pcadapt_manhattan.sh
Submitted batch job 11151432


awk -F ',' '$5<0.01 {print "NC_"$1".1", $2}' all.pcadapt_stats_matrix.csv > all.pcadapt_significantsnps.txt


# top snps 

#!/bin/bash

#SBATCH --job-name=relevantgenes
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.relevantgenes.%j

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca

while read -r line;
do

chr=`awk 'BEGIN {FS = " "} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = " "} {print $2}' <<<"${line}"`

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$midpos' && $5 > '$midpos')' >> all.relevantgenes_pcadapt.txt

done < all.pcadapt_significantsnps.txt


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' all.relevantgenes_pcadapt.txt | sed 's/Name\=//g' | sort -u > all.relevantgenenames_pcadapt.txt


sbatch genelist.sh 
Submitted batch job 11097475



wc -l relevantgenenames_pcadapt.txt

grep 'CARCAR_R09680' /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff 


("GO:0008045","GO:0097720","GO:0006112","GO:0046470","GO:0005977","GO:0045637","GO:0006390","GO:0033962","GO:0009225","GO:0034982","GO:0006509","GO:0044042","GO:0006817","GO:0045840","GO:0035329","GO:0045931","GO:0033619","GO:0050806","GO:0006611","GO:0019985","GO:0006301","GO:0038127","GO:0031623","GO:0060078","GO:0005976","GO:0007173","GO:0000288","GO:0000731","GO:0000959","GO:0061157","GO:0050779","GO:0007224","GO:0034314","GO:0015698","GO:0090068","GO:0031146","GO:0006413","GO:0018345","GO:0006821","GO:0061014","GO:0016237","GO:0071897","GO:0034727","GO:0006997","GO:2000045","GO:0045787","GO:0030490","GO:0006998","GO:0098661","GO:0044773","GO:0050773","GO:0045765","GO:1901342","GO:1902476","GO:0033045","GO:0000462","GO:0042274","GO:0044772","GO:0044770","GO:1902115","GO:0006099","GO:0044843","GO:0000082","GO:0006282","GO:0007030","GO:0007416","GO:0007088","GO:0006820","GO:0000422","GO:0030010","GO:0002181","GO:1903313","GO:0007219","GO:0070972","GO:0043543","GO:0033044","GO:0061025","GO:0140053","GO:0098656","GO:0051783","GO:0006417","GO:0090630","GO:0000209","GO:0090174","GO:0051168","GO:0006906","GO:0002028","GO:0050808","GO:1901990","GO:0000725","GO:0000724","GO:0034248","GO:0000956","GO:0048488","GO:0072665","GO:0042254","GO:0043122","GO:0010608","GO:0043161","GO:0098771","GO:0006364","GO:0042158","GO:0009451","GO:0006497","GO:0043487","GO:0006400","GO:0048284","GO:0018105","GO:0055074","GO:1903050","GO:0043488","GO:0006511","GO:0000045","GO:0018209","GO:0022613","GO:0006914","GO:0061919","GO:0051169","GO:0006913","GO:0034470","GO:0008033","GO:0061013","GO:1901987","GO:0006897","GO:0019941","GO:0016236","GO:0006403","GO:0006470","GO:0016485","GO:0042273","GO:0007033","GO:0010720","GO:0010564","GO:0031503","GO:0043414","GO:0032259","GO:0006874","GO:0010948","GO:0034504","GO:0015748","GO:0043632","GO:0001510","GO:0071826","GO:0006302","GO:0051668","GO:0006281","GO:0042157","GO:0010498","GO:0048193","GO:0006310","GO:0006974","GO:0016311","GO:0006888","GO:0050658","GO:0050657","GO:0051236","GO:0007346","GO:0072657","GO:0030111","GO:0072594","GO:0006351","GO:0006402","GO:1901988","GO:0006886","GO:0016072","GO:0032774","GO:0022618","GO:0006260","GO:1990778","GO:0006261","GO:0033365","GO:0051603","GO:0006605","GO:1903008","GO:0009894","GO:0001934","GO:0030163","GO:0030258","GO:0060627","GO:0051347","GO:0009057","GO:0006486","GO:0043413","GO:0006396","GO:0034660","GO:0015931","GO:0099504","GO:0031032","GO:0099003","GO:0033674","GO:0006259","GO:0061024","GO:0031329","GO:0044270","GO:0046700","GO:0016050","GO:0090304","GO:0006869","GO:0016310","GO:0042327","GO:0016070","GO:0019439","GO:0070727","GO:0008104","GO:0009100","GO:0009059","GO:0010562","GO:0045937","GO:0070085","GO:1901136","GO:0009896","GO:0016071","GO:0033036","GO:0043412","GO:0018193","GO:0006399","GO:0031401","GO:0120035","GO:0044087","GO:0022411","GO:0051094","GO:0016192","GO:0006468","GO:0050801","GO:0010467","GO:0036211","GO:1903311","GO:0006508","GO:0034641","GO:0010256","GO:0055080","GO:1901361","GO:0051604","GO:0010876","GO:1901565","GO:0006725","GO:0045184","GO:0043170","GO:0051338","GO:0022603","GO:0043687","GO:0006139","GO:0033554","GO:0046483","GO:0009101","GO:0070647","GO:0032446","GO:0016567","GO:0034655","GO:0043085","GO:0015031","GO:0019725","GO:0051641","GO:1901360","GO:0019538","GO:0046907","GO:2000026","GO:0051649","GO:0006807","GO:0051247","GO:0043549","GO:0006325","GO:0045859","GO:0042592","GO:0051128","GO:0006397","GO:0042325","GO:0051174","GO:0019220","GO:0044249","GO:0051726","GO:0071824","GO:0051640","GO:0043043","GO:0007005","GO:1901576","GO:0009058","GO:0044238","GO:0001932","GO:0048878","GO:0034654","GO:0044093","GO:0007169","GO:0051130","GO:0019438","GO:0009056","GO:0044271","GO:0000398","GO:0000377","GO:0000375","GO:0006412","GO:1901564","GO:0050793","GO:0051179","GO:0051234","GO:0071704","GO:0035556","GO:1901575","GO:0044237","GO:0006810","GO:0051246","GO:0006457","GO:0031399","GO:0006518","GO:1901362","GO:0098657","GO:0006793","GO:0018130","GO:0008152","GO:0006796","GO:0034330","GO:0044085","GO:0007018","GO:0006338","GO:0044248","GO:0033043","GO:0023057","GO:0071705","GO:0043933","GO:0141124","GO:0010648","GO:0007015","GO:1901566","GO:0051173","GO:0050790","GO:0071702","GO:0048522","GO:0043604","GO:0070925","GO:0031325","GO:0051049","GO:0065009","GO:0010604","GO:0006644","GO:0006811","GO:0030036","GO:0006950","GO:0009968","GO:0009893","GO:0071840","GO:0022607","GO:1901135","GO:0043603","GO:0048585","GO:0007167","GO:0045935","GO:0098660","GO:0048518","GO:0016043","GO:0010646","GO:0065003","GO:0009966","GO:0023051","GO:0051254","GO:0032879","GO:0030029","GO:0055085","GO:0006996","GO:0034220","GO:0007166","GO:0044281","GO:0019637","GO:0048523","GO:0009987","GO:0044255","GO:0048519","GO:0008150","GO:0065008","GO:0031323","GO:0065007","GO:0051716","GO:0050794","GO:0051171","GO:0050789","GO:0019222","GO:0080090","GO:0060255","GO:0031326","GO:0009889","GO:0048513","GO:0007186","GO:0002682","GO:0007417","GO:0002376","GO:0007389","GO:0007218","GO:0007423","GO:0060322","GO:0007420","GO:0006955","GO:0009952","GO:0050776","GO:0045109","GO:0042742","GO:0060326","GO:0002764","GO:0002252","GO:0002250","GO:0006182","GO:0019730","GO:0003073","GO:0015986","GO:0006754","GO:0016064","GO:0002460","GO:0002449","GO:0002443","GO:0002440","GO:0070542","GO:0086002","GO:0071398","GO:0002443","GO:0002440","GO:0070542","GO:0086002","GO:0071398")



# GO analyses Gene ontology


# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

df<-read.csv("relevantgenenames_snps_top.95.txt", header=TRUE)

# top 0.05 of snps
go_list <- c("GO:0044550","GO:0048009","GO:0019748","GO:0006509","GO:0090257","GO:0043491","GO:0001822","GO:0035418","GO:0006937","GO:1905475","GO:1902414","GO:0051209","GO:0007157","GO:0048514","GO:0001525","GO:0035335","GO:0099072","GO:0051283","GO:0001944","GO:0035239","GO:0035295","GO:0001568","GO:0090263","GO:0051282","GO:0032922","GO:0098661","GO:1902476","GO:0030177","GO:0098656","GO:0007411","GO:0097485","GO:0006821","GO:0097553","GO:0046856","GO:0006820","GO:0015698","GO:0006470","GO:0043122","GO:0098742","GO:0090630","GO:0006469","GO:0099177","GO:0050804","GO:0033674","GO:0051347","GO:0007169","GO:0030335","GO:0098609","GO:0043547","GO:2000147","GO:0040017","GO:0048667","GO:0007409","GO:0050808","GO:0051056","GO:0000902","GO:0030334","GO:0070588","GO:2000145","GO:0040012","GO:0072359","GO:0051241","GO:0016311","GO:0048812","GO:0048858","GO:0120039","GO:0051338","GO:0048666","GO:0043408","GO:0043549","GO:0061564","GO:0042327","GO:0043085","GO:1902531","GO:0010562","GO:0045937","GO:0043484","GO:0031175","GO:0051345","GO:0007155","GO:0042325","GO:0031503","GO:0044093","GO:0007167","GO:0034330","GO:0043087","GO:0006816","GO:0008284","GO:0051174","GO:0019220","GO:0098771","GO:0023051","GO:0010646","GO:0032989","GO:0050801","GO:0055080","GO:0030182","GO:0009653","GO:0048699","GO:1902533","GO:0009966","GO:0006468","GO:0048646","GO:0050790","GO:0001932","GO:0009101","GO:0010647","GO:0023056","GO:0099536","GO:0048468","GO:0120036","GO:0022008","GO:0030030","GO:0032879","GO:0016310","GO:0065009","GO:0007275","GO:0051239","GO:0031399","GO:0009967","GO:0007399","GO:0099537","GO:0050793","GO:0007268","GO:0098916","GO:0048731","GO:0016477","GO:0051247","GO:0007166","GO:0098660","GO:0048870","GO:0048856","GO:0035556","GO:0007267","GO:0032501","GO:0048583","GO:0006796","GO:0006793","GO:0048584","GO:0051246","GO:0034220","GO:0030001","GO:0032502","GO:0006811","GO:0036211","GO:0048522","GO:0023052","GO:0048518","GO:0007154","GO:0043412","GO:0007165","GO:0009893","GO:0031325","GO:0065007","GO:0050794","GO:0050789","GO:0051179","GO:0051641","GO:0030154","GO:0048869","GO:0006810","GO:0019538","GO:0016043","GO:0051234","GO:0071840","GO:0051716","GO:0008150","GO:0051171","GO:0080090","GO:0031323","GO:0009987","GO:1901564","GO:0060255","GO:0019222","GO:0050896","GO:0002682","GO:0002376","GO:0002764","GO:0050776","GO:0006955","GO:0002252","GO:0002376","GO:0002764","GO:0050776","GO:0006955","GO:0002252")



# run locally
# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

# df<-read.csv("relevantgenenames_snps_top.95.txt", header=TRUE)

# top 0.05 of snps
go_list <- c("GO:0044550","GO:0048009","GO:0019748","GO:0006509","GO:0090257","GO:0043491","GO:0001822","GO:0035418","GO:0006937","GO:1905475","GO:1902414","GO:0051209","GO:0007157","GO:0048514","GO:0001525","GO:0035335","GO:0099072","GO:0051283","GO:0001944","GO:0035239","GO:0035295","GO:0001568","GO:0090263","GO:0051282","GO:0032922","GO:0098661","GO:1902476","GO:0030177","GO:0098656","GO:0007411","GO:0097485","GO:0006821","GO:0097553","GO:0046856","GO:0006820","GO:0015698","GO:0006470","GO:0043122","GO:0098742","GO:0090630","GO:0006469","GO:0099177","GO:0050804","GO:0033674","GO:0051347","GO:0007169","GO:0030335","GO:0098609","GO:0043547","GO:2000147","GO:0040017","GO:0048667","GO:0007409","GO:0050808","GO:0051056","GO:0000902","GO:0030334","GO:0070588","GO:2000145","GO:0040012","GO:0072359","GO:0051241","GO:0016311","GO:0048812","GO:0048858","GO:0120039","GO:0051338","GO:0048666","GO:0043408","GO:0043549","GO:0061564","GO:0042327","GO:0043085","GO:1902531","GO:0010562","GO:0045937","GO:0043484","GO:0031175","GO:0051345","GO:0007155","GO:0042325","GO:0031503","GO:0044093","GO:0007167","GO:0034330","GO:0043087","GO:0006816","GO:0008284","GO:0051174","GO:0019220","GO:0098771","GO:0023051","GO:0010646","GO:0032989","GO:0050801","GO:0055080","GO:0030182","GO:0009653","GO:0048699","GO:1902533","GO:0009966","GO:0006468","GO:0048646","GO:0050790","GO:0001932","GO:0009101","GO:0010647","GO:0023056","GO:0099536","GO:0048468","GO:0120036","GO:0022008","GO:0030030","GO:0032879","GO:0016310","GO:0065009","GO:0007275","GO:0051239","GO:0031399","GO:0009967","GO:0007399","GO:0099537","GO:0050793","GO:0007268","GO:0098916","GO:0048731","GO:0016477","GO:0051247","GO:0007166","GO:0098660","GO:0048870","GO:0048856","GO:0035556","GO:0007267","GO:0032501","GO:0048583","GO:0006796","GO:0006793","GO:0048584","GO:0051246","GO:0034220","GO:0030001","GO:0032502","GO:0006811","GO:0036211","GO:0048522","GO:0023052","GO:0048518","GO:0007154","GO:0043412","GO:0007165","GO:0009893","GO:0031325","GO:0065007","GO:0050794","GO:0050789","GO:0051179","GO:0051641","GO:0030154","GO:0048869","GO:0006810","GO:0019538","GO:0016043","GO:0051234","GO:0071840","GO:0051716","GO:0008150","GO:0051171","GO:0080090","GO:0031323","GO:0009987","GO:1901564","GO:0060255","GO:0019222","GO:0050896","GO:0002682","GO:0002376","GO:0002764","GO:0050776","GO:0006955","GO:0002252","GO:0002376","GO:0002764","GO:0050776","GO:0006955","GO:0002252")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

scatterPlot(simMatrix, reducedTerms)

# "Anatomical structure involved in morphogenesis" seems relevant to bill morphology
heatmapPlot(reducedTerms,
            annotateParent=TRUE,
            fontsize=6)

treemapPlot(reducedTerms)

