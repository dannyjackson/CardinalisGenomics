# running interactively because it is a quick script
chmod +x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh

# FST
# Define species, environments, and window sizes
# Identify genes in windows
species=( "noca" "pyrr" )
window_sizes=( 50000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_fst.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst_${win}.outlier.csv \
    -n ${sp} \
    -m fst \
    -w ${win}

    done
done

# Identify nonsynonymous and missense snps in outlier windows
# filter snpEff annotated genome with outlier regions
species=( "noca" "pyrr" )
window_sizes=( 50000 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst_${win}.outlier.csv ${win} /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst.${win}.outlier.regions.bed /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

        IN_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst.${win}.outlier.regions.bed
        VCF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf
        VCF_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst.${win}.outlier.ann.snps.vcf
        GENE_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/${sp}urban_${sp}rural.fst.${win}.outlier.snpeff.genes.txt

        grep '#' ${VCF} > ${VCF_OUT}
        bedtools intersect -a ${VCF} -b ${IN_FILE} -wa >> ${VCF_OUT}

        grep -e 'missense' -e 'nonsense' ${VCF_OUT} | awk '{print $8}' | awk -F"|" '{print $5}' | sort -u > ${GENE_OUT}
    done
done

# further annotate by high fst snps

# find snps with more than 2 individuals varying
bcftools view -i 'AC>0' pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf.gz 

-Oz -o pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.filtered.vcf.gz


bcftools +fill-tags pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf.gz -- -t AC | bcftools view -i 'AC>2 & AC<24' -Oz -o pyrrurban_pyrrrural.filtered.vcf.gz


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT ]\n' pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf.gz | \
awk '{count=0; for (i=5; i<=NF; i++) if ($i ~ /1/) count++; if (count > 2) print $0}' | bcftools view -i 'AC>2 & AC<24' 

-Oz -o pyrrurban_pyrrrural.filtered.vcf.gz




# 



# just look at maf > 0.2 or something
# pyrr
bcftools view -e'MAF>0.3' pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf.gz | bcftools view -e'MAF>0.8' > pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.maf.vcf

grep -e 'missense' -e 'nonsense' pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.maf.vcf | awk '{print $8}' | awk -F"|" '{print $5}' | grep -v 'LOC' | sort -u > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrr.fst.50kb.snpeff.txt


python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_1.outlier.csv 1 /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst.1.outlier.snps.bed /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

awk 'BEGIN {OFS="\t"}; {print $1, $2 - 1, $2}' pyrrurban_pyrrrural.fst.1.outlier.snps.bed > pyrrurban_pyrrrural.fst.1.outlier.regions.bed

awk '{if(/#/)print;else exit}' /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf > pyrrurban_pyrrrural.fst.50000.outlier.ann.vcf

bedtools intersect -a /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf -b pyrrurban_pyrrrural.fst.50000.outlier.regions.bed -wa >> pyrrurban_pyrrrural.fst.50000.outlier.ann.vcf

awk '{if(/#/)print;else exit}' pyrrurban_pyrrrural.fst.50000.outlier.ann.vcf > pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf

bedtools intersect -a tmp.vcf -b pyrrurban_pyrrrural.fst.1.outlier.regions.bed -wa > pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf 


grep -e 'MODERATE' -e 'HIGH' pyrrurban_pyrrrural.fst.50000.outlier.ann.snps.vcf | awk '{print $8}' | awk -F"|" '{print $5}' | grep -v 'LOC' | sort -u > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrr.fst.50kb.snpeff.txt

grep -e 'chromosome_number_variation' -e 'exon_loss_variant' -e 'frameshift_variant' -e 'rare_amino_acid_variant' -e 'splice_acceptor_variant' -e 'splice_donor_variant' -e 'start_lost' -e 'stop_gained' -e 'stop_lost' -e 'transcript_ablation' -e '3_prime_UTR_truncation & exon_loss' -e '5_prime_UTR_truncation & exon_loss_variant' -e 'coding_sequence_variant' -e 'conservative_inframe_deletion' -e 'conservative_inframe_insertion' -e 'disruptive_inframe_deletion' -e 'disruptive_inframe_insertion' -e 'missense_variant' -e 'regulatory_region_ablation' -e 'splice_region_variant' -e 'TFBS_ablation' pyrrurban_pyrrrural.fst.50000.outlier.ann.vcf | awk '{print $8}' | awk -F"|" '{print $5}' | grep -v 'LOC' | sort -u > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrr.fst.50kb.snpeff.txt

cat  pyrrurban_pyrrrural.fst.50000.outlier.ann.vcf | awk '{print $8}' | awk -F"|" '{print $5}' | grep -v 'LOC' | sort -u | wc -l 





# FDR adjust p values 











pyrr.fst.50kb.snpeff.txt (raw P-value)

pvals <- c(0.001, 0.02, 0.05, 0.10, 0.15, 0.25, 0.40)

# FDR correction
qvals <- p.adjust(pvals, method = "fdr")

# Filter q-values greater than 0.2
qvals_filtered <- qvals[qvals > 0.2]

# Print results
qvals_filtered

# DXY'/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_500000.snps.outlier.csv
# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 50000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_dxy.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy_${win}.outlier.csv \
    -n ${sp} \
    -m dxy \
    -w ${win}

    done
done


# Identify nonsynonymous and missense snps in outlier windows
# filter snpEff annotated genome with outlier regions
species=( "noca" "pyrr" )
window_sizes=( 50000 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy_${win}.outlier.csv ${win} /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy.${win}.outlier.regions.bed /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

        IN_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy.${win}.outlier.regions.bed
        VCF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf
        VCF_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy.${win}.outlier.ann.snps.vcf
        GENE_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/${sp}urban_${sp}rural.dxy.${win}.outlier.snpeff.genes.txt

        grep '#' ${VCF} > ${VCF_OUT}
        bedtools intersect -a ${VCF} -b ${IN_FILE} -wa >> ${VCF_OUT}

        grep -e 'missense' -e 'nonsense' ${VCF_OUT} | awk '{print $8}' | awk -F"|" '{print $5}' | sort -u > ${GENE_OUT}
    done
done


# RAiSD

species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 50 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_raisd.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${sp}/${sp}.raisd_${win}.outlier.csv \
    -n ${sp} \
    -m raisd \
    -w ${win}

    done
done

# RAiSD
chmod +x ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh
species=( "noca" "pyrr" )
window_sizes=( 50 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh \
    -p ~/programs/CardinalisGenomics/params_base.sh \
    -i ${sp}urban \
    -q ${sp}rural \
    -m raisd \
    -w ${win}

    done
done


# Identify nonsynonymous and missense snps in outlier windows
# filter snpEff annotated genome with outlier regions
species=( "noca" "pyrr" )
window_sizes=( 50 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${sp}urban/${sp}urban.raisd_${win}.outlier.csv ${win} /xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${sp}urban/${sp}urban_.raisd.${win}.outlier.regions.bed /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

        IN_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/raisd/${sp}urban/${sp}urban_.raisd.${win}.outlier.regions.bed
        VCF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf
        VCF_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban.raisd.${win}.outlier.ann.snps.vcf
        GENE_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/${sp}urban_.raisd.${win}.outlier.snpeff.genes.txt

        grep '#' ${VCF} > ${VCF_OUT}
        bedtools intersect -a ${VCF} -b ${IN_FILE} -wa >> ${VCF_OUT}

        grep -e 'missense' -e 'nonsense' ${VCF_OUT} | awk '{print $8}' | awk -F"|" '{print $5}' | sort -u > ${GENE_OUT}
    done
done






















# 
ls /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/noca* | grep -v 'raisd' > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/filelist.txt
ls /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/pyrr* | grep -v 'raisd' >> /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/filelist.txt
ls *raisd*unique* >> /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/filelist.txt

# Define the file containing genes to exclude
EXCLUDE_FILE="excludedgenes.3depth28.prop90.txt"


# Define the file that contains the list of gene files
FILE_LIST="filelist.txt"

# Check if the required files exist
if [[ ! -f "$EXCLUDE_FILE" ]]; then
    echo "Error: Exclude file '$EXCLUDE_FILE' not found."
    exit 1
fi

if [[ ! -f "$FILE_LIST" ]]; then
    echo "Error: File list '$FILE_LIST' not found."
    exit 1
fi

# Process each file listed in file_list.txt
while IFS= read -r FILE; do
    if [[ -f "$FILE" ]]; then
        OUTPUT_FILE="${FILE%.txt}_filtered.txt"
        cp "$FILE" "$OUTPUT_FILE"  # Make a copy to modify

        # Iterate through each gene in exclude_genes.txt and remove it from the file
        while IFS= read -r GENE; do
            sed -i "/^$GENE$/d" "$OUTPUT_FILE"
        done < "$EXCLUDE_FILE"

        echo "Filtered file created: $OUTPUT_FILE"
    else
        echo "Warning: File not found - $FILE"
    fi
done < "$FILE_LIST"

mv *_filtered.txt ../final_gene_lists

cd ../final_gene_lists

cat pyrr* | sort -u > pyrr.all.genenames_filtered.txt
cat noca* | sort -u > noca.all.genenames_filtered.txt

rm /home/u15/dannyjackson/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/*
cp * /home/u15/dannyjackson/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists

# Analyze overlap between methods
library(ggvenn)

df_noca_fst<-read.csv("noca.fst.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_noca_raisd<-read.csv("nocaurban.raisd.50kb.unique.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_pyrr_fst<-read.csv("pyrr.fst.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_pyrr_raisd<-read.csv("pyrrurban.raisd.50kb.unique.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_owl_tmp<-read.csv("BurrowingOwlGenes.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_owl<-unique(df_owl_tmp)
df_lizard<-read.csv("LizardUrbanGenes.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_lizardmorph<-read.csv("LizardMorphologyGenes.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_tit<-read.csv("GreatTitGenes.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_noca_fstraisd <- unique(c(df_noca_fst, df_noca_raisd))
df_pyrr_fstraisd <- unique(c(df_pyrr_fst, df_pyrr_raisd))

# compare with owl and great tit
gene_sets <- list(
  noca = df_noca_fstraisd,
  pyrr = df_pyrr_fstraisd,
  owl = df_owl,
  tit = df_tit
)

intersect_owl <- Reduce(intersect, list(df_noca_fstraisd, df_pyrr_fstraisd, df_owl))
write(intersect_owl, file="intersection.owl.csv")

intersect_owl_noca = intersect(df_noca_fstraisd, df_owl)
write(intersect_owl_noca, file="intersection.owl_noca.csv")

intersect_owl_pyrr = intersect(df_pyrr_fstraisd, df_owl)
write(intersect_owl_pyrr, file="intersection.owl_pyrr.csv")

pdf(file = "intersection.owl.tit.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("noca", "pyrr", "owl", "tit"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


# identify genes in overlapping bits
# noca owl
noca_owl = intersect(df_noca_fstraisd, df_owl)
write(noca_owl, file="intersection.noca_owl.csv")
# noca owl tit
noca_owl_tit <- Reduce(intersect, list(df_noca_fstraisd, df_owl, df_tit))
write(noca_owl_tit, file="intersection.noca_owl_tit.csv")
# noca tit
noca_tit = intersect(df_noca_fstraisd, df_tit)
write(noca_tit, file="intersection.noca_tit.csv")
# noca pyrr owl
noca_pyrr_owl <- Reduce(intersect, list(df_noca_fstraisd, df_owl, df_pyrr_fstraisd))
write(noca_pyrr_owl, file="intersection.noca_pyrr_owl.csv")
# pyrr owl
pyrr_owl = intersect(df_pyrr_fstraisd, df_owl)
write(pyrr_owl, file="intersection.pyrr_owl.csv")
# pyrr tit
pyrr_tit = intersect(df_pyrr_fstraisd, df_tit)
write(pyrr_tit, file="intersection.pyrr_tit.csv")
# owl tit
owl_tit = intersect(df_owl, df_tit)
write(owl_tit, file="intersection.owl_tit.csv")


# genes in 2/4 lists
fst = intersect(df_noca_fst, df_pyrr_fst)

# genes in 3/4 lists
intersect_nocaall_pyrrfst <- Reduce(intersect, list(df_noca_fst, df_noca_raisd, df_pyrr_fst))


# genes in either cardinal species + one other species
overlapping<-unique(c(noca_owl, noca_tit, pyrr_owl, pyrr_tit))
write(overlapping, file="intersection.overlapping.csv")











# compare between cardinals
gene_sets <- list(
  noca_fst = df_noca_fst,
  noca_raisd = df_noca_raisd,
  pyrr_fst = df_pyrr_fst,
  pyrr_raisd = df_pyrr_raisd
)




# Compute overlaps
noca_fst_raisd = intersect(df_noca_fst, df_noca_raisd)
pyrr_fst_raisd = intersect(df_pyrr_fst, df_pyrr_raisd)
noca_pyrr = intersect(df_noca_fstraisd, df_pyrr_fstraisd)
noca_pyrr_allanalyses = intersect(noca_fst_raisd, pyrr_fst_raisd)
raisd = intersect(df_noca_raisd, df_pyrr_raisd)
fst = intersect(df_noca_fst, df_pyrr_fst)

# genes in 3/4 lists
intersect_nocaall_pyrrfst <- Reduce(intersect, list(df_noca_fst, df_noca_raisd, df_pyrr_fst))
intersect_nocaall_pyrrraisd <- Reduce(intersect, list(df_noca_fst, df_noca_raisd, df_pyrr_raisd))
intersect_nocafst_pyrrall <- Reduce(intersect, list(df_pyrr_fst, df_pyrr_raisd, df_noca_fst))
intersect_nocaraisd_pyrrall <- Reduce(intersect, list(df_pyrr_fst, df_pyrr_raisd, df_noca_raisd))

intersect_3of4 <- unique(vctrs::vec_c(intersect_nocaall_pyrrfst, intersect_nocaall_pyrrraisd, intersect_nocafst_pyrrall, intersect_nocaraisd_pyrrall))
write(intersect_3of4, file="intersection.3of4.csv")

write(raisd, file="intersection.raisd.csv")
write(fst, file="intersection.fst.csv")

write(noca_fst_raisd, file="intersection.noca_fst_raisd.csv")
write(pyrr_fst_raisd, file="intersection.pyrr_fst_raisd.csv")
write(noca_pyrr, file="intersection.noca_pyrr.csv")
write(noca_pyrr_allanalyses, file="intersection.noca_pyrr_allanalyses.csv")

pdf(file = "intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("noca_fst", "noca_raisd", "pyrr_raisd",  "pyrr_fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

# test for significance in overlap
# noca vs pyrr
background<-read.csv("background_genelist.full.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
N=length(background) # 16264
k = 132
m = length(df_noca_fstraisd) # 667
n = length(df_pyrr_fstraisd) # 895

# fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))

fisher.test(matrix(c(132, 684-132, 924-132, 16264-684-924+132), nrow = 2))

Fisher's Exact Test for Count Data

data:  matrix(c(132, 684 - 132, 924 - 132, 16264 - 684 - 924 + 132), nrow = 2)
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
3.615506 5.483880
sample estimates:
odds ratio 
4.464239 

# noca vs owl
N=length(background) # 16264
k = length(intersect(df_noca_fstraisd, df_owl)) # 21

m = length(df_noca_fstraisd) # 684
n = length(df_owl) # 174

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 9.369e-06
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.960531 5.238269
sample estimates:
odds ratio 
  3.280927 

# noca vs tit
N=length(background) # 16264
k = length(intersect(df_noca_fstraisd, df_tit)) # 8

m = length(df_noca_fstraisd) # 684
n = length(df_tit) # 73

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 0.00987
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.196662 6.094150
sample estimates:
odds ratio 
  2.900475 


# pyrr vs owl
N=length(background) # 16264
k = length(intersect(df_pyrr_fstraisd, df_owl)) # 19

m = length(df_pyrr_fstraisd) # 895
n = length(df_owl) # 174

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))

        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 0.003943
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.241715 3.458386
sample estimates:
odds ratio 
  2.128775 

# pyrr vs tit
N=length(background) # 16264
k = length(intersect(df_pyrr_fstraisd, df_tit)) # 5

m = length(df_pyrr_fstraisd) # 895
n = length(df_tit) # 73

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 0.6012
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.396738 3.108486
sample estimates:
odds ratio 
  1.264106 


# owl vs tit
N=length(background) # 16264
k = length(intersect(df_owl, df_tit)) # 5

m = length(df_owl) # 895
n = length(df_tit) # 73

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


# p adjust
# noca_pyrr, noca_owl, noca_tit, pyrr_owl, pyrr_tit
p = c(2.2e-16, 9.369e-06, 0.00987, 0.003943, 0.6012)
options(scipen = 999)
p.adjust(p, method = "fdr")
0.0000000000000011 # noca pyrr
0.0000234225000000 # noca owl
0.0123375000000000 # noca tit
0.0065716666666667 # pyrr owl
0.6012000000000000 # pyrr tit

# noca_pyrr, noca_owl, noca_tit, pyrr_owl, pyrr_tit, owl_tit
p = c(2.2e-16, 9.369e-06, 0.00987, 0.003943, 0.6012, 0.007766)
options(scipen = 999)
p.adjust(p, method = "fdr")
[1] 0.00000000000000132 0.00002810700000000 0.01184400000000000
[4] 0.00788600000000000 0.60119999999999996 0.01164900000000000
options(scipen = 0)
p.adjust(p, method = "fdr")
# GO Term analysis
# Remove extraneous information from Panther output and only leave the GO terms

for file 
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.noca_all.nocorrection.txt > intersection.noca_all.nocorrection.GO.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.pyrr_all.nocorrection.txt > intersection.pyrr_all.nocorrection.GO.txt

grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.pyrr_all.nocorrection.txt > intersection.pyrr_all.nocorrection.GO.txt



grep -o 'GO:[0-9]\{7\}' pantheroutput/both.all.GO.nocorrection.txt > both.all.GO.nocorrection.txt