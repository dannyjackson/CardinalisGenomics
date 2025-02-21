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
# NOCA
df_noca_fst<-read.csv("noca.fst.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_noca_dxy<-read.csv("noca.dxy.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_noca_raisd<-read.csv("nocaurban.raisd.50kb.unique.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  noca_fst = df_noca_fst,
  noca_dxy = df_noca_dxy,
  noca_raisd = df_noca_raisd
)

# Compute overlaps
noca_fst_dxy = intersect(df_noca_fst, df_noca_dxy)
noca_fst_raisd = intersect(df_noca_fst, df_noca_raisd)
noca_dxy_raisd = intersect(df_noca_dxy, df_noca_raisd)
noca_fst_dxy_raisd <- Reduce(intersect, list(df_noca_fst, df_noca_dxy, df_noca_raisd))



write(noca_fst_dxy, file="intersection.noca_fst_dxy.csv")
write(noca_fst_raisd, file="intersection.noca_fst_raisd.csv")
write(noca_dxy_raisd, file="intersection.noca_dxy_raisd.csv")
write(noca_fst_dxy_raisd, file="intersection.noca_fst_dxy_raisd.csv")

pdf(file = "noca_windows.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("noca_fst", "noca_dxy", "noca_raisd"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

# PYRR

df_pyrr_fst<-read.csv("pyrr.fst.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_pyrr_dxy<-read.csv("pyrr.dxy.50000kb.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_pyrr_raisd<-read.csv("pyrrurban.raisd.50kb.unique.genenames_filtered.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  pyrr_fst = df_pyrr_fst,
  pyrr_dxy = df_pyrr_dxy,
  pyrr_raisd = df_pyrr_raisd
)

# Compute overlaps
pyrr_fst_dxy = intersect(df_pyrr_fst, df_pyrr_dxy)
pyrr_fst_raisd = intersect(df_pyrr_fst, df_pyrr_raisd)
pyrr_dxy_raisd = intersect(df_pyrr_dxy, df_pyrr_raisd)
pyrr_fst_dxy_raisd <- Reduce(intersect, list(df_pyrr_fst, df_pyrr_dxy, df_pyrr_raisd))



write(pyrr_fst_dxy, file="intersection.pyrr_fst_dxy.csv")
write(pyrr_fst_raisd, file="intersection.pyrr_fst_raisd.csv")
write(pyrr_dxy_raisd, file="intersection.pyrr_dxy_raisd.csv")
write(pyrr_fst_dxy_raisd, file="intersection.pyrr_fst_dxy_raisd.csv")

pdf(file = "pyrr_windows.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("pyrr_fst", "pyrr_dxy", "pyrr_raisd"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


# subsets, both species
cat intersection.pyrr* | sort -u > intersection.pyrr_all.csv
cat intersection.noca* | sort -u > intersection.noca_all.csv

library(ggvenn)
df_pyrr_all<-read.csv("intersection.pyrr_all.csv", header=TRUE, stringsAsFactors = FALSE)[[1]]
df_noca_all<-read.csv("intersection.noca_all.csv", header=TRUE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  pyrr_all = df_pyrr_all,
  noca_all = df_noca_all
)

# Compute overlaps
both = intersect(df_pyrr_all, df_noca_all)

write(both, file="intersection.both.csv")

pdf(file = "both_windows.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("pyrr_all", "noca_all"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()






# GO Term analysis
# Remove extraneous information from Panther output and only leave the GO terms
grep -o 'GO:[0-9]\{6\}' pantheroutput/intersection.noca_all.nocorrection.txt > intersection.noca_all.nocorrection.GO.txt
grep -o 'GO:[0-9]\{6\}' pantheroutput/intersection.pyrr_all.nocorrection.txt > intersection.pyrr_all.nocorrection.GO.txt