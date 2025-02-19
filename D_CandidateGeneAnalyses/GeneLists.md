# running interactively because it is a quick script
chmod +x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh

# FST
# Define species, environments, and window sizes
# Identify genes in windows
species=( "noca" "pyrr" )
window_sizes=( 25000 )

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
for TAXA in pyrr noca; do
    python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${TAXA}urban_${TAXA}rural/${TAXA}urban_${TAXA}rural.fst_25000.outlier.csv 25000 /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${TAXA}urban_${TAXA}rural/${TAXA}urban_${TAXA}rural.fst.25000.outlier.regions.bed /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

    IN_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${TAXA}urban_${TAXA}rural/${TAXA}urban_${TAXA}rural.fst.25000.outlier.regions.bed
    GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
    VCF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf
    VCF_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${TAXA}urban_${TAXA}rural/${TAXA}urban_${TAXA}rural.fst.25000.outlier.ann.snps.vcf
    GENE_OUT=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/gene_names/${TAXA}urban_${TAXA}rural.fst.25000.outlier.snpeff.genes.txt

    grep '#' ${VCF} > ${VCF_OUT}
    bedtools intersect -a ${VCF} -b ${IN_FILE} -wa >> ${VCF_OUT}

    grep -e 'missense' -e 'nonsense' ${VCF_OUT} | awk '{print $8}' | awk -F"|" '{print $5}' | sort -u > ${GENE_OUT}
done



# DXY'/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_500000.snps.outlier.csv
# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 25000 5000 )

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

# Tajima
# Define species, environments, and window sizes
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 50000 10000 5000 1000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_tajima.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/Tajima/${sp}/${sp}.Tajima_${win}.outlier.csv \
    -n ${sp} \
    -m Tajima \
    -w ${win}

    done
done


chmod +x ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh
species=( "noca" "pyrr" )
window_sizes=( 50000 10000 5000 1000 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh \
    -p ~/programs/CardinalisGenomics/params_base.sh \
    -i ${sp}urban \
    -q ${sp}rural \
    -m Tajima \
    -w ${win}

    done
done

# RAiSD

species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 10 20 50 100 500 1000 )

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
window_sizes=( 10 20 50 100 500 1000 )

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

