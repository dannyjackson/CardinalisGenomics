# running interactively because it is a quick script
chmod +x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh

# FST
# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 500000 10000 1000 )

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
NC_044582.1:
19879791-19882847
        19879791-19882847
19878000-19879000
NC_044582.1:19878000-19883000
19852994
19878927

python ${SCRIPTDIR}/Genomics-Main/general_scripts/outlier_to_bed.py ${IN_FILE} ${WIN} ${OUT_FILE} ${CHR_FILE}

bedtools intersect -a /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff -b pyrr.fst.1000kb.bed -wa | grep 'NC_044582'
# plot chr 12 pyrrurban_pyrrrural.fst.25000.Ztransformed.csv
head -n 1 pyrrurban_pyrrrural.fst.25000.Ztransformed.csv > chr12.25000.Ztrfm.csv

awk -F'\t' '$2 == 12' pyrrurban_pyrrrural.fst.50000.Ztransformed.csv >> chr12.50000.Ztrfm.csv

awk -F'\t' '$2 == 12' pyrrurban_pyrrrural.fst.25000.Ztransformed.csv >> chr12.25000.Ztrfm.csv

awk -F'\t' '$2 == 12' pyrrurban_pyrrrural.fst.10000.Ztransformed.csv >> chr12.10000.Ztrfm.csv


grep "\"12\"" pyrrurban_pyrrrural.fst_.outlier.csv 
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"

12388,19823114)(19800000,19825000)     12      19812500        5       0.016027        -0.3674201366344778     0.19155460291816487
(15414,15430)(19828908,19849057)(19825000,19850000)     12      19837500        18      0.184089        2.17035355893555        1.8241973225707662
(15431,15438)(19852979,19870063)(19850000,19875000)     12      19862500        9       0.340516        4.532436158672151       5.535307056197194
(15439,15448)(19875595,19899392)(19875000,19900000)     12      19887500        11      0.295764        3.856670974251225       4.240552337322197
(15449,15461)(19901097,19921304)(19900000,19925000)     12      19912500        14      0.087106        0.7058886214870607      0.6195559265322832


# DXY'/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_500000.snps.outlier.csv
# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 500000 10000 1000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_dxy.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy_${win}.snps.outlier.csv \
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

