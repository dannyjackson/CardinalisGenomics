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
    -p ~/programs/CardinalisGenomics/noca_params_fst.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst_${win}.snps.outlier.csv \
    -n ${sp} \
    -m fst \
    -w ${win}

    done
done


# DXY
# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 500000 10000 1000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/noca_params_dxy.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/${sp}urban_${sp}rural/${sp}urban_${sp}rural.dxy_${win}.snps.outlier.csv \
    -n ${sp} \
    -m dxy \
    -w ${win}

    done
done

# Tajima
# Define species, environments, and window sizes
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 50000 10000 1000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/noca_params_dxy.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/Tajima/${sp}/${sp}.Tajima_${win}.snps.outlier.csv \
    -n ${sp} \
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
    -p ~/programs/CardinalisGenomics/params_base.sh
    -i ${sp}urban \
    -q ${sp}rural \
    -m raisd \
    -w ${win}

    done
done

