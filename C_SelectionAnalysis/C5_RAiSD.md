# RAiSD.md 

# merge vcfs

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocaurban.txt

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocarural.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrurban.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrrural.txt




sbatch --account=mcnew \
    --job-name=raisd.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.raisd.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=10:00:00 \
    ~/programs/raisd_preprocessing.sh 

# 
SPECIES=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
WINDOW_SIZES=( 1000 )

# Iterate over each combination
for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    sbatch --account=mcnew \
        --job-name=raisd_${WIN}_${SP} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.raisd_${WIN}_${SP}.%j \
        --nodes=1 \
        --time=1:00:00 \
        ~/programs/CardinalisGenomics/Genomics-Main/raisd/raisd.sh -p ~/programs/CardinalisGenomics/Genomics-Main/raisd/params_raisd.sh -n ${SP} -w ${WIN}
    done 
done

~/programs/CardinalisGenomics/Genomics-Main/raisd/raisd.sh -p ~/programs/CardinalisGenomics/pyrrurban_params_raisd.sh -n pyrrurban -w 10000

Rscript ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/manhattanplot.r "/xdisk/mcnew/dannyjackson/cardinals" "#2d6a4f" "#74c69d" "0.001" "pyrrurban/pyrrurban.raisd_10000.Ztransformed.csv" "10000" "raisd" "pyrrurban"
