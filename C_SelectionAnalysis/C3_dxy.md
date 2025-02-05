# C3_dxy.md

# Run once per species to generate all files

# Define species, environments, and window sizes
declare -A time_limits=( [500000]=1:00:00 [10000]=1:00:00 [1000]=1:00:00 [1]=10:00:00 )
species=( "noca" "pyrr" )
window_sizes=( 500000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        [ $win -eq 1 ] && mem_limit="100gb"

        sbatch --account=mcnew \
               --job-name=dxy_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.dxy_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=$time_limit \
               --mem=$mem_limit \
               ~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_dxy.sh \
               -w $win -s $win
    done
done

# Iterate over several window sizes

# Define species, environments, and window sizes
declare -A time_limits=( [500000]=1:00:00 [10000]=1:00:00 [1000]=1:00:00 [1]=10:00:00 )
declare -A mem_limits=( [500000]="10gb" [10000]="10gb" [1000]="10gb" [1]="100gb" )
species=( "noca" "pyrr" )
window_sizes=( 10000 1000 1 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}

        sbatch --account=mcnew \
               --job-name=dxy_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.dxy_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=$time_limit \
               --mem=$mem_limit \
               ~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_dxy.sh \
               -w $win -s $win
    done
done


















# compute 500,000bp windows and snp analysis 
# northern cardinals
sbatch --account=mcnew \
--job-name=dxy_50k_noca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_50k_noca.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/noca_params_dxy.sh -w 500000 -s 500000 
# pyrrhuloxia
sbatch --account=mcnew \
--job-name=dxy_50k_pyrr \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_50k_pyrr.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/pyrr_params_dxy.sh -w 500000 -s 500000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt


git push https://dannyjackson@github.com/dannyjackson/CardinalisGenomics.git

# compute 10,000bp windows 
# northern cardinals
# 12091102
sbatch --account=mcnew \
--job-name=dxy_10k_noca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_10k_noca.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/noca_params_dxy.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# pyrrhuloxia
# 12091103
sbatch --account=mcnew \
--job-name=dxy_10k_pyrr \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_10k_pyrr.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/pyrr_params_dxy.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt


# compute 1,000bp windows 
# northern cardinals
# 12091150
sbatch --account=mcnew \
--job-name=dxy_1k_noca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1k_noca.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/noca_params_dxy.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt

# pyrrhuloxia
# 12091151
sbatch --account=mcnew \
--job-name=dxy_1k_pyrr \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1k_pyrr.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/dxy/dxy.sh -p ~/programs/CardinalisGenomics/pyrr_params_dxy.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt