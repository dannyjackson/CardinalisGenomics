# C4_tajima.md

# run once per population to generate all files

# Define species, environments, and window sizes
declare -A time_limits=( [500000]=1:00:00 [10000]=1:00:00 [1000]=1:00:00 [1]=10:00:00 )
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 500000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        [ $win -eq 1 ] && mem_limit="100gb"

        sbatch --account=mcnew \
               --job-name=tajima_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.tajima_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=$time_limit \
               --mem=$mem_limit \
               ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_tajima.sh \
               -w $win -s $win
    done
done


# iterate over additional window sizes

# Define species, environments, and window sizes
declare -A time_limits=( [500000]=1:00:00 [10000]=1:00:00 [1000]=1:00:00 [1]=10:00:00 )
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 10000 1000 1 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        [ $win -eq 1 ] && mem_limit="100gb"

        sbatch --account=mcnew \
               --job-name=tajima_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.tajima_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=$time_limit \
               --mem=$mem_limit \
               ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_tajima.sh \
               -w $win -s $win
    done
done



# troubleshoot snps sbatch settings

sbatch --account=mcnew \
        --job-name=tajima_1_nocaurban \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1_nocaurban.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=48:00:00 \
        --mem=470gb \
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/nocaurban_params_tajima.sh \
        -w 1 -s 1


sbatch --account=mcnew \
        --job-name=tajima_1_nocarural \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1_nocarural.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=48:00:00 \
        --mem=470gb \
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/nocarural_params_tajima.sh \
        -w 1 -s 1

sbatch --account=mcnew \
        --job-name=tajima_1_pyrrurban \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1_pyrrurban.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=48:00:00 \
        --mem=470gb \
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/pyrrurban_params_tajima.sh \
        -w 1 -s 1


sbatch --account=mcnew \
        --job-name=tajima_1_pyrrrural \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1_pyrrrural.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=48:00:00 \
        --mem=470gb \
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/pyrrrural_params_tajima.sh \
        -w 1 -s 1

# 12115213-12115216