# C4_tajima.md

# run once per population to generate all files
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )

# Iterate over each combination
for sp in "${species[@]}"; do
sbatch --account=mcnew \
        --job-name=tajima_500000_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_500000_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=10:00:00 \
        --mem=50gb \
        ~/programs/CardinalisGenomics/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/${sp}_params_tajima.sh \
        -w 500000 -s 500000
done


# Define species, environments, and window sizes
declare -A time_limits=( [50000]=1:00:00 [10000]=1:00:00 [5000]=1:00:00 [1000]=1:00:00 [1]=10:00:00 )
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
window_sizes=( 50000 25000 10000 5000 1000 1)

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}

        sbatch --account=mcnew \
               --job-name=tajima_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.tajima_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=$time_limit \
               --mem=50gb \
               ~/programs/CardinalisGenomics/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_tajima.sh \
               -w $win -s 1
    done
done





# troubleshoot snps sbatch settings

sbatch --account=mcnew \
        --job-name=tajima_1_nocaurban \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1000_nocaurban.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=5:00:00 \
        --mem=100gb \
        ~/programs/CardinalisGenomics/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/CardinalisGenomics/nocaurban_params_tajima.sh \
        -w 1000 -s 1000 -m Tajima

# 12129439

sbatch --account=mcnew \
        --job-name=tajima_1_nocarural \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_1_nocarural.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=48:00:00 \
        --mem=470gb \
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima_plot.sh \
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
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima_plot.sh \
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
        ~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima_plot.sh \
        -p ~/programs/CardinalisGenomics/pyrrrural_params_tajima.sh \
        -w 1 -s 1

# 12117077-12117080