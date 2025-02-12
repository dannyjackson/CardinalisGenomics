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
declare -A time_limits=( [500000]=1:00:00 [10000]=1:00:00 [1000]=1:00:00 [1]=48:00:00 )
declare -A mem_limits=( [500000]="10gb" [10000]="10gb" [1000]="10gb" [1]="100gb" )
species=( "noca" "pyrr" )
window_sizes=( 10000 1000 1 )

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


