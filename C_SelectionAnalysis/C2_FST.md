Notes on how I analyzed the *Cardinalis* genomes using F<sub>ST</sub>

First, I copied the following scripts from my [Genomics-Main](https://github.com/dannyjackson/Genomics-Main) repository into my CardinalisGenomics repository: 
 - fst_params.sh
 - fst_1.sh
 - fst_2.sh
 - fst_3.sh
 - fst_4.sh
 - fst_window.r
 - fst_snps.r

Most of these remained unedited, but I copied them to this repository so that the scripts used for the analyses in this paper can be archived upon publication.

I edited the fst_params.sh script to reflect my desired parameters for each run (now saved as noca_params.sh and pyrr_params.sh). I then performed the following steps:

Index the reference genome:
```
module load samtools/1.19.2
samtools faidx /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
```


Set the working directory to fst so that all relevant slurm output files appear here as well:
```
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst
```
# Run once per species to generate all files

# Define species, environments, and window sizes
species=( "noca" "pyrr" )

# Iterate over each combination

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=fst_500000_${sp} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.fst_500000_${sp}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/CardinalisGenomics/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
            -p ~/programs/CardinalisGenomics/${sp}_params_fst.sh \
            -w 500000 -s 500000
done






# Iterate over several window sizes

# Define species, environments, and window sizes
species=( "noca" "pyrr" )
window_sizes=( 25000 5000 )


# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        step=$(expr $win / 2)

        sbatch --account=mcnew \
               --job-name=fst_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.fst_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=1:00:00 \
               --mem=50gb \
               ~/programs/CardinalisGenomics/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
               -p ~/programs/CardinalisGenomics/${sp}_params_fst.sh \
               -w $win -s $step 
    done
done

