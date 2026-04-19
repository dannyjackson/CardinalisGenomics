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
window_sizes=( 50000 5000 )


# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        step=$(expr $win / 4)

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


window_sizes=( 1 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}

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
               -w $win -s $win 
    done
done




# Make violin plots
library(ggplot2)
library(dplyr)

# Read the files
file1 <- read.table("nocaurban_nocarural/50000/nocaurban_nocarural.50000.fst", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table("pyrrurban_pyrrrural/50000/pyrrurban_pyrrrural.50000.fst", header = TRUE, stringsAsFactors = FALSE)

file1 <- read.table("nocaurban_nocarural/1/nocaurban_nocarural.1.fst", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table("pyrrurban_pyrrrural/1/pyrrurban_pyrrrural.1.fst", header = TRUE, stringsAsFactors = FALSE)

file1$dataset <- "file1"
file2$dataset <- "file2"

combined <- bind_rows(file1, file2)

combined$dataset <- factor(combined$dataset, 
                           levels = c("file1", "file2"), 
                           labels = c("noca", "pyrr"))

wilcox_test <- wilcox.test(fst ~ dataset, data = combined)
print(wilcox_test)

# snp data
data:  fst by dataset
W = 3.7979e+11, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0

# 50000 data
data:  fst by dataset
W = 2613407156, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0

# Create and save the plot

p <- ggplot(combined, aes(x = dataset, y = fst, fill = dataset)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  theme_minimal() +
  labs(x = "", y = "FST value", title = "FST Comparison: noca vs pyrr") +
  theme(legend.position = "none")

# Save to file
ggsave("fst_violin_noca_pyrr.1.png", plot = p, width = 6, height = 4, dpi = 300)
# ggsave("fst_violin_noca_pyrr.50000.png", plot = p, width = 6, height = 4, dpi = 300)


# compute FST comparison between taxa
# pyrr urban vs noca urban
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/pyrrurban.saf.idx  /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/nocaurban.saf.idx > /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrurban_nocaurban.ml

POP1="pyrrurban"
POP2="nocaurban"
OUTDIR="/xdisk/mcnew/dannyjackson/cardinals/"
ANGSD="~/programs/angsd/"
# Compute FST index
FST_INDEX="${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx"

~/programs/angsd/misc/realSFS fst index "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" \
    -sfs /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrurban_nocaurban.ml -fstout "${OUTDIR}/analyses/fst/${POP1}_${POP2}"


# Compute global FST estimate
GLOBAL_FST_FILE="${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt"
if [ -f "$GLOBAL_FST_FILE" ]; then
    echo "Global FST estimate already exists, skipping computation."
else
    echo "Computing global FST estimate..."
    echo -e "FST.Unweight\tFST.Weight" > "$GLOBAL_FST_FILE"
    ~/programs/angsd/misc/realSFS fst stats "$FST_INDEX" >> "$GLOBAL_FST_FILE"
fi
## -> FST.Unweight[nObs:895129]:0.414326 Fst.Weight:0.638210

# pyrr urban vs noca rural
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/pyrrurban.saf.idx  /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/nocarural.saf.idx > /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrurban_nocarural.ml

POP1="pyrrurban"
POP2="nocarural"
OUTDIR="/xdisk/mcnew/dannyjackson/cardinals/"
ANGSD="~/programs/angsd/"
# Compute FST index
FST_INDEX="${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx"

~/programs/angsd/misc/realSFS fst index "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" \
    -sfs /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrurban_nocarural.ml -fstout "${OUTDIR}/analyses/fst/${POP1}_${POP2}"

# Compute global FST estimate
GLOBAL_FST_FILE="${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt"
if [ -f "$GLOBAL_FST_FILE" ]; then
    echo "Global FST estimate already exists, skipping computation."
else
    echo "Computing global FST estimate..."
    echo -e "FST.Unweight\tFST.Weight" > "$GLOBAL_FST_FILE"
    ~/programs/angsd/misc/realSFS fst stats "$FST_INDEX" >> "$GLOBAL_FST_FILE"
fi
## -> FST.Unweight[nObs:895133]:0.414907 Fst.Weight:0.637908


# pyrr rural vs noca urban

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/pyrrrural.saf.idx  /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/nocaurban.saf.idx > /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrrural_nocaurban.ml

POP1="pyrrrural"
POP2="nocaurban"
OUTDIR="/xdisk/mcnew/dannyjackson/cardinals/"
ANGSD="~/programs/angsd/"
# Compute FST index
FST_INDEX="${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx"

~/programs/angsd/misc/realSFS fst index "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" \
    -sfs /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrrural_nocaurban.ml -fstout "${OUTDIR}/analyses/fst/${POP1}_${POP2}"

# Compute global FST estimate
GLOBAL_FST_FILE="${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt"
if [ -f "$GLOBAL_FST_FILE" ]; then
    echo "Global FST estimate already exists, skipping computation."
else
    echo "Computing global FST estimate..."
    echo -e "FST.Unweight\tFST.Weight" > "$GLOBAL_FST_FILE"
    ~/programs/angsd/misc/realSFS fst stats "$FST_INDEX" >> "$GLOBAL_FST_FILE"
fi
## -> FST.Unweight[nObs:895135]:0.415734 Fst.Weight:0.637540

# pyrr rural vs noca rural

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/pyrrrural.saf.idx  /xdisk/mcnew/dannyjackson/cardinals/datafiles/safs/nocarural.saf.idx > /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrrural_nocarural.ml

POP1="pyrrrural"
POP2="nocarural"
OUTDIR="/xdisk/mcnew/dannyjackson/cardinals/"
ANGSD="~/programs/angsd/"
# Compute FST index
FST_INDEX="${OUTDIR}/analyses/fst/${POP1}_${POP2}.fst.idx"

~/programs/angsd/misc/realSFS fst index "${OUTDIR}/datafiles/safs/${POP1}.saf.idx" "${OUTDIR}/datafiles/safs/${POP2}.saf.idx" \
    -sfs /xdisk/mcnew/dannyjackson/cardinals/datafiles/pyrrrural_nocarural.ml -fstout "${OUTDIR}/analyses/fst/${POP1}_${POP2}"

# Compute global FST estimate
GLOBAL_FST_FILE="${OUTDIR}/analyses/fst/${POP1}_${POP2}.globalFST.txt"
if [ -f "$GLOBAL_FST_FILE" ]; then
    echo "Global FST estimate already exists, skipping computation."
else
    echo "Computing global FST estimate..."
    echo -e "FST.Unweight\tFST.Weight" > "$GLOBAL_FST_FILE"
    ~/programs/angsd/misc/realSFS fst stats "$FST_INDEX" >> "$GLOBAL_FST_FILE"
fi

##        -> FST.Unweight[nObs:895049]:0.415577 Fst.Weight:0.636954

