# Assess DNA damage

mkdir -p /xdisk/mcnew/dannyjackson/cardinals/analyses/mapdamage
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/mapdamage

module load python

pip install mapdamage

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/mapdamage
module load python
module load micromamba
micromamba activate r_elgato

BAMDIR=/xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/
REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna # path to reference genome


BAM=${BAMDIR}/MSB25201.realigned.bam
/usr/bin/time -v mapDamage \
    -i "$BAM" \
    -r "$REF" \
    -d "$OUTDIR" \
    --merge-reference-sequences \
    -n 0.01 

mapDamage -i ${BAMDIR}/NOCA_003.realigned.bam -r $REF -n 0.1 --merge-reference-sequences
mapDamage -i ${BAMDIR}/NOCA_004.realigned.bam -r $REF -n 0.1 --merge-reference-sequences
mapDamage -i ${BAMDIR}/UWBM100619.realigned.bam -r $REF -n 0.1 --merge-reference-sequences


#!/bin/bash
#SBATCH --job-name=mapdamage
#SBATCH --output=logs/mapdamage_%A_%a.out
#SBATCH --error=logs/mapdamage_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=50G
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=1
#SBATCH --array=0-24   # adjust if number of BAMs changes

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/mapdamage

module load python
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate r_elgato

BAMDIR=/xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment
REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna

bam_list.txt  results_MSB25201.realigned  results_NOCA003.realigned  results_NOCA004.realigned  results_NOCA006.realigned  results_NOCA008.realigned  run_mapdamage_array.sh

ls ${BAMDIR}/*.realigned.bam | sort > bam_list.txt

while read -r BAM; do
    echo "Processing $BAM"

    BASE=$(basename "$BAM" .bam)
    OUTDIR="results_${BASE}"

    mkdir -p "$OUTDIR"

    mapDamage \
        -i "$BAM" \
        -r "$REF" \
        -d "$OUTDIR" \
        -n 0.01 \
        --merge-reference-sequences

    echo "Finished $BAM"
    echo "--------------------------------------"
done < bam_list.txt

BAMDIR=/xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment
REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna

# Build list of BAMs
mapfile -t BAMFILES < <(ls ${BAMDIR}/*.realigned.bam | sort)

BAM=${BAMFILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing $BAM"

mapDamage \
    -i "$BAM" \
    -r "$REF" \
    -n 0.01 \
    --merge-reference-sequences


# Combine data output
echo -e "sample\tCtoT_5p\tGtoA_3p" > damage_summary.tsv

for d in results_*.realigned; do
    sample=$(basename "$d" | sed 's/results_//')

    ctot=$(awk 'NR>1 && NR<=4 {sum+=$2; n++} END {if(n>0) print sum/n; else print "NA"}' "$d/5pCtoT_freq.txt")
    gtoa=$(awk 'NR>1 && NR<=4 {sum+=$2; n++} END {if(n>0) print sum/n; else print "NA"}' "$d/3pGtoA_freq.txt")

    echo -e "${sample}\t${ctot}\t${gtoa}" >> damage_summary.tsv
done

sed -i 's/\.realigned//g' damage_summary.tsv

cat <<EOF > metadata.tsv
sample	type    species
MSB25201	Museum  PYRR
NOCA003	Field   NOCA
NOCA004	Field   NOCA
NOCA006	Field   NOCA
NOCA008	Field   NOCA
NOCA012	Field   NOCA    
NOCA013	Field   NOCA
PYRR003	Field   PYRR
PYRR004	Field   PYRR
PYRR006	Field   PYRR
PYRR007	Field   PYRR
PYRR009	Field   PYRR
PYRR011	Field   PYRR
UWBM100619	Museum  NOCA
UWBM100620	Museum  NOCA
UWBM100621	Museum  NOCA
UWBM103345	Museum  NOCA
UWBM103346	Museum  PYRR
UWBM77548	Museum  PYRR
UWBM77718	Museum  PYRR
UWBM77780	Museum  PYRR
UWBM77781	Museum  PYRR
UWBM77856	Museum  NOCA
UWBM77978	Museum  NOCA
EOF


R
library(ggplot2)
library(dplyr)

damage <- read.table("damage_summary.tsv", header=TRUE)
meta   <- read.table("metadata.tsv", header=TRUE)

df_full <- merge(damage, meta, by="sample")
df_full$damage_combined <- (df_full$CtoT_5p + df_full$GtoA_3p)/2


df_full %>%
  group_by(species) %>%
  summarise(
    field_mean  = mean(damage_combined[type == "Field"]),
    museum_mean = mean(damage_combined[type == "Museum"]),
    p_value     = wilcox.test(damage_combined ~ type)$p.value
  )

p <- ggplot(df_full, aes(x = type, y = GtoA_3p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = sample), width = 0.1, size = 2) +
  facet_wrap(~ species) +
  theme_classic()

ggsave("GtoA_3p_boxplot_by_species.pdf", plot = p, width = 7, height = 4)

  species field_mean museum_mean p_value
  <chr>        <dbl>       <dbl>   <dbl>
1 NOCA        0.0243      0.0251   0.240
2 PYRR        0.0243      0.0255   0.180

(Wilcox test: Northern cardinals field mean = 0.024, museum mean = 0.025, p-value = 0.240; Pyrrhuloxia field mean = 0.024, museum mean = 0.026, p-value = 0.180)