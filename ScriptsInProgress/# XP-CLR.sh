# XP-CLR

/xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/nocaurban.${SCAFFOLD}.phased.vcf
/xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/nocaurban.NC_044571.1.phased.vcf

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/xpclr/noca


# Set working directory
VCF_DIR="/xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs"
OUT_DIR="$VCF_DIR/merged_by_species_chrom"

mkdir -p "$OUT_DIR"

# bzip all files
cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs

for f in *.vcf; do
    bgzip -c "$f" > "$f.gz"
    tabix -p vcf "$f.gz"
done

# Get unique chromosomes from filenames
chroms="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

for chrom in $(cat "$chroms"); do
    for species in noca pyrr; do
        # Find all files for this species and chromosome
        files=($(ls $VCF_DIR/${species}*.$chrom.phased.vcf.gz 2>/dev/null))
        
        # Proceed only if there are at least two files to merge
        if [ "${#files[@]}" -ge 1 ]; then
            out_file="$OUT_DIR/${species}.${chrom}.merged.vcf"
            echo "Merging ${#files[@]} files for $species $chrom into $out_file"
            bcftools merge -Oz -o "$out_file" "${files[@]}"
            tabix -p vcf "$out_file"
        fi
    done
done

for chrom in $(cat "$chroms"); do
    for species in noca pyrr; do
        out_file="$OUT_DIR/${species}.${chrom}.merged.vcf"
        mv $out_file $out_file.gz
        zcat $out_file.gz > $out_file
    done
done

# 

xpclr --help


CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"


xpclr -I /xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/merged_by_species_chrom/noca.NC_044600.1.merged.vcf \
-Sa /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocaurban.txt \
-Sb /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocarural.txt \
-O cra.NC_044600 -C NC_044600.1 --size 50000 --step 12500 -V 50

#!/bin/bash

CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r CHROM; do

echo "starting $CHROM"

xpclr -I /xdisk/mcnew/dannyjackson/cardinals/datafiles/mergedvcfs/merged_by_species_chrom/noca.${CHROM}.merged.vcf \
-Sa /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocaurban.txt \
-Sb /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocarural.txt \
-O noca.${CHROM} -C ${CHROM} --size 50000 --step 12500 -V 50

echo "finished $CHROM"

done < "$CHROM"



sbatch --account=mcnew \
--job-name=xpclr_noca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/xpclr_noca%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=72:00:00 \
xpclr.noca.sh


