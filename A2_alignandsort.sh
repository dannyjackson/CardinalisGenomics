dfinch_A2_alignandsort.sh

# align 

for i in `cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=align_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=output.align_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=12 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/align.sh $i
done

# 11123375 - 11123398

#!/bin/bash
IND=$1

module load bwa
module load samtools
bwa mem -t 12 /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    /xdisk/mcnew/dannyjackson/cardinals/trimmed_fastas/${IND}_trimmed_1P.fq.gz \
    /xdisk/mcnew/dannyjackson/cardinals/trimmed_fastas/${IND}_trimmed_2P.fq.gz | \
    samtools view -b -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/bamfiles/${IND}.bam -S 






# sort with a slurm array

for i in `cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=sort_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=output.sort_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=10:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/sort.sh $i
done

# 11123978 - 11124001

# this one is still working so leave it be and rerun sort.sh on it solo once it's done # 11123390 UWBM100621

sbatch --account=mcnew \
	--job-name=sort_UWBM100621 \
    --partition=standard \
	--mail-type=ALL \
	--output=output.sort_UWBM100621.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=10:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/sort.sh UWBM100621

# 3419887


#!/bin/bashcd
IND=$1

module load bowtie2
module load picard
module load samtools
module load parallel

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles

# Specify the path to the config file
mkdir /xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}

samtools sort -T temp -@ 12 \
-o /xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted.bam \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/bamfiles/${IND}.bam

picard AddOrReplaceReadGroups \
I=/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted.bam \
O=/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${IND}

picard MarkDuplicates \
I=/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted_RGadded.bam \
O=/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam \
M=/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}.duplicate.metrics.txt

samtools index \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam

# some alignment stats 

while read -r SAMPLE;
do

echo "$SAMPLE"

samtools flagstat /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

