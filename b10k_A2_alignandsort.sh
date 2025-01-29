# align and sort

#!/bin/bash

#SBATCH --job-name=alignsort
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j


cd /xdisk/mcnew/dannyjackson/cardinals/bamfiles


module load bowtie2
module load picard
module load samtools
module load parallel

align () {
    file=$1  # Get the file name from the argument

    echo "aligning $file" >> bwa_alignment_log.txt

    bwa mem -t 12 /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna \
    /xdisk/mcnew/dannyjackson/cardinals/trimmed_fastas/"$file"_trimmed_1P.fq.gz \
    /xdisk/mcnew/dannyjackson/cardinals/trimmed_fastas/"$file"_trimmed_2P.fq.gz | \
    samtools view -b -o /xdisk/mcnew/dannyjackson/cardinals/bamfiles/"$file".bam -S 

    echo "sam file piped into samtools view to convert to .bam for $file" >> bwa_alignment_log.txt
}


export -f align

parallel align ::: $(cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt)

cd /xdisk/mcnew/dannyjackson/cardinals/scripts/output

sbatch ../2_align.sh 
Submitted batch job 2152309



# sort 


cd /xdisk/mcnew/dannyjackson/cardinals/
mkdir sortedbamfiles





# I'm going to try it with a slurm array
#!/bin/bash 

#SBATCH --job-name=sort_array
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --array=1-24

module load bowtie2
module load picard
module load samtools
module load parallel

cd /xdisk/mcnew/dannyjackson/cardinals/sorted_array

# Specify the path to the config file
config=/xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
SAMPLE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $1}' $config)



samtools sort -T temp -@ 12 \
-o /xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted.bam \
/xdisk/mcnew/dannyjackson/cardinals/bamfiles/"$SAMPLE".bam

picard AddOrReplaceReadGroups \
I=/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$SAMPLE"

picard MarkDuplicates \
I=/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted_RGadded.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
M=/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE".duplicate.metrics.txt

samtools index \
/xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted_RGadded_dupmarked.bam

rm /xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted.bam
rm /xdisk/mcnew/dannyjackson/cardinals/sorted_array/"$SAMPLE"_sorted_RGadded.bam

sbatch ../3_sort.sh 
Submitted batch job 10662364

# revised with new independent directory
sbatch sort_array.sh 
Submitted batch job 10669037


#!/bin/bash 

#SBATCH --job-name=sort
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j


module load bowtie2
module load picard
module load samtools
module load parallel

cd /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv

while read -r SAMPLE;
do
samtools sort -T temp -@ 12 \
-o /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted.bam \
/xdisk/mcnew/dannyjackson/cardinals/bamfiles/"$SAMPLE".bam

picard AddOrReplaceReadGroups \
I=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$SAMPLE"

picard MarkDuplicates \
I=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
M=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE".duplicate.metrics.txt

samtools index \
/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded_dupmarked.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt


sbatch sort_indv.sh 
Submitted batch job 2153222



# redo RG and dup marking

# For some reason, this is impossible to run in slurm. Run interactively.
# trying again, 10669394

interactive -a mcnew -t 24:00:00

module load bowtie2
module load picard
module load samtools
module load parallel

cd /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv

while read -r SAMPLE;
do
picard AddOrReplaceReadGroups \
I=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE

picard MarkDuplicates \
I=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam \
O=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
M=/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE".duplicate.metrics.txt

samtools index \
/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded_dupmarked.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt


# I can't mark duplicates...

while read -r SAMPLE;
do

picard MarkDuplicates \
-I /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam \
-O /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
-M /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE".duplicate.metrics.txt

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt


wget https://github.com/broadinstitute/picard/archive/refs/tags/3.3.0.tar.gz

tar -xf 3.3.0.tar.gz

./gradlew shadowJar

java -jar picard.jar MarkDuplicates \
-I /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/MSB25201_sorted_RGadded.bam \
-O /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/MSB25201_sorted_RGadded_dupmarked.bam \
-M /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/MSB25201.duplicate.metrics.txt


# I can't mark duplicates... giving up and moving forward without it while also looking into it separately...

while read -r SAMPLE;
do

samtools index \
/xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt



# some alignment stats 

while read -r SAMPLE;
do

echo "$SAMPLE"

samtools flagstat /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$SAMPLE"_sorted_RGadded.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

