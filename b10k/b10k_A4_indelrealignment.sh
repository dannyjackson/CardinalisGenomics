## Indel realignment

module load samtools 


cd /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1

cp GCA_013397215.1_ASM1339721v1_genomic.fna GCA_013397215.1_ASM1339721v1_genomic.fa

samtools faidx GCA_013397215.1_ASM1339721v1_genomic.fa

samtools dict GCA_013397215.1_ASM1339721v1_genomic.fa -o GCA_013397215.1_ASM1339721v1_genomic.dict


indexbams () {
samtools index /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$@".all.sorted.marked.clipped.bam
}

export -f indexbams 

parallel -j 12 indexbams :::: /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt


#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=1:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools
module load parallel

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/MSB25201.all.sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/cardinals/indelmaps/MSB25201.test.intervals

sbatch indelmap_test.sh 
Submitted batch job 2163042


#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

while read -r bird;
do

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$bird".all.sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/cardinals/indelmaps/"$bird".intervals

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch indelmap.sh 
Submitted batch job 2163045



# didn't finish in time, rerun with last two
#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=3:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

while read -r bird;
do

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$bird".all.sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/cardinals/indelmaps/"$bird".intervals

done < tmp.sampleids.txt

sbatch indelmap_tmp.sh 
Submitted batch job 10691104

#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$bird".all.sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/cardinals/indelmaps/"$bird".intervals

# test
#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=2:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j


apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
--consensusDeterminationModel USE_READS \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/MSB25201.all.sorted.marked.clipped.bam \
--targetIntervals /xdisk/mcnew/dannyjackson/cardinals/indelmaps/MSB25201.intervals \
-o /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/MSB25201.realigned.bam

Submitted batch job 2163055



#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j

while read -r bird;
do

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
-R /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fa \
--consensusDeterminationModel USE_READS \
-I /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$bird".all.sorted.marked.clipped.bam \
--targetIntervals /xdisk/mcnew/dannyjackson/cardinals/indelmaps/"$bird".intervals \
-o /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/"$bird".realigned.bam

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch indelrealign.sh 
Submitted batch job 10703001
