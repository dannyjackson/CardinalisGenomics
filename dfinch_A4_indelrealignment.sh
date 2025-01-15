## Indel realignment
## Index bams

for i in `cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=clip_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.clip_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/index_clipped.sh $i
done

#!/bin/bash
IND=$1
module load samtools
samtools index /xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/$IND.all.sorted.marked.clipped.bam 
echo "done " ${IND} >>/xdisk/mcnew/dannyjackson/cardinals_dfinch/index_clippedstats.txt 


for i in `cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=indelmap_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.indelmap_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/indelmap.sh $i
done

#!/bin/bash
IND=$1

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
-I /xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/${IND}.all.sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelmaps/${IND}.intervals

# 11126461 - 11126484




for i in `cat /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=indelrealignment_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.indelrealignment_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment.sh $i
done

#!/bin/bash
IND=$1

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
-R /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
--consensusDeterminationModel USE_READS \
-I /xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/${IND}.all.sorted.marked.clipped.bam \
--targetIntervals /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelmaps/${IND}.intervals \
-o /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment/${IND}.realigned.bam


# 3467823 - 3467870