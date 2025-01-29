## Clip overlapping read pairs

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
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/clip.sh $i
done


# 11125990 - 11126013

#!/bin/bash
IND=$1

echo "clipping" "${IND}" >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/clippingstats.txt 

/home/u15/dannyjackson/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/cardinals_dfinch/sortedbamfiles/${IND}/${IND}_sorted_RGadded_dupmarked.bam --out /xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/${IND}.all.sorted.marked.clipped.bam --stats --params

echo "done " ${IND} >>/xdisk/mcnew/dannyjackson/cardinals_dfinch/clipoverlap/clippingstats.txt 

sbatch --account=mcnew \
	--job-name=clip_UWBM100621 \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.clip_UWBM100621.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/cardinals_dfinch/clip.sh UWBM100621

# Submitted batch job 11126328




