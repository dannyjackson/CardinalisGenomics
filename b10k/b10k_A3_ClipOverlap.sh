## Clip overlapping read pairs
#!/bin/bash

#SBATCH --job-name=clipoverlap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.clipping.%j

cd /xdisk/mcnew/dannyjackson/cardinals/clipoverlap

module load parallel

clipping () {
echo "clipping" "$@" >> /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/clippingstats.txt 


/home/u15/dannyjackson/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/cardinals/sortedbamfiles_indv/"$@"_sorted_RGadded.bam --out /xdisk/mcnew/dannyjackson/cardinals/clipoverlap/"$@".all.sorted.marked.clipped.bam --stats --params


echo "done " $@ >>/xdisk/mcnew/dannyjackson/cardinals/clipoverlap/clippingstats.txt 

}


export -f clipping 

parallel -j 12 clipping :::: /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch clipoverlap.sh 
Submitted batch job 2162945