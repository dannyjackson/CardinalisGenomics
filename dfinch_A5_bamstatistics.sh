# bam statistics 


## Compute statistics on bam files 

#!/bin/bash
module load samtools
while read -r bird; do
echo $bird >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats/"$bird"_depthstats.txt 
samtools depth /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment/"$bird".realigned.bam >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats/"$bird"_depthstats.txt 
done <  /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch --account=mcnew \
--job-name=bamstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.bamstats.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats.sh 
# Submitted batch job 12031591

cd /xdisk/mcnew/dannyjackson/cardinals/bamstats

while read -r cardinal;
do 
  echo $cardinal
  echo $cardinal >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$cardinal".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt



# depth stats aligned to b10k genome
Sample  Average Stdev
MSB25201	7.93225	44.6036
NOCA003	4.94421	24.7466
NOCA004	4.30927	23.3997
NOCA006	6.75042	30.0491
NOCA008	5.75371	26.9686
NOCA012	5.1381	26.5608
NOCA013	5.93655	30.8417
PYRR003	4.91073	26.0414
PYRR004	6.63732	34.9726
PYRR006	5.30687	27.212
PYRR007	5.04058	27.263
PYRR009	5.16903	29.4502
PYRR011	5.47466	32.114
UWBM100619	7.50557	32.8815
UWBM100620	6.55938	23.1015
UWBM100621	7.731	32.225
UWBM103345	5.81635	28.8372
UWBM103346	4.90758	24.5868
UWBM77548	5.67293	20.6727
UWBM77718	5.65916	25.5391
UWBM77780	6.51966	33.7783
UWBM77781	4.35542	39.769
UWBM77856	5.5126	32.4115
UWBM77978	6.47198	29.9459
mean	5.833972083	29.498825
