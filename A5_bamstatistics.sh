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


#!/bin/bash

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats

while read -r cardinal;
do 
  echo $cardinal
  echo $cardinal >> depthstats.txt

  # average and standard deviation
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$cardinal"_depthstats.txt >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch --account=mcnew \
--job-name=calcdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.calcdepth.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/calcdepth.sh

Submitted batch job 11131337

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


# depth stats aligned to  darwin's finch genome
Sample  Average Stdev
MSB25201	6.97814	51.5322
NOCA003	4.49975	35.926
NOCA004	3.50371	32.1905
NOCA006	5.80631	48.1065
NOCA008	4.85283	41.6011
NOCA012	4.58593	40.8785
NOCA013	5.25671	50.7062
PYRR003	4.08382	22.4905
PYRR004	5.57957	35.043
PYRR006	4.62345	24.0894
PYRR007	4.20723	23.6736
PYRR009	4.5434	27.004
PYRR011	4.83675	31.0408
UWBM100619	6.6082	39.3782
UWBM100620	5.66758	24.0897
UWBM100621	6.64143	30.1046
UWBM103345	5.12688	40.3668
UWBM103346	4.42484	19.0454
UWBM77548	4.97307	15.239
UWBM77718	4.85499	19.9508
UWBM77780	5.73054	28.0262
UWBM77781	4.00269	43.9651
UWBM77856	4.82667	30.0542
UWBM77978	5.66109	28.2475