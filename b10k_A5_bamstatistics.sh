# bam statistics 


## Compute statistics on bam files 
# some alignment stats 

#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/cardinals/bamstats

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/"$@".realigned.bam > "$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/"$@".realigned.bam -o "$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch bamstats.sh 
Submitted batch job 10707464


# compute average and sd depth per individual # sorted marked bam files
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/cardinals/bamstats


while read -r cardinal;
do 
  echo $cardinal
  echo $cardinal >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$cardinal".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt


sbatch avgsd_bamstats.sh 
Submitted batch job 2184368




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


# python script for summarizing depth stats from ANGSD
python
import numpy as np
import pandas as pd 

df = pd.read_csv('stats.NC_044571.depthSample', sep='\t', header=None)
cols = df.shape[1]
multiplier = list(range(1, cols+1))

df_adj = df.mul(multiplier, axis = 1)

counts = df.sum(axis = 1)
totaldepth = df_adj.sum(axis = 1)
avgdepth = totaldepth / counts

df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = pd.concat([df_names, avgdepth], axis=1)

df_final = df_final.rename({0: 'avgdepth'}, axis=1)

df_final.groupby(['species', 'treatment']).mean()

df_final.sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'CRA'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'FOR'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'PAR'].sort_values(by=['avgdepth']).round(decimals=2)

print(df_final.groupby(['species', 'treatment']).size())

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')



