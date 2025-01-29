Notes on how I ran Fst on my Cardinalis genomes

Index the reference genome
```
module load samtools/1.19.2
samtools faidx /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
```

Make chromosome conversion file
```
echo '1,NC_044571.1' > /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '2,NC_044572.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '3,NC_044573.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '4,NC_044574.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '5,NC_044575.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '6,NC_044576.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '7,NC_044577.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '8,NC_044578.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '9,NC_044579.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '10,NC_044580.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '11,NC_044581.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '12,NC_044582.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '13,NC_044583.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '14,NC_044584.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '15,NC_044585.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '1A,NC_044586.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '17,NC_044587.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '18,NC_044588.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '19,NC_044589.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '20,NC_044590.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '21,NC_044591.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '22,NC_044592.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '23,NC_044593.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '24,NC_044594.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '25,NC_044595.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '26,NC_044596.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '27,NC_044597.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '28,NC_044598.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '29,NC_044599.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo '4A,NC_044600.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
echo 'Z,NC_044601.1' >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
```


Set the working directory to fst so that all relevant slurm output files appear here as well
```
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst
```


Run fst 1 on both northern cardinals and pyrrhuloxia
```
sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh
# 3694147

sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_1.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh
# 3694148
```


Run fst 2 on both northern cardinals and pyrrhuloxia
```
sbatch --account=mcnew \
--job-name=fst_2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_2.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_2.sh -p ~/programs/CardinalisGenomics/noca_params.sh -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt

sbatch --account=mcnew \
--job-name=fst_2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_2.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_2.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
```


Run fst 3 on both northern cardinals and pyrrhuloxia
First with 10,000 bp windows:
```
sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/Genomics-Main/noca_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt

sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/Genomics-Main/pyrr_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
```

And next with 1,000 bp windows:
```
sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/Genomics-Main/noca_params.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt

sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/Genomics-Main/pyrr_params.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/GCF_901933205_chromconversion.txt
```