Notes on how I analyzed the *Cardinalis* genomes using F<sub>ST</sub>

First, I copied the following scripts from my [Genomics-Main](https://github.com/dannyjackson/Genomics-Main) repository into my CardinalisGenomics repository: 
 - fst_params.sh
 - fst_1.sh
 - fst_2.sh
 - fst_3.sh
 - fst_4.sh
 - fst_window.r
 - fst_snps.r

Most of these remained unedited, but I copied them to this repository so that the scripts used for the analyses in this paper can be archived upon publication.

I edited the fst_params.sh script to reflect my desired parameters for each run (now saved as noca_params.sh and pyrr_params.sh). I then performed the following steps:

Index the reference genome:
```
module load samtools/1.19.2
samtools faidx /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
```


Set the working directory to fst so that all relevant slurm output files appear here as well:
```
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst
```


Run fst 1 on both northern cardinals and pyrrhuloxia:
```
sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/fst_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh
# 3696285

sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/fst_1.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh
# 3696282
# 3696445
```


Run fst 2 on both northern cardinals and pyrrhuloxia
```
sbatch --account=mcnew \
--job-name=fst_2_noca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_2.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/fst_2.sh -p ~/programs/CardinalisGenomics/noca_params.sh -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# Submitted batch job 3696587

sbatch --account=mcnew \
--job-name=fst_2_pyrr \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_2.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/fst_2.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# 3696593
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
~/programs/Genomics-Main/fst_3.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# Submitted batch job 3696588

sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# Submitted batch job 3696594
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
~/programs/Genomics-Main/fst_3.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# Submitted batch job 3696595

sbatch --account=mcnew \
--job-name=fst_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_3.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/fst_3.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -w 1000 -s 1000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# Submitted batch job 3696596
```