dfinch_B2_FST.sh


# update git code
cd ~/programs/Genomics-Main/
chmod -x fst_1.sh
git pull

# test new code
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/fst
chmod +x ~/programs/Genomics-Main/fst_1.sh


# index reference genome
module load samtools/1.19.2
samtools faidx /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna


sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fst_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=40:00:00 \
~/programs/Genomics-Main/fst_1.sh -p ~/programs/Genomics-Main/fst_params.sh

# Submitted batch job 3691902