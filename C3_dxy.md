# C3_dxy.md

# script
sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/C3_dxy_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -s 10000

sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/C3_dxy_1.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -w 10000 -s 10000


# testing

cd ~/programs/CardinalisGenomics/
chmod -x *
git pull
cd ~/programs/Genomics-Main/
chmod -x *
git pull
chmod +x dxy_1.sh
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy

~/programs/Genomics-Main/dxy_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -s 10000

