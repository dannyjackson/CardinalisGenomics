# C3_dxy.md
# testing

cd ~/programs/CardinalisGenomics/
chmod -x *
git pull
cd ~/programs/Genomics-Main/
chmod -x *
git pull
chmod +x dxy_1.sh
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy

~/programs/CardinalisGenomics/Genomics-Main/dxy_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt



# script
sbatch --account=mcnew \
--job-name=dxy_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/dxy_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# 3696632

sbatch --account=mcnew \
--job-name=dxy_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/Genomics-Main/dxy_1.sh -p ~/programs/CardinalisGenomics/pyrr_params.sh -w 10000 -s 10000 -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
# 3696633

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

