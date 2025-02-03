# C4_tajima.md


# compute 500,000bp windows and snp analysis 
# northern cardinals
# 12089968
sbatch --account=mcnew \
--job-name=tajima_50k_nocaurban \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.tajima_50k_nocaurban.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/tajima/tajima.sh -p ~/programs/CardinalisGenomics/nocaurban_params_tajima.sh -w 500000 -s 500000 

# -c /xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt