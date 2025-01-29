# C3_dxy.md

sbatch --account=mcnew \
--job-name=fst_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.dxy_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=5:00:00 \
~/programs/CardinalisGenomics/dxy_1.sh -p ~/programs/CardinalisGenomics/noca_params.sh -w 10000 -s 10000