A1_preprocessing.md

cd /xdisk/mcnew/dannyjackson/cardinals

# A1.3 site allele frequency
species=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
    --job-name=saf \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.saf.%j \
    --nodes=1 \
    --ntasks-per-node=4 \
    --time=5:00:00 \
    ~/programs/CardinalisGenomics/Genomics-Main/A_Preprossesing/A1.3_siteallelefrequency.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -n ${sp}
done

# A2.1 call variants
sbatch --account=mcnew \
--job-name=callvariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.callvariants.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=240:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprossesing/A2.1_callvariants.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -r cardinalis
# 12102489

# A2.2 generate mask
sbatch --account=mcnew \
--job-name=generate_mask \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.generate_mask.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=20:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprossesing/A2.2_generate_mask.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh
# 12102822