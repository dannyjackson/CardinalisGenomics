A1_preprocessing.md

cd /xdisk/mcnew/dannyjackson/cardinals
# A0.5 Bam statistics
sbatch --account=mcnew \
--job-name=bamstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.bamstats.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=10:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A0.5_bamstatistics.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -r cardinalis

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
    ~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -n ${sp}
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
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.1_callvariants.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -r cardinalis
# 12102489

# A2.2 filter vcf
sbatch --account=mcnew \
--job-name=filtervcf \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.filtervcf.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=240:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.2_filtervcf.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh 
# 12107182

# A2.3 generate mask
sbatch --account=mcnew \
--job-name=generate_mask \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.generate_mask.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=20:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.3_generate_mask.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh
# 12102865


# A2.4 individual mask vcf 
# test
sbatch --account=mcnew \
--job-name=ind_mask_vcf \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.ind_mask_vcf.%j \
--nodes=1 \
--ntasks-per-node=8 \
--time=48:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -i NOCA003
# 12107151

# run in a slurm array
for i in `cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_subset.txt`;
	do echo $i
	IND=$i
    sbatch --account=mcnew \
    --job-name=mask_vcf_${i} \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.mask_vcf_${i}.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -i $i
done


# phasing

# test on one ind
sbatch --account=mcnew \
--job-name=phase_NOCA003 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.phase_NOCA003.%j \
--nodes=1 \
--ntasks-per-node=8 \
--time=48:00:00 \
~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.5_phasing.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -i NOCA003
# 12113101

# run in a slurm array
for i in `cat /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_subset.txt`;
	do echo $i
	IND=$i
    sbatch --account=mcnew \
    --job-name=phase_${i} \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.phase_${i}.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/A2.5_phasing.sh -p ~/programs/CardinalisGenomics/params_preprocessing.sh -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/indelrealignment/ -i $i
done

# 12113266 - 12113289