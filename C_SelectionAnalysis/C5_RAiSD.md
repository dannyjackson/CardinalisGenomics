# RAiSD.md 

# merge vcfs

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocaurban.txt

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocarural.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrurban.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrrural.txt




sbatch --account=mcnew \
    --job-name=raisd.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.raisd.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=10:00:00 \
    ~/programs/raisd_preprocessing.sh 



for GROUP in nocaurban nocarural pyrrurban pyrrrural; do
    sbatch --account=mcnew \
        --job-name=raisd.test \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.raisd.test.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=10:00:00 \
        ~/programs/CardinalisGenomics/Genomics-Main/raisd/raisd.sh -p ~/programs/CardinalisGenomics/Genomics-Main/raisd/params_raisd.sh -n ${GROUP} -w 1000
done 
# 12114965-12114968