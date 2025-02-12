# gene visualization

find /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/* | awk -F"[./]" '{ OFS="." ; print $(NF-5), $(NF-4), $(NF-3) }' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/genelistnames.txt


Rscript ~/ /xdisk/mcnew/dannyjackson/cardinals/referencelists/genelistnames.txt

while read -r file;
do

    sbatch --account=mcnew \
            --job-name=visualization_${file} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.visualization_${file}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=1:00:00
            ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_visualization.sh \
            -p ~/programs/CardinalisGenomics/params_base.sh \
            -f $file 

done < filenames.txt