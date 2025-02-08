
# snps 
# /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_10000.snps.outlier.csv
# /xdisk/mcnew/dannyjackson/cardinals/analyses/dxy/pyrrurban_pyrrrural/pyrrurban_pyrrrural.dxy_10000.snps.outlier.csv

# windowed
IN_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrrurban_pyrrrural.fst_10000.snps.outlier.csv
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
OUT_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/pyrr.fst.10kb.bed
GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
GENES_FILE=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/pyrr.fst.10kb.genes.txt
GENENAMES=/xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/pyrr.fst.10kb.genenames.txt


# sed -i 's/\"//g' ${IN_FILE}

module load python/3.11/3.11.4

python ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py ${IN_FILE} 10000 ${OUT_FILE} ${CHROM_CONVERSION}

# has to be done on elgato for now
module load gnu8/8.3.0
module load bedtools2/2.29.2

bedtools intersect -a ${GFF} -b ${OUT_FILE} -wa > ${GENES_FILE}

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${GENES_FILE} | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}






# snps
awk 'BEGIN {FS = ","} {$1=""}1' ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv.tmp

mv ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv.tmp ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv

sed -i 's/\"//g' ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv

tail -n +2 ${OUTDIR}/analyses/fst/pyrr.outlierfst.csv  > ${OUTDIR}/analyses/fst/pyrr.outlierfst.headless.csv


awk -F',' 'NR>1 {print "NC_0"$1".1" "\t" $2-1 "\t" $2}' ${OUTDIR}/analyses/fst/pyrr.outlierfst.headless.csv > ${OUTDIR}/analyses/fst/pyrr.outlierfst.bed

bedtools intersect -a /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff -b ${OUTDIR}/analyses/fst/pyrr.outlierfst.bed -wa > ${OUTDIR}/analyses/genelist/relevantgenes_snps_top.95.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' ${OUTDIR}/analyses/genelist/relevantgenes_snps_top.95.txt | sed 's/ID\=gene\-//g' | sort -u > ${OUTDIR}/analyses/genelist/relevantgenenames_snps_top.95.txt
