# filter snpeff vcf to just sites where urban genomes have no variation

# noca

# noca

#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural
# Query the genotypes for the urban samples and filter for sites where all urban genotypes are the same

bcftools view -HH -s NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz | awk '{if (($10 ~ /^0\/0/ || $10 ~ /^1\/1/) && ($11 ~ /^0\/0/ || $11 ~ /^1\/1/) && ($12 ~ /^0\/0/ || $12 ~ /^1\/1/) && ($13 ~ /^0\/0/ || $13 ~ /^1\/1/) && ($14 ~ /^0\/0/ || $14 ~ /^1\/1/) && ($15 ~ /^0\/0/ || $15 ~ /^1\/1/)) print $0}' > noca_urban_same.txt



sbatch --account=mcnew \
            --job-name=noca_urban_same \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.noca_urban_same$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/noca_urban_same.sh 

#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural


bcftools view -HH -s UWBM100619,UWBM100620,UWBM100621,UWBM103345,UWBM77856,UWBM77978 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz | \
awk '{  
    # Store genotypes in an array, ignoring missing ones  
    delete gt; count = 0;  
    for (i = 10; i <= 14; i++) {  
        if ($i != "./.") {  
            gt[$i]++; count++;  
        }  
    }  

    # If more than one unique non-missing genotype exists, print line  
    if (length(gt) > 1) print $0;  
}' > noca_rural_diff.txt


sbatch --account=mcnew \
            --job-name=noca_rural_diff \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.noca_rural_diff$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/noca_rural_diff.sh 

# Extract the intersection of both urban same and rural differing sites
grep -F -f noca_urban_same.txt noca_rural_diff.txt > nocaurbansame_ruraldiff.txt

awk -v OFS="\t" '{print $1,$2}' noca_urban_same.txt > noca_urban_same.sites.txt
awk -v OFS="\t" '{print $1,$2}' noca_rural_diff.txt > noca_rural_diff.sites.txt

sort noca_urban_same.sites.txt > noca_urban_same_sorted.txt
sort noca_rural_diff.sites.txt > noca_rural_diff_sorted.txt
comm -12 noca_urban_same_sorted.txt noca_rural_diff_sorted.txt > nocaurbansame_ruraldiff.txt


#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural

# Use the combined sites list to extract those variants from the original VCF
bcftools view -R nocaurbansame_ruraldiff.txt /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz -o noca.urbanfixedsites.vcf

sbatch --account=mcnew \
            --job-name=noca_vcf \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.noca_vcf$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/nocaurban_nocarural/filterurbansites.sh 


# see which genes are in these files
VCF_FILE="noca.urbanfixedsites.vcf"

# fst 
 while read -r GENE; do
    grep -E "($GENE\b|#)" "$VCF_FILE" | grep -E "#|HIGH|MODERATE" > "${GENE}.vcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 > "${GENE}.impact.tsv"
    COUNT=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 | wc -l)
    echo -e "$GENE\t$COUNT" >> genestats.txt
    rm "${GENE}.vcf"
 done < /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/noca.fst.50000kb.genenames_filtered.txt

awk '0 < $2 {print $1}' genestats.txt > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/noca.fst.50000kb.snpeff.genenames_filtered.txt
awk '0 < $2 {print $1}' genestats.txt > ~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/noca.fst.50000kb.snpeff.genenames_filtered.txt

# raisd 
 while read -r GENE; do
    grep -E "($GENE\b|#)" "$VCF_FILE" | grep -E "#|HIGH|MODERATE" > "${GENE}.vcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 > "${GENE}.impact.tsv"
    COUNT=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 | wc -l)
    echo -e "$GENE\t$COUNT" >> genestats.txt
    rm "${GENE}.vcf"
 done < /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/nocaurban.raisd.50kb.unique.genenames_filtered.txt


awk '0 < $2 {print $1}' genestats.txt > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/noca.raisd.50000kb.snpeff.genenames_filtered.txt
awk '0 < $2 {print $1}' genestats.txt > ~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/noca.raisd.50000kb.snpeff.genenames_filtered.txt


awk '0 < $2 {print $1}' genestats.txt  | wc -l
wc -l genestats.txt
wc -l  /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/nocaurban.raisd.50kb.unique.genenames_filtered.txt



# pyrr

#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural
# Query the genotypes for the urban samples and filter for sites where all urban genotypes are the same

bcftools view -HH -s PYRR003,PYRR004,PYRR006,PYRR007,PYRR009,PYRR011 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz | awk '{if (($10 ~ /^0\/0/ || $10 ~ /^1\/1/) && ($11 ~ /^0\/0/ || $11 ~ /^1\/1/) && ($12 ~ /^0\/0/ || $12 ~ /^1\/1/) && ($13 ~ /^0\/0/ || $13 ~ /^1\/1/) && ($14 ~ /^0\/0/ || $14 ~ /^1\/1/) && ($15 ~ /^0\/0/ || $15 ~ /^1\/1/)) print $0}' > pyrr_urban_same.txt



sbatch --account=mcnew \
            --job-name=pyrr_urban_same \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.pyrr_urban_same$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrr_urban_same.sh 

#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural


bcftools view -HH -s MSB25201,UWBM103346,UWBM77548,UWBM77780,UWBM77781 /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz | \
awk '{  
    # Store genotypes in an array, ignoring missing ones  
    delete gt; count = 0;  
    for (i = 10; i <= 14; i++) {  
        if ($i != "./.") {  
            gt[$i]++; count++;  
        }  
    }  

    # If more than one unique non-missing genotype exists, print line  
    if (length(gt) > 1) print $0;  
}' > pyrr_rural_diff.txt


sbatch --account=mcnew \
            --job-name=pyrr_rural_diff \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.pyrr_rural_diff$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/pyrr_rural_diff.sh 

# Extract the intersection of both urban same and rural differing sites
grep -F -f pyrr_urban_same.txt pyrr_rural_diff.txt > pyrrurbansame_ruraldiff.txt

awk -v OFS="\t" '{print $1,$2}' pyrr_urban_same.txt > pyrr_urban_same.sites.txt
awk -v OFS="\t" '{print $1,$2}' pyrr_rural_diff.txt > pyrr_rural_diff.sites.txt

sort pyrr_urban_same.sites.txt > pyrr_urban_same_sorted.txt
sort pyrr_rural_diff.sites.txt > pyrr_rural_diff_sorted.txt
comm -12 pyrr_urban_same_sorted.txt pyrr_rural_diff_sorted.txt > pyrrurbansame_ruraldiff.txt


#!/bin/bash

module load bcftools

cd /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural

# Use the combined sites list to extract those variants from the original VCF
bcftools view -R pyrrurbansame_ruraldiff.txt /xdisk/mcnew/dannyjackson/cardinals/datafiles/genotype_calls/cardinals_filtered.snpEffann.vcf.gz -o pyrr.urbanfixedsites.vcf

sbatch --account=mcnew \
            --job-name=pyrr_vcf \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.pyrr_vcf$.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=50gb \
            /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/pyrrurban_pyrrrural/filterurbansites.sh 


# see which genes are in these files
VCF_FILE="pyrr.urbanfixedsites.vcf"

# fst 
 while read -r GENE; do
    grep -E "($GENE\b|#)" "$VCF_FILE" | grep -E "#|HIGH|MODERATE" > "${GENE}.vcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 > "${GENE}.impact.tsv"
    COUNT=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 | wc -l)
    echo -e "$GENE\t$COUNT" >> genestats.txt
    rm "${GENE}.vcf"
 done < /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrr.fst.50000kb.genenames_filtered.txt

awk '0 < $2 {print $1}' genestats.txt > /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrr.fst.50000kb.snpeff.genenames_filtered.txt

# raisd 
 while read -r GENE; do
    grep -E "($GENE\b|#)" "$VCF_FILE" | grep -E "#|HIGH|MODERATE" > "${GENE}.vcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 > "${GENE}.impact.tsv"
    COUNT=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 | wc -l)
    echo -e "$GENE\t$COUNT" >> genestats.txt
    rm "${GENE}.vcf"
 done < /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrrurban.raisd.50kb.unique.genenames_filtered.txt


awk '0 < $2 {print $1}' genestats.txt  | wc -l
wc -l genestats.txt
wc -l  /xdisk/mcnew/dannyjackson/cardinals/analyses/genelist/final_gene_lists/pyrrurban.raisd.50kb.unique.genenames_filtered.txt




# Loop through each gene and count occurrences in the VCF file
for GENE in "${GENES[@]}"; do
    grep -E "($GENE\b|#)" "$VCF_FILE" | grep -E "#|HIGH|MODERATE" > "${GENE}.vcf"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 > "${GENE}.impact.tsv"
    COUNT=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' "${GENE}.vcf" | cut -d '|' -f 2,3,4,5,10 | wc -l)
    echo -e "$GENE\t$COUNT" >> genestats.txt
done

for GENE in "${GENES[@]}"; do
    rm "${GENE}.vcf"
done


grep '#' pyrr.urbanfixedsites.vcf > rho.vcf
grep 'RHO' pyrr.urbanfixedsites.vcf >> rho.vcf

bcftools view -i 'INFO/ANN~"HIGH" || INFO/ANN~"MODERATE"' rho.vcf > RHO_highmoderate.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ANN]\n' RHO_highmoderate.vcf | cut -d '|' -f 2,3,4,5,10 > RHO_impact.tsv

missense_variant|MODERATE|RHO|RHO|c.409G>A
missense_variant|MODERATE|RHO|RHO|c.613A>G
missense_variant|MODERATE|RHO|RHO|c.649G>A

grep 'CHROM' pyrr.urbanfixedsites.vcf
cat rho.vcf | awk '{printf "%s %s", $1, $2; for (i=10; i<=33; i++) printf " %s", substr($i,1,3); print ""}'
cat rho.vcf | awk '{printf "%s %s", $1, $2; for (i=10; i<=32; i++) printf " %s", substr($i,1,3); print ""}'

# 409G
NC_044582.1     19881243        .       G       A       3532.93 .       ANN=A|missense_variant|MODERATE|RHO|RHO|transcript|XP_030812390.1|protein_coding|2/5|c.409G>A|p.Val137Ile|409/1056|409/1056|137/351||,A|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XP_030812856.1|protein_coding||c.-4404G>A|||||4404|,A|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XM_030956996.1|pseudogene||n.-4404G>A|||||4404|,A|downstream_gene_variant|MODIFIER|IFT122|IFT122|transcript|XM_030956834.1|pseudogene||n.*4965G>A|||||4965|,A|downstream_gene_variant|MODIFIER|IFT122|IFT122|transcript|XM_030956835.1|pseudogene||n.*4965G>A|||||4965|,A|downstream_gene_variant|MODIFIER|IFT122|IFT122|transcript|XM_030956836.1|pseudogene||n.*4965G>A|||||4965|,A|non_coding_transcript_exon_variant|MODIFIER|RHO|RHO|transcript|XM_030956530.1|pseudogene|2/5|n.531G>A||||||       GT:PL:DP:SP:AD  1/1:158,21,0:7:0:0,7       1/1:104,9,0:3:0:0,3     ./.:0,0,0:0:0:0,0       1/1:248,30,0:10:0:0,10     1/1:135,12,0:4:0:0,4    1/1:104,9,0:3:0:0,3     1/1:67,6,0:2:0:0,2      1/1:160,18,0:6:0:0,6       1/1:255,33,0:11:0:0,11  1/1:182,18,0:6:0:0,6    1/1:74,6,0:2:0:0,2 1/1:118,12,0:4:0:0,4    1/1:135,12,0:4:0:0,4    1/1:135,12,0:4:0:0,4    1/1:207,21,0:7:0:0,7       1/1:183,21,0:7:0:0,7    1/1:143,18,0:6:0:0,6    1/1:161,15,0:5:0:0,5       1/1:185,18,0:6:0:0,6    1/1:94,9,0:3:0:0,3      1/1:216,27,0:9:0:0,9       1/1:95,9,0:3:0:0,3      1/1:225,24,0:8:0:0,8    1/1:172,18,0:6:0:0,6
# 613A
NC_044582.1     19881560        .       A       G       3177.93 .       ANN=G|missense_variant|MODERATE|RHO|RHO|transcript|XP_030812390.1|protein_coding|3/5|c.613A>G|p.Ile205Val|613/1056|613/1056|205/351||,G|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XP_030812856.1|protein_coding||c.-4087A>G|||||4087|,G|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XM_030956996.1|pseudogene||n.-4087A>G|||||4087|,G|non_coding_transcript_exon_variant|MODIFIER|RHO|RHO|transcript|XM_030956530.1|pseudogene|3/5|n.735A>G||||||       GT:PL:DP:SP:AD  1/1:255,39,0:13:0:0,13     1/1:129,15,0:5:0:0,5    ./.:0,0,0:0:0:0,0       1/1:179,18,0:6:0:0,6       1/1:159,18,0:6:0:0,6    1/1:67,6,0:2:0:0,2      1/1:73,9,0:3:0:0,31/1:57,6,0:2:0:0,2       1/1:113,12,0:4:0:0,4    1/1:161,15,0:5:0:0,5    1/1:67,6,0:2:0:0,2 1/1:149,15,0:5:0:0,5    1/1:165,18,0:6:0:0,6    1/1:160,15,0:5:0:0,5    1/1:165,18,0:6:0:0,6       1/1:235,30,0:10:0:0,10  1/1:160,15,0:5:0:0,5    1/1:179,18,0:6:0:0,6       1/1:179,18,0:6:0:0,6    1/1:67,6,0:2:0:0,2      1/1:208,24,0:8:0:0,8       1/1:37,3,0:1:0:0,1      1/1:37,3,0:1:0:0,1      1/1:200,21,0:7:0:0,7
# 649G
NC_044582.1     19881596        .       G       A       3200.93 .       ANN=A|missense_variant|MODERATE|RHO|RHO|transcript|XP_030812390.1|protein_coding|3/5|c.649G>A|p.Ala217Thr|649/1056|649/1056|217/351||,A|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XP_030812856.1|protein_coding||c.-4051G>A|||||4051|,A|upstream_gene_variant|MODIFIER|LOC115908437|LOC115908437|transcript|XM_030956996.1|pseudogene||n.-4051G>A|||||4051|,A|non_coding_transcript_exon_variant|MODIFIER|RHO|RHO|transcript|XM_030956530.1|pseudogene|3/5|n.771G>A||||||       GT:PL:DP:SP:AD  1/1:255,39,0:13:0:0,13     1/1:149,15,0:5:0:0,5    1/1:37,3,0:1:0:0,1      1/1:152,15,0:5:0:0,5       1/1:178,21,0:7:0:0,7    1/1:92,9,0:3:0:0,3      1/1:92,9,0:3:0:0,31/1:67,6,0:2:0:0,2       1/1:161,18,0:6:0:0,6    1/1:129,12,0:4:0:0,4    1/1:67,6,0:2:0:0,2 1/1:149,15,0:5:0:0,5    1/1:129,12,0:4:0:0,4    1/1:160,15,0:5:0:0,5    1/1:143,18,0:6:0:0,6       1/1:209,27,0:9:0:0,9    1/1:160,15,0:5:0:0,5    1/1:205,21,0:7:0:0,7       1/1:145,15,0:5:0:0,5    1/1:74,6,0:2:0:0,2      1/1:195,21,0:7:0:0,7       1/1:74,6,0:2:0:0,2      1/1:67,6,0:2:0:0,2      1/1:135,12,0:4:0:0,4