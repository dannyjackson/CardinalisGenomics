# rho investigation

/xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff
grep 'RHO' /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff

VYXE01006626.1:115846-118575

module load bcftools

~/programs/angsd/angsd -bam angsd.both.pops -out RPE65 -setMinDepthInd 2 -setMaxDepthInd 12 -doBcf 1 -doCounts 1 -doGeno 4 -doPost 1 -gl 1 -doMajorMinor 1 -doMaf 2 -r NC_087518.1:40880137-40888431 -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna

bcftools mpileup -Ou -f /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -a FORMAT/AD,DP,INFO/AD,SP -r VYXE01006626.1:115846-118575 /xdisk/mcnew/dannyjackson/cardinals/indelrealignment/*bam | bcftools call -mv > genocalls.vcf

bcftools index angsdput.bcf

bcftools convert angsdput.bcf > angsdput.vcf

bcftools query -f '%CHROM %POS %REF %ALT{0} GTs:[ %GT]\n' genocalls.vcf | grep '0/1'

bcftools query -f '%CHROM %POS %REF %ALT{0} GTs:[ %GT]\n' genocalls.vcf | grep '1/1' > fixedsnps_rho.txt

mkdir /xdisk/mcnew/dannyjackson/copythis

cp fixedsnps_rho.txt /xdisk/mcnew/dannyjackson/copythis

MSB25201    PYRR      
NOCA003 NOCA
NOCA004 NOCA
NOCA006 NOCA
NOCA008 NOCA
NOCA012 NOCA
NOCA013 NOCA
PYRR003 PYRR
PYRR004 PYRR
PYRR006 PYRR
PYRR007 PYRR
PYRR009 PYRR
PYRR011 PYRR
BM100619    NOCA
BM100620    NOCA
BM100621    NOCA
BM103345    NOCA
BM103346    PYRR
BM77548 PYRR
BM77718 PYRR
BM77780 PYRR
BM77781 PYRR
BM77856 NOCA
BM77978 NOCA