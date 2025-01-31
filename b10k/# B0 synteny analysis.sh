# B0 synteny analysis

conda activate ncbi_datasets
module load samtools
# parus
# datasets download genome accession GCF_001522545.3 --include gff3,rna,cds,protein,genome,seq-report

# zebra finch
datasets download genome accession GCF_003957565.2 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip

# zebra finch z = 	NC_044241.2; w = NC_045028.1
# (parus) z = 	NC_031799.1

#index a genome
# samtools faidx /xdisk/mcnew/dannyjackson/cardinals/synteny/ncbi_dataset/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna

samtools faidx /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/ncbi_dataset/data/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

#select chromosomes or regions
 
# samtools faidx /xdisk/mcnew/dannyjackson/cardinals/synteny/ncbi_dataset/data/GCF_001522545.3/GCF_001522545.3_Parus_major1.1_genomic.fna NC_031799.1 > zchrom_parusmajor.fa
samtools faidx /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/ncbi_dataset/data/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna NC_045028.1 > /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/wchrom_zebra.fa

samtools faidx /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/ncbi_dataset/data/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna NC_044241.2 > /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/zchrom_zebra.fa


# synteny with zebra z chrom

#!/bin/bash

#SBATCH --job-name=synteny1
#SBATCH --ntasks=47
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny1.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/wchrom_zebra.fa -o /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/output/wchrom_cardinals_zebra1

sbatch synteny1.sh 
Submitted batch job 11112958


#!/bin/bash

#SBATCH --job-name=synteny1.5
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny1.5.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/wchrom_zebra.fa -o /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/output/wchrom_cardinals_zebra15

sbatch synteny15.sh 
Submitted batch job 11112959

#!/bin/bash

#SBATCH --job-name=synteny1.25
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny1.25.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/wchrom_zebra.fa -o /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/output/wchrom_cardinals_zebra125

sbatch synteny125.sh
Submitted batch job 3413618

#!/bin/bash

#SBATCH --job-name=chrmap125
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.chrmap125.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Chromosemble -t -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/ncbi_dataset/data/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna -o /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/chromosomemap/cardinals_zebra125

Submitted batch job 3413653




#!/bin/bash

#SBATCH --job-name=synteny2
#SBATCH --ntasks=94
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny2.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -t /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -q /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/wchrom_zebra.fa -o /xdisk/mcnew/dannyjackson/cardinals/synteny/zebra/output/wchrom_cardinals_zebra2

sbatch synteny2.sh 
Submitted batch job 11112961




# synteny with parus z chrom

#!/bin/bash

#SBATCH --job-name=synteny
#SBATCH --ntasks=12
#SBATCH --nodes=10
#SBATCH --time=100:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/ncbi_dataset/data/GCF_001522545.3/zchrom_parusmajor.fa -o /xdisk/mcnew/dannyjackson/cardinals/synteny/output2/zchrom_cardinals_parusmajor

sbatch zchromsynt.sh 
Submitted batch job 3413523

















# if needed, full genome synteny


#!/bin/bash

#SBATCH --job-name=synteny
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.synteny.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Satsuma -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna -t /xdisk/mcnew/dannyjackson/cardinals/synteny/ncbi_dataset/data/GCA_014549065.1/GCA_014549065.1_CCar_1.0_genomic.fna -o /xdisk/mcnew/dannyjackson/cardinals/synteny/output/xcorr_cardinals

sbatch synteny.sh 
Submitted batch job 3413514


awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna.fai > /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/chromosomelist.txt


#!/bin/bash

#SBATCH --job-name=chrmap
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=100:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.chrmap.%j

/xdisk/mcnew/dannyjackson/cardinals/synteny/satsuma-code/Chromosemble -t /xdisk/mcnew/dannyjackson/cardinals/synteny/ncbi_dataset/data/GCA_014549065.1/GCA_014549065.1_CCar_1.0_genomic.fna -q /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/chromosomelist.txt -o /xdisk/mcnew/dannyjackson/cardinals/synteny/chromosomemap

sbatch chrmap.sh 
Submitted batch job 3413516