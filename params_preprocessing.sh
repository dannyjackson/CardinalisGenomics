module load bwa/0.7.17
module load samtools/1.19.2
module load bowtie2
module load picard
module load samtools
module load parallel
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9

source ~/programs/CardinalisGenomics/params_base.sh
source ~/programs/CardinalisGenomics/Genomics-Main/A_Preprocessing/preprocessing_setup.sh

# across all preprocessing
THREADS=8

# trimming
FASTAS=/xdisk/mcnew/dannyjackson/cardinals/condensed_fastas/ # format must be samplename_
TRIMJAR=~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar
LEAD=20 # value to trim from leading strand, often 20
TRAIN=20 # value to trim from trailing strand, often 20
SLIDE=4:20 # threshold and windlow length, often 4:20
MINREADLEN=90 # minimum length for a read to be kept, often 90 for 150bp sequencing

# clipping
BAMUTILBAM=/home/u15/dannyjackson/programs/bamUtil-master/bin/bam

# bam statistics

# snp ID
ANGSD=~/programs/angsd/ # path to directory with angsd executables
SNPPVAL=1e-6 # max p-value for snp to be considered significant, often 1e-6
MINDEPTHIND=4 # minimum depth per individual required for a site to be kept
MININD=20 # minimum number of individuals required for a site to be kept
MINQ=30 # minimum quality score required for a site to be kept
MINMAF=0.05 # minimum minor allele frequency required for a site to be kept
MINMAPQ=30 # minimum mapping quality score required for a site to be kept

# generate mask
k=150
prefix=GCF_901933205
