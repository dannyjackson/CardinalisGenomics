# Load required modules
module load R/4.4.0
module load htslib/1.19.1
module load bedtools2/2.29.2
module load python/3.11/3.11.4

# Define paths
MSMCTOOLS=~/programs/msmc-tools  # Directory with MSMC-tools binaries
ANGSD=~/programs/angsd/           # Path to ANGSD executables
SCRIPTDIR=~/programs/CardinalisGenomics/

# Update PATH
export PATH="$PATH:$MSMCTOOLS:$SCRIPTDIR"

# Define output directory
OUTDIR=/xdisk/mcnew/dannyjackson/cardinals/

# Define population names
POP1=nocaurban
POP2=nocarural

# Create necessary directories
mkdir -p ${OUTDIR}/analyses/{fst,genelist,dxy}
mkdir -p ${OUTDIR}/datafiles/{safs,mls}
mkdir -p ${OUTDIR}/analyses/fst/${WIN}
mkdir -p ${OUTDIR}/analyses/genelist/${WIN}
mkdir -p ${OUTDIR}/analyses/dxy/${POP1}_${POP2}

# Define colors
COLOR1="#7570B3"
COLOR2="#CAB2D6"

# Chromosome-related definitions
CHRLEAD="NC_0"
SEXCHR="NC_044601"

# Reference genome and annotation paths
REF=${OUTDIR}/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
GFF=${OUTDIR}/datafiles/referencegenome/ncbi_dataset/data/GCA_013397215.1/genomic.gff

# Create reference chromosome lengths file
LENGTHS_FILE=${OUTDIR}/referencelists/autosomes_lengths.txt
if [ ! -f "$LENGTHS_FILE" ]; then
    awk 'BEGIN {OFS="\t"} {print $1, $2}' ${REF}.fai | grep ${CHRLEAD} | grep -v ${SEXCHR} > $LENGTHS_FILE
    while IFS=',' read -r first second; do
        sed -i "s/$second/$first/g" $LENGTHS_FILE
    done <<< "$CHROM"
else
    echo "Chromosome lengths table already complete, moving on!"
fi

# Create chromosome conversion file
CONVERSION_FILE=${OUTDIR}/referencelists/GCF_901933205_chromconversion.txt
if [ ! -f "$CONVERSION_FILE" ]; then
    cat << EOF > $CONVERSION_FILE
1,NC_044571.1
2,NC_044572.1
3,NC_044573.1
4,NC_044574.1
5,NC_044575.1
6,NC_044576.1
7,NC_044577.1
8,NC_044578.1
9,NC_044579.1
10,NC_044580.1
11,NC_044581.1
12,NC_044582.1
13,NC_044583.1
14,NC_044584.1
15,NC_044585.1
1A,NC_044586.1
17,NC_044587.1
18,NC_044588.1
19,NC_044589.1
20,NC_044590.1
21,NC_044591.1
22,NC_044592.1
23,NC_044593.1
24,NC_044594.1
25,NC_044595.1
26,NC_044596.1
27,NC_044597.1
28,NC_044598.1
29,NC_044599.1
4A,NC_044600.1
Z,NC_044601.1
EOF
else
    echo "Chromosome conversion table already complete, moving on!"
fi