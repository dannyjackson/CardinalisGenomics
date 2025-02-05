# params_tajima
source /home/u15/dannyjackson/programs/CardinalisGenomics/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=0.001

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1="#6247aa"
COLOR2="#b185db"

# define the name of the population that will be analyzed
POP=nocarural

# source the setup file for tajima
source ${SCRIPTDIR}/Genomics-Main/tajima/setup_tajima.sh

