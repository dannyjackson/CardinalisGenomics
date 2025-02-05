# params_tajima
source /home/u15/dannyjackson/programs/CardinalisGenomics/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=0.001

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1="#7570B3"
COLOR2="#CAB2D6"

# define the names of the two populations that will be compared
POP=nocaurban

# source the setup file for tajima
source ${SCRIPTDIR}/Genomics-Main/tajima/setup_tajima.sh

