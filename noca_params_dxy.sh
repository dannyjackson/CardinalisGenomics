# params_dxy
source /home/u15/dannyjackson/programs/CardinalisGenomics/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=0.001

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1="#6247aa"
COLOR2="#b185db"

# define the names of the two populations that will be compared
POP1=nocaurban
POP2=nocarural

# source the setup file for dxy
source ${SCRIPTDIR}/Genomics-Main/dxy/setup_dxy.sh
