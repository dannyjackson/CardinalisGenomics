# params_fst
source /home/u15/dannyjackson/programs/CardinalisGenomics/params_base.sh

ANGSD=~/programs/angsd/ # path to directory with angsd executables

CUTOFF=0.001

# define two colors to be used (alternating chromosomes in manhattan plots)
COLOR1="#2d6a4f"
COLOR2="#74c69d"

# define the names of the two populations that will be compared
POP1=pyrrurban
POP2=pyrrrural

# source the setup file for fst
source ${SCRIPTDIR}/Genomics-Main/fst/setup_fst.sh
