# gene visualization
cd /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms

mkdir -p plot


# local testing

ls pantheroutput/pyrr*txt | xargs -n 1 basename > /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/GOlistnames.txt


while read -r file;
do

grep -o 'GO:[0-9]\{7\}' pantheroutput/${file} > ${file}

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}

done < /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/GOlistnames.txt















# I could not get rrvgo installed on the hpc and decided to just run this locally for now. Keep troubleshooting, see below for notes.
/
find ~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/* | awk -F"[./]" '{ OFS="." ; print $(NF-6), $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1) }' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/GOlistnames.txt




while read -r file;
do

    sbatch --account=mcnew \
            --job-name=visualization_${file} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.visualization_${file}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=1:00:00 \
            ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_visualization.sh \
            -p ~/programs/CardinalisGenomics/params_base.sh \
            -f $file 

done < /xdisk/mcnew/dannyjackson/cardinals/referencelists/GOlistnames.txt



# testing remote

source ~/programs/CardinalisGenomics/params_base.sh


chmod +x ~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_visualization.sh

~/programs/CardinalisGenomics/Genomics-Main/D_GeneVisualization/D2_visualization.sh -p ~/programs/CardinalisGenomics/params_base.sh -f nocaurban.raisd.50kb.unique.GOterms.fdr

install.packages("igraph", type = "source")
install.packages("xml2", type = "source")

pak::pak("igraph/rigraph")
pak::pak("r-lib/xml2")

install.packages("tm", type = "source")

rrvgo

* DONE (umap)
ERROR: dependency ‘igraph’ is not available for package ‘treemap’
* removing ‘/home/u15/dannyjackson/R/x86_64-pc-linux-gnu-library/4.4/treemap’
ERROR: dependency ‘xml2’ is not available for package ‘tm’
* removing ‘/home/u15/dannyjackson/R/x86_64-pc-linux-gnu-library/4.4/tm’
ERROR: dependencies ‘treemap’, ‘tm’ are not available for package ‘rrvgo’
* removing ‘/home/u15/dannyjackson/R/x86_64-pc-linux-gnu-library/4.4/rrvgo’

# Load required packages, installing if necessary

required_packages <- c("BiocManager", "limma", "rrvgo")

installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
scriptdir <- args [2]
file <- args[3]
