# gene visualization
cd /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms

mkdir -p plot

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r both.all.GO.nocorrection

# local testing

ls pantheroutput/*b.txt | xargs -n 1 basename > /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/GOlistnames.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/pyrr.fst.50kb.snpeff
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/pyrr.raisd.50kb.snpeff
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/noca.fst.50kb.snpeff
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/noca.raisd.50kb.snpeff

grep -o 'GO:[0-9]\{7\}' pantheroutput/noca.raisd.50kb.snpeff.fdr.txt > noca.raisd.50kb.snpeff.fdr.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r noca.fstraisd.50kb.fdr.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r noca.fst.50kb.fdr.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r noca.raisd.50kb.fdr.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r pyrr.fst.50kb.fdr.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r pyrr.raisd.50kb.fdr.GO


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/pyrr.fst.50kb
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/bothspecies
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/fst
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/raisd

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/overlappinggenes
grep -o 'GO:[0-9]\{7\}' pantheroutput/overlappinggenes.fdr.txt > overlappinggenes.fdr.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r overlappinggenes.fdr

grep -o 'GO:[0-9]\{7\}' pantheroutput/noca.raisd.50kb.snpeff.fdr.txt
 > ${file}

while read -r file;
do

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/${file%.txt}

grep -o 'GO:[0-9]\{7\}' pantheroutput/${file%.txt}.fdr.txt > ${file%.txt}.fdr.GO.txt

done < /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/GOlistnames.txt


while read -r file;
do

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO

done < /Users/danjack/Documents/Github_local/CardinalisGenomics/D_CandidateGeneAnalyses/GOterms/GOlistnames.txt

cd plot/
mv *heatmap.png heatmap/
mv *scatter.png scatter/
mv *treemap.png tree/

mv *heatmap.pdf heatmap/
mv *scatter.pdf scatter/
mv *treemap.pdf tree/

grep -o 'GO:[0-9]\{7\}' pantheroutput/noca.fstraisd.50kb.fdr.txt > noca.fstraisd.50kb.fdr.GO.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r noca.fstraisd.50kb.fdr.GO


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r pyrr.fst.50kb.all.GO

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r both.fst.50kb.all.GO.fdr02

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r noca.fst.50kb.all.GO












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





# visualize overlapping GO terms across analyses

# Analyze overlap between methods
library(ggvenn)

df_noca_fstraisd<-read.csv("noca.fstraisd.50kb.fdr.GO.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_noca_fst<-read.csv("noca.fst.50kb.fdr.GO.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_noca_raisd<-read.csv("noca.raisd.50kb.fdr.GO.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

# 

df_pyrr_fst<-read.csv("pyrr.fst.50kb.fdr.GO.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_pyrr_raisd<-read.csv("pyrr.raisd.50kb.fdr.GO.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]


# Compute overlaps
# none across noca analyses
noca_fst_raisd = intersect(df_noca_fst, df_noca_raisd) # none
noca_fst_raisd_all <- Reduce(intersect, list(df_noca_fst, df_noca_fstraisd, df_noca_raisd))

pyrr_fst_raisd = intersect(df_pyrr_fst, df_pyrr_raisd) # [1] "GO:0031099" "GO:0042246" # regeneration and tissue regeneration


gene_sets <- list(
  pyrr_fst = df_pyrr_fst,
  pyrr_raisd = df_pyrr_raisd
)


# pdf(file = "noca_windows.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("pyrr_fst", "pyrr_raisd"),
  stroke_size = 0.5, set_name_size = 4
  )

# dev.off()


# compare between noca and pyrr
fst = intersect(df_noca_fst, df_pyrr_fst) # none
raisd = intersect(df_noca_raisd, df_pyrr_raisd)
# "GO:0008033" "GO:0006399" tRNA processing, tRNA metabolic process