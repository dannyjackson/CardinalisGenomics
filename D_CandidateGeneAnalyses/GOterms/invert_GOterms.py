# flip GO term table
import pandas as pd

# Load data from file, ensuring internal commas within quotes are handled correctly
df = pd.read_csv('pyrr.raisd.50kb.fdr.GO.genelist.csv', sep=',', quotechar='"')

# Splitting the 'Genes' column into lists
df["Genes"] = df["Genes"].str.split(",")

# Exploding the dataframe so each gene gets its own row
df = df.explode("Genes")

# Pivoting the table
df_pivot = df.groupby("Genes")["GO_Term"].apply(lambda x: ", ".join(x)).reset_index()

# Save to file
df_pivot.to_csv('pyrr.raisd.50kb.fdr.GO.genelist.transformed.csv', sep=',', index=False)



# Define the VCF file path
VCF_FILE="pyrrurban_pyrrrural.fst.50000.outlier.snpeff.genes_filtered.txt"

# List of genes
GENES=(
ABI2
ANKIB1
ARIH1
BABAM2
CACTIN
CDK14
CDKAL1
CHEK1
CSGALNACT1
CTSC
CTSO
DCUN1D1
DMD
DROSHA
DTNBP1
DYRK1A
EI24
EIF3F
ENTPD1
EPHA10
ERBB4
FUT8
G6PC2
GAP43
GGT7
GIGYF2
GPM6A
HEXA
IGF1
IGF2BP2
KIRREL3
KYNU
LOXL4
MAN1A2
MTA1
MTMR4
NAV2
NCF4
NDUFA9
NLGN3
NOC3L
NTN4
NUP85
NYAP2
OTUD7A
PBRM1
PDPK1
PIP5K1C
PLCE1
PPP2R3B
PPP6R2
PPP6R3
PRKCH
RAB19
RBFOX1
RNF13
RNF43
RNGTT
RUNX2
SDK1
SEMA6D
SENP2
SLIT3
SNX7
STK4
STT3A
SUMO2
TASP1
TDO2
TENM3
TRAK1
TRIP12
UBA6
UBE2I
WNT16
)

# Loop through each gene and count occurrences in the VCF file
for GENE in "${GENES[@]}"; do
    COUNT=$(grep -w "$GENE" "$VCF_FILE" | wc -l)
    echo -e "$GENE\t$COUNT"
done




# Manual analysis of nonsignificant after FDr but sig before GO terms
# flip GO term table
import pandas as pd

# Load data from file, ensuring internal commas within quotes are handled correctly
df = pd.read_csv('pyrr.fstraisd.50kb.manual.GO.flipped.csv', sep=',', quotechar='"')

# Splitting the 'Genes' column into lists
df["Genes"] = df["Genes"].str.split(", ")

# Exploding the dataframe so each gene gets its own row
df = df.explode("Genes")

# Pivoting the table
df_pivot = df.groupby("Genes")["GO_Term"].apply(lambda x: ", ".join(x)).reset_index()

# Save to file
df_pivot.to_csv('pyrr.fstraisd.50kb.manual.GO.flipped.transformed.csv', sep=',', index=False)

# Display result
display(df_pivot)


#!/bin/bash

# Define the VCF file path
VCF_FILE="pyrrurban_pyrrrural.fst.50000.outlier.ann.maf.vcf"

# List of genes
GENES=(
    AKNA CCDC80 NECTIN3 PLXND1 RELL2 TEN1 TENM3 CRY1 FSHB ID4 MTA1 RHO WWTR1 YAP1
    GIGYF2 IGF1 BANK1 IL15 FAM49B CASQ2 DMD FHOD3 POPDC2 SGCG ABI2 DCDC2 KIAA0319
    KIRREL3 NLGN3 NTN3 SDK1 GAP43
)

# Loop through each gene and count occurrences in the VCF file
for GENE in "${GENES[@]}"; do
    COUNT=$(grep -w "$GENE" "$VCF_FILE" | wc -l)
    echo -e "$GENE\t$COUNT"
done




# noca
# flip GO term table
import pandas as pd

# Load data from file, ensuring internal commas within quotes are handled correctly
df = pd.read_csv('noca.fstraisd.50kb.manual.GO.flipped.csv', sep=',', quotechar='"')

# Splitting the 'Genes' column into lists
df["Genes"] = df["Genes"].str.split(",")

# Exploding the dataframe so each gene gets its own row
df = df.explode("Genes")

# Pivoting the table
df_pivot = df.groupby("Genes")["GO_Term"].apply(lambda x: ", ".join(x)).reset_index()

# Save to file
df_pivot.to_csv('noca.fstraisd.50kb.manual.GO.flipped.transformed.csv', sep=',', index=False)


ARFIP1 CACTIN CDH11 CDH18 CDH20 CDH9 CLDN23 DMD ECM2 FGA FGB IPO8 IPO9 MYF5 MYF6 NLGN3 PPARG RAB32 RAB34 SDK1 SGCZ SLC6A1 SPATA7 SUSD4 TJP3 TLR5 VIPAS39


# Define the VCF file path
VCF_FILE="~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/noca.fst.50000kb.snpeff.genenames_filtered.txt"

# List of genes
GENES=(
    ARFIP1 CACTIN CDH11 CDH18 CDH20 CDH9 CLDN23 DMD ECM2 FGA FGB IPO8 IPO9 MYF5 MYF6 NLGN3 PPARG RAB32 RAB34 SDK1 SGCZ SLC6A1 SPATA7 SUSD4 TJP3 TLR5 VIPAS39
)

# Loop through each gene and count occurrences in the VCF file
for GENE in "${GENES[@]}"; do
    COUNT=$(grep -w "$GENE" ~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/noca.fst.50000kb.snpeff.genenames_filtered.txt | wc -l)
    echo -e "$GENE\t$COUNT"
done

cd ~/programs/CardinalisGenomics/D_CandidateGeneAnalyses/GeneLists/

for GENE in "${GENES[@]}"; do
    echo -e "$GENE"
    grep -w "$GENE" *
done



