# B3_dxy
# after all of this, it's only identifying sex chromosome differences. Which really sucks because I have batch biased sex groupings. I'll look into satsuma again i guess and try to remove sex chromosomes.

# noca

module load R/4.4.0

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.mafs -q /xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural.mafs -t 2618621

mv Dxy_persite.txt Dxy_persite_noca.txt

[1] "Created Dxy_persite.txt"
[1] "Global dxy is: 437265.040009416"
[1] "Global per site Dxy is: 0.166982942552365"

awk '{print $2}' Dxy_persite_noca.txt  | tail -n +2 > noca_sites.txt

tail -n +2 Dxy_persite_noca.txt > Dxy_persite_noca.forplot.txt
sed -i 's/dxy/fst/g' Dxy_persite_noca.forplot.txt

Rscript plotDXY.R -i Dxy_persite_noca.forplot.txt -o noca_windowed -p noca_sites.txt -w 50000 -s 1 

# sed -i 's/fst/dxy/g' Dxy_persite_noca.txt
awk '{print $1}' Dxy_persite_noca.txt | sort -u | tail -n +2 > chroms.txt
awk 'BEGIN {OFS = "\t"} {print $1,$2}' /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/GCA_013397215.1_ASM1339721v1_genomic.fna.fai > allchroms.txt

while read -r chrom;
do
grep -m 1 "$chrom" allchroms.txt >> chroms.len.txt

done < chroms.txt


# write the following python script to file 50kbwin_noca.py
import pandas as pd

# Load the chromosome lengths
chroms_len = pd.read_csv("chroms.len.txt", sep="\t", header=None, names=["chromosome", "length"])

# Load the dxy data
dxy_data = pd.read_csv("Dxy_persite_noca.txt", sep="\t")

# Initialize an empty list to store the window data
window_data = []

# Define the window size
window_size = 50000

# Process each chromosome
for _, chrom_row in chroms_len.iterrows():
    chrom = chrom_row["chromosome"]
    chrom_length = chrom_row["length"]
    
    # Filter the dxy data for the current chromosome
    chrom_dxy = dxy_data[dxy_data["chromo"] == chrom]
    
    # Divide the chromosome into windows
    for start in range(0, chrom_length, window_size):
        end = min(start + window_size, chrom_length)
        mid = (start + end) // 2
        win_len = (end - start)
        # Filter dxy data for the current window
        window_sites = chrom_dxy[(chrom_dxy["position"] >= start) & (chrom_dxy["position"] < end)]
        
        # Calculate statistics
        total_sites = len(window_sites)
        average_dxy = window_sites["dxy"].mean() if total_sites > 3 else None
        
        # Append the result
        window_data.append([chrom, start, end, mid, win_len, total_sites, average_dxy])

# Convert the result into a DataFrame
output_df = pd.DataFrame(window_data, columns=["chromosome", "start", "end", "mid", "win_size", "total_sites", "average_dxy"])

# Save to a new file
output_df.to_csv("average_dxy_50kb_windows.noca.txt", sep="\t", index=False, na_rep="NA")




# write to file 50kbwin_pyrr.py
import pandas as pd

# Load the chromosome lengths
chroms_len = pd.read_csv("chroms.len.txt", sep="\t", header=None, names=["chromosome", "length"])

# Load the dxy data
dxy_data = pd.read_csv("Dxy_persite_pyrr.txt", sep="\t")

# Initialize an empty list to store the window data
window_data = []

# Define the window size
window_size = 50000

# Process each chromosome
for _, chrom_row in chroms_len.iterrows():
    chrom = chrom_row["chromosome"]
    chrom_length = chrom_row["length"]
    
    # Filter the dxy data for the current chromosome
    chrom_dxy = dxy_data[dxy_data["chromo"] == chrom]
    
    # Divide the chromosome into windows
    for start in range(0, chrom_length, window_size):
        end = min(start + window_size, chrom_length)
        mid = (start + end) // 2
        win_len = (end - start)
        # Filter dxy data for the current window
        window_sites = chrom_dxy[(chrom_dxy["position"] >= start) & (chrom_dxy["position"] < end)]
        
        # Calculate statistics
        total_sites = len(window_sites)
        average_dxy = window_sites["dxy"].mean() if total_sites > 3 else None
        
        # Append the result
        window_data.append([chrom, start, end, mid, win_len, total_sites, average_dxy])

# Convert the result into a DataFrame
output_df = pd.DataFrame(window_data, columns=["chromosome", "start", "end", "mid", "win_size", "total_sites", "average_dxy"])

# Save to a new file
output_df.to_csv("average_dxy_50kb_windows.pyrr.txt", sep="\t", index=False, na_rep="NA")




#!/bin/bash

#SBATCH --job-name=noca_win_dxy
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.noca_win_dxy.%j

cd /xdisk/mcnew/dannyjackson/cardinals/dxy
module load python/3.11/3.11.4
python 50kbwin_noca.py

sbatch submit_noca_win_dxy.sh 
Submitted batch job 3413441


#!/bin/bash

#SBATCH --job-name=pyrr_win_dxy
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pyrr_win_dxy.%j

cd /xdisk/mcnew/dannyjackson/cardinals/dxy
module load python/3.11/3.11.4
python 50kbwin_pyrr.py


sbatch submit_pyrr_win_dxy.sh 
Submitted batch job 3413442






# Plot dxy noca
install.packages("readr")
install.packages("ggrepel")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("RColorBrewer")

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('average_dxy_50kb_windows.noca.txt', sep ='\t')
dxy_noNA <- na.omit(dxy)
nrow(dxy) - nrow(dxy_noNA)
# 2999
dxy = dxy_noNA

dxy = dxy[dxy$win_size>40000,]

max(dxy$average_dxy)
# 0.5090998
min(dxy$average_dxy)
# 0

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(average_dxy)) 

nrow(dxy)
# 17501 windows, so top 0.1% would be 18

outlier_dxy_disorder <- ordered_dxy[1:18,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromosome, mid)

min(outlier_dxy_disorder$average_dxy)
# 0.3808062
max(outlier_dxy_disorder$average_dxy)
# 0.5090998
write.csv(outlier_dxy, "outlierdxy_noca_windows.csv")

colnames(dxy) <- c("chromo", "start", "end", "position", "win_size", "total_sites", "dxy")

# draw it with cutoff line 


blues <- c("#4EAFAF", "#082B64")
df.tmp <- dxy %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dxy, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("noca.dxy.windowed.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 6766 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.3808062) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()

# plot dxy pyrr

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('average_dxy_50kb_windows.pyrr.txt', sep ='\t')
dxy_noNA <- na.omit(dxy)
nrow(dxy) - nrow(dxy_noNA)
# 2999
dxy = dxy_noNA

dxy = dxy[dxy$win_size>40000,]

max(dxy$average_dxy)
# 0.4923265
min(dxy$average_dxy)
# 0

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(average_dxy)) 

nrow(dxy)
# 17501 windows, so top 0.1% would be 18

outlier_dxy_disorder <- ordered_dxy[1:18,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromosome, mid)

min(outlier_dxy_disorder$average_dxy)
# 0.4062524
max(outlier_dxy_disorder$average_dxy)
# 0.4923265
write.csv(outlier_dxy, "outlierdxy_pyrr_windows.csv")

colnames(dxy) <- c("chromo", "start", "end", "position", "win_size", "total_sites", "dxy")

# draw it with cutoff line 


blues <- c("#4EAFAF", "#082B64")
df.tmp <- dxy %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dxy, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("pyrr.dxy.windowed.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 6766 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.4062524) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()







# make relevant gene lists

awk 'BEGIN {FS = ","} {$1=""}1' outlierdxy_noca_windows.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > outlierdxy_noca.csv.tmp

mv outlierdxy_noca.csv.tmp outlierdxy_noca_windows.csv

awk 'BEGIN {FS = ","} {$1=""}1' outlierdxy_pyrr_windows.csv | awk 'BEGIN {OFS = ","} {$1=$1}1 ' > outlierdxy_pyrr.csv.tmp

mv outlierdxy_pyrr.csv.tmp outlierdxy_pyrr_windows.csv

sed -i 's/\"//g' outlierdxy_noca_windows.csv
sed -i 's/\"//g' outlierdxy_pyrr_windows.csv


# noca

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $2}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_noca_windowed_dxy.txt

done < outlierdxy_noca_windows.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_noca_windowed_dxy.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_noca_windowed_dxy.txt

# pyrr

while read -r line;
do

chr=`awk 'BEGIN {FS = ","} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = ","} {print $2}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))


grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $5 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_pyrr_windowed_dxy.txt

done < outlierdxy_pyrr_windows.csv


awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' relevantgenes_pyrr_windowed_dxy.txt | sed 's/Name\=//g' | sort -u > relevantgenenames_pyrr_windowed_dxy.txt