# change in allele frequency

cd /xdisk/mcnew/dannyjackson/cardinals/deltaAF
/xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.mafs 
/xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural.mafs 
"chromo", "position", "major", "minor", "anc", "knownEM", "nInd"

module load R/4.4.0
# Load necessary library
library(dplyr) 

# Read the input files
file1 <- read.csv("/xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural.mafs", sep = '\t')
file2 <- read.csv("/xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.mafs", sep = "\t")

# Perform the calculation
result <- file1 %>%
  left_join(file2, by = c("chromo", "position"), suffix = c(".file1", ".file2")) %>%
  mutate(
    diff = case_when(
      major.file1 == major.file2 ~ knownEM.file1 - knownEM.file2,
      major.file1 == minor.file2 ~ (knownEM.file1 - knownEM.file2) * -1,
      TRUE ~ NA_real_
    )
  )

# View the result
# print(result)

# Optionally save to file
write.csv(result, "noca_deltaAF.csv", row.names = FALSE)


# plot it 
dAF <- read.csv("noca_deltaAF.csv")
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dAF_noNA <- na.omit(dAF)
nrow(dAF) - nrow(dAF_noNA)
# 2999
dAF = dAF_noNA

max(dAF$diff)
# 0.915149
min(dAF$diff)
# -0.915149

ordered_dAF <- dAF %>% 
 # desc orders from largest to smallest
 arrange(desc(diff)) 

nrow(dAF)
# 2618621 snps, so top 0.1% would be 2619

outlier_dAF_disorder <- ordered_dAF[1:2619,]

outlier_dAF <- outlier_dAF_disorder %>% arrange(chromo, position)

min(outlier_dAF_disorder$diff)
# 0.570472
max(outlier_dAF_disorder$diff)
# 1
write.csv(outlier_dAF, "outlierdAF_noca_snps.csv")

# colnames(dAF) <- c("chromo", "start", "end", "position", "win_size", "total_sites", "dAF")

# draw it with cutoff line 


blues <- c("#4EAFAF", "#082B64")
df.tmp <- dAF %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dAF, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("noca.dAF.snps.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(diff))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1309310)) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dAF") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.611649) +
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











# windowed
awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.mafs | sort -u | tail -n +2 > chroms.txt
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
dxy_data = pd.read_csv("noca_deltaAF.csv", sep=",")

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
        average_deltaAF = window_sites["diff"].mean() if total_sites > 3 else None
        
        # Append the result
        window_data.append([chrom, start, end, mid, win_len, total_sites, average_deltaAF])

# Convert the result into a DataFrame
output_df = pd.DataFrame(window_data, columns=["chromosome", "start", "end", "mid", "win_size", "total_sites", "average_deltaAF"])

# Save to a new file
output_df.to_csv("average_deltaAF_50kb_windows.noca.txt", sep="\t", index=False, na_rep="NA")






# plot it 
dAF <- read.csv("average_deltaAF_50kb_windows.noca.txt", sep = '\t')
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dAF_noNA <- na.omit(dAF)
nrow(dAF) - nrow(dAF_noNA)
# 2999
dAF = dAF_noNA

max(dAF$average_deltaAF)
# 0.3071835
min(dAF$average_deltaAF)
# -0.2310722

ordered_dAF <- dAF %>% 
 # desc orders from largest to smallest
 arrange(desc(average_deltaAF)) 

nrow(dAF)
# 27529 snps, so top 0.1% would be 28

outlier_dAF_disorder <- ordered_dAF[1:28,]

outlier_dAF <- outlier_dAF_disorder %>% arrange(chromosome, mid)

min(outlier_dAF_disorder$average_deltaAF)
# 0.1439908
max(outlier_dAF_disorder$average_deltaAF)
# 1
write.csv(outlier_dAF, "outlierdAF_noca_windows.csv")

colnames(dAF) <- c("chromo", "start", "end", "position", "win_size", "total_sites", "dAF")

# draw it with cutoff line 


blues <- c("#4EAFAF", "#082B64")
df.tmp <- dAF %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dAF, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("noca.dAF.window.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dAF))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1309310)) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dAF") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.1439908) +
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
