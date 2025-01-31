# B4_tajimasd

module load R/4.4.0

cd /xdisk/mcnew/dannyjackson/cardinals/tajima

# noca urban
~/programs/angsd/misc/realSFS -P 24 /xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.saf.idx > noca_urban.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/cardinals/fst/noca_urban.saf.idx -outname noca_urban -sfs noca_urban.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat noca_urban.thetas.idx

# 50 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat noca_urban.thetas.idx -win 50000 -step 1  -outnames nocaurban.theta.thetasWindow

# noca rural
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural.saf.idx > noca_rural.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/cardinals/fst/noca_rural.saf.idx -outname noca_rural -sfs noca_rural.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat noca_rural.thetas.idx

# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat noca_rural.thetas.idx -win 50000 -step 1  -outnames nocarural.theta.thetasWindow





# pyrr


# pyrr urban
~/programs/angsd/misc/realSFS -P 24 /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_urban.saf.idx > pyrr_urban.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_urban.saf.idx -outname pyrr_urban -sfs pyrr_urban.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat pyrr_urban.thetas.idx

# 50 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat pyrr_urban.thetas.idx -win 50000 -step 1  -outnames pyrrurban.theta.thetasWindow

# pyrr rural
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_rural.saf.idx > pyrr_rural.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/cardinals/fst/pyrr_rural.saf.idx -outname pyrr_rural -sfs pyrr_rural.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat pyrr_rural.thetas.idx

# 50 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat pyrr_rural.thetas.idx -win 50000 -step 1  -outnames pyrrrural.theta.thetasWindow
























# NOCA
# noca urban 

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > nocaurban.thetasWindow.foroutliers
grep 'VYXE' nocaurban.theta.thetasWindow.pestPG >> nocaurban.thetasWindow.foroutliers

cp nocaurban.thetasWindow.foroutliers nocaurban.thetasWindow.forplot

sed -i 's/VYXE//g' nocaurban.thetasWindow.forplot
sed -i 's/\.1\t/\t/g' nocaurban.thetasWindow.forplot




# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('nocaurban.thetasWindow.foroutliers', sep ='\t')

df <- read.csv('nocaurban.thetasWindow.forplot', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )




nrow(df_outliers)
# 13723 rows
# 14 sig



max(df_outliers$Tajima)
# 2.679703
min(df_outliers$Tajima)
# -1.71403

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:14,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.71403
max(outlier_taj_disorder$Tajima)
# -1.271379
write_tsv(outlier_taj, "nocaurban.outliertaj_50kb.tsv")




png("nocaurban.tajima.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1806 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.271379) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> nocaurban.taj.relevantgenes_50kb.txt


done < nocaurban.outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' nocaurban.taj.relevantgenes_50kb.txt | sed 's/Name\=//g' | sort -u > nocaurban.taj.relevantgenenames_50kb.txt


# noca rural

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > nocarural.thetasWindow.foroutliers
grep 'VYXE' nocarural.theta.thetasWindow.pestPG >> nocarural.thetasWindow.foroutliers

cp nocarural.thetasWindow.foroutliers nocarural.thetasWindow.forplot

sed -i 's/VYXE//g' nocarural.thetasWindow.forplot
sed -i 's/\.1\t/\t/g' nocarural.thetasWindow.forplot




# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('nocarural.thetasWindow.foroutliers', sep ='\t')

df <- read.csv('nocarural.thetasWindow.forplot', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )




nrow(df_outliers)
# 13723 rows
# 14 sig



max(df_outliers$Tajima)
# 2.64898
min(df_outliers$Tajima)
# -1.8341

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:14,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.8341
max(outlier_taj_disorder$Tajima)
# -1.385304
write_tsv(outlier_taj, "nocarural.outliertaj_50kb.tsv")




png("nocarural.tajima.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1806 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.385304) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> nocarural.taj.relevantgenes_50kb.txt


done < nocarural.outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' nocarural.taj.relevantgenes_50kb.txt | sed 's/Name\=//g' | sort -u > nocarural.taj.relevantgenenames_50kb.txt























# PYRR
# pyrr urban

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > pyrrurban.thetasWindow.foroutliers
grep 'VYXE' pyrrurban.theta.thetasWindow.pestPG >> pyrrurban.thetasWindow.foroutliers

cp pyrrurban.thetasWindow.foroutliers pyrrurban.thetasWindow.forplot

sed -i 's/VYXE//g' pyrrurban.thetasWindow.forplot
sed -i 's/\.1\t/\t/g' pyrrurban.thetasWindow.forplot




# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('pyrrurban.thetasWindow.foroutliers', sep ='\t')

df <- read.csv('pyrrurban.thetasWindow.forplot', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )




nrow(df_outliers)
# 13723 rows
# 14 sig



max(df_outliers$Tajima)
# 2.771468
min(df_outliers$Tajima)
# -1.411386

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:14,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.411386
max(outlier_taj_disorder$Tajima)
# -1.175647
write_tsv(outlier_taj, "pyrrurban.outliertaj_50kb.tsv")




png("pyrrurban.tajima.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1806 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.175647) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> pyrrurban.taj.relevantgenes_50kb.txt


done < pyrrurban.outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' pyrrurban.taj.relevantgenes_50kb.txt | sed 's/Name\=//g' | sort -u > pyrrurban.taj.relevantgenenames_50kb.txt


# pyrr rural

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > pyrrrural.thetasWindow.foroutliers
grep 'VYXE' pyrrrural.theta.thetasWindow.pestPG >> pyrrrural.thetasWindow.foroutliers

cp pyrrrural.thetasWindow.foroutliers pyrrrural.thetasWindow.forplot

sed -i 's/VYXE//g' pyrrrural.thetasWindow.forplot
sed -i 's/\.1\t/\t/g' pyrrrural.thetasWindow.forplot




# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('pyrrrural.thetasWindow.foroutliers', sep ='\t')

df <- read.csv('pyrrrural.thetasWindow.forplot', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )




nrow(df_outliers)
# 13723 rows
# 14 sig



max(df_outliers$Tajima)
# 2.788987
min(df_outliers$Tajima)
# -1.512961

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:14,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.512961
max(outlier_taj_disorder$Tajima)
# -1.320125
write_tsv(outlier_taj, "pyrrrural.outliertaj_50kb.tsv")




png("pyrrrural.tajima.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 1806 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.320125) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 25000))
maxpos=$((midpos + 25000))

grep "$chr" /xdisk/mcnew/dannyjackson/cardinals/reference_dataset/ncbi_dataset/data/GCA_013397215.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> pyrrrural.taj.relevantgenes_50kb.txt


done < pyrrrural.outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[2])}' pyrrrural.taj.relevantgenes_50kb.txt | sed 's/Name\=//g' | sort -u > pyrrrural.taj.relevantgenenames_50kb.txt

