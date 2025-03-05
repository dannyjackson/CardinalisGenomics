# RAiSD.md 

# merge vcfs

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocaurban.txt

grep 'NOCA' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_nocarural.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Urban' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrurban.txt

grep 'PYRR' /xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset.txt | grep 'Rural' | awk '{print $1}' > /xdisk/mcnew/dannyjackson/cardinals/referencelists/samplenames_pyrrrural.txt




sbatch --account=mcnew \
    --job-name=raisd.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.raisd.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=10:00:00 \
    ~/programs/raisd_preprocessing.sh 


SPECIES=( "nocaurban" "nocarural" "pyrrurban" "pyrrrural" )
WINDOW_SIZES=( 10 20 50 100 500 1000 )

# Iterate over each combination
for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    sbatch --account=mcnew \
        --job-name=raisd_${WIN}_${SP} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.raisd_${WIN}_${SP}.%j \
        --nodes=1 \
        --time=1:00:00 \
        ~/programs/CardinalisGenomics/Genomics-Main/raisd/raisd.sh -p ~/programs/CardinalisGenomics/${SP}_params_raisd.sh -n ${SP} -w ${WIN}
    done 
done

~/programs/CardinalisGenomics/Genomics-Main/raisd/raisd.sh -p ~/programs/CardinalisGenomics/pyrrurban_params_raisd.sh -n pyrrurban -w 10000

Rscript ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/manhattanplot.r "/xdisk/mcnew/dannyjackson/cardinals" "#2d6a4f" "#74c69d" "0.001" "pyrrurban/pyrrurban.raisd_10000.Ztransformed.csv" "10000" "raisd" "pyrrurban"

# plot overlapping

# plot it in R

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/dannyjackson/cardinals"
color1 <- "#2d6a4f"
color2 <- "#74c69d"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "pyrrurban/pyrrurban.raisd_50.Ztransformed.csv"
input2 <- "pyrrrural/pyrrrural.raisd_50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "pyrrurban"
pop2 <- "pyrrrural"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = ",", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = ",", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)



# Prepare data for plotting
# data 1
cat("Preparing data for plotting...\n")
data1$chromo <- factor(data1$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data1 <- data1 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data1, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data1 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# data 2
cat("Preparing data for plotting...\n")
data2$chromo <- factor(data2$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data2 <- data2 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf2 <- plot_data2 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)

  geom_point(data = plot_data1, aes(x = BPcum, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = plot_data2, aes(x = BPcum, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data1$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "pyrr.raisd.50.pdf", 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")




# Plot just chromosome 2

# pyrr 
head -n 1 pyrrurban/pyrrurban.raisd_50.Ztransformed.csv > pyrrurban/pyrrurban.chr2.raisd_50.Ztransformed.csv
awk -F',' '$1 == "2"' pyrrurban/pyrrurban.raisd_50.Ztransformed.csv >> pyrrurban/pyrrurban.chr2.raisd_50.Ztransformed.csv
head -n 1 pyrrrural/pyrrrural.raisd_50.Ztransformed.csv > pyrrrural/pyrrrural.chr2.raisd_50.Ztransformed.csv
awk -F',' '$1 == "2"' pyrrrural/pyrrrural.raisd_50.Ztransformed.csv >> pyrrrural/pyrrrural.chr2.raisd_50.Ztransformed.csv



# plot it in R

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/dannyjackson/cardinals"
"0.01" "pyrrurban/pyrrurban.raisd_10000.Ztransformed.csv" "50" "raisd" "pyrrurban"
color1 <- "#2d6a4f"
color2 <- "#74c69d"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "pyrrurban/pyrrurban.chr2.raisd_50.Ztransformed.csv"
input2 <- "pyrrrural/pyrrrural.chr2.raisd_50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "pyrrurban"
pop2 <- "pyrrrural"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = ",", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = ",", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)




# Prepare data for plotting

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)
  geom_point(data = data1, aes(x = position, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = data2, aes(x = position, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  
    # Vertical dashed line at 57805281
  geom_vline(xintercept = 57805281, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 61205541, linetype = "dashed", color = "black") +

  
  # X-axis with chromosome labels and tick marks every 2M bp
  scale_x_continuous(
    labels = data1$chromo, 
    breaks = data1$center, 
    sec.axis = dup_axis(
      breaks = seq(0, 152000000, 
                   by = 3.8e+07), 
      labels = scales::comma
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = file.path("pyrr.chr2.raisd.pdf"), 
       width = 5, height = 5, units = "in")

cat("Script completed successfully!\n")


# noca

# noca 
head -n 1 nocaurban/nocaurban.raisd_50.Ztransformed.csv > nocaurban/nocaurban.chr2.raisd_50.Ztransformed.csv
awk -F',' '$1 == "2"' nocaurban/nocaurban.raisd_50.Ztransformed.csv >> nocaurban/nocaurban.chr2.raisd_50.Ztransformed.csv
head -n 1 nocarural/nocarural.raisd_50.Ztransformed.csv > nocarural/nocarural.chr2.raisd_50.Ztransformed.csv
awk -F',' '$1 == "2"' nocarural/nocarural.raisd_50.Ztransformed.csv >> nocarural/nocarural.chr2.raisd_50.Ztransformed.csv



# plot it in R

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/dannyjackson/cardinals"
"0.01" "nocaurban/nocaurban.raisd_10000.Ztransformed.csv" "50" "raisd" "nocaurban"
color1 <- "#6247aa"
color2 <- "#b185db"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "nocaurban/nocaurban.chr2.raisd_50.Ztransformed.csv"
input2 <- "nocarural/nocarural.chr2.raisd_50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "nocaurban"
pop2 <- "nocarural"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = ",", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = ",", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)




# Prepare data for plotting

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)
  geom_point(data = data1, aes(x = position, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = data2, aes(x = position, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  

  
  # X-axis with chromosome labels and tick marks every 2M bp
  scale_x_continuous(
    labels = data1$chromo, 
    breaks = data1$center, 
    sec.axis = dup_axis(
      breaks = seq(0, 152000000, 
                   by = 3.8e+07), 
      labels = scales::comma
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = file.path("noca.chr2.raisd.pdf"), 
       width = 5, height = 5, units = "in")

cat("Script completed successfully!\n")
