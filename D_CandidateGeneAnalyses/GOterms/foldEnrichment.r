# Compare GO enrichment between species

# Load required libraries
library(tidyverse)
library(ggplot2)

# Function to convert "<X" values to numeric
convert_numeric <- function(x) {
  x <- as.character(x)  # Ensure it's character
  x <- gsub("<", "", x) # Remove "<" symbols
  as.numeric(x)         # Convert to numeric
}

# Read in PANTHER GO enrichment results for both species
noca.fstraisd <- read.csv("noca.fstraisd.50kb.txt", sep="\t", header=TRUE)
noca.fst <- read.csv("noca.fst.50kb.txt", sep="\t", header=TRUE)
noca.raisd <- read.csv("noca.raisd.50kb.txt", sep="\t", header=TRUE)
noca <- full_join(noca.fstraisd, noca.fst, noca.raisd, by = "PANTHER.GO.Slim.Biological.Process")


pyrr.fstraisd <- read.csv("pyrr.fstraisd.50kb.txt", sep="\t", header=TRUE)
pyrr.fst <- read.csv("pyrr.fst.50kb.txt", sep="\t", header=TRUE)
pyrr.raisd <- read.csv("pyrr.raisd.50kb.txt", sep="\t", header=TRUE)
pyrr <- full_join(pyrr.fstraisd, pyrr.fst, pyrr.raisd, by = "PANTHER.GO.Slim.Biological.Process")

noca.fstraisd.fdr <- read.csv("noca.fstraisd.50kb.fdr.txt", sep="\t", header=TRUE)
pyrr.fstraisd.fdr <- read.csv("pyrr.fstraisd.50kb.fdr.txt", sep="\t", header=TRUE)
noca.fst.fdr <- read.csv("noca.fst.50kb.fdr.txt", sep=",", header=TRUE)
pyrr.fst.fdr <- read.csv("pyrr.fst.50kb.fdr.txt", sep=",", header=TRUE)
noca.raisd.fdr <- read.csv("noca.raisd.50kb.fdr.txt", sep=",", header=TRUE)
pyrr.raisd.fdr <- read.csv("pyrr.raisd.50kb.fdr.txt", sep=",", header=TRUE)

# Ensure column names are correct (modify if needed)
colnames(noca) <- c("GO_term", "Background", "Observed", "Expected", "Direction", "Fold_Enrichment", "P_value")
colnames(pyrr) <- c("GO_term", "Background", "Observed", "Expected", "Direction", "Fold_Enrichment", "P_value")

# Convert problematic columns to numeric
noca$Fold_Enrichment <- convert_numeric(noca$Fold_Enrichment)
noca$P_value <- convert_numeric(noca$P_value)
pyrr$Fold_Enrichment <- convert_numeric(pyrr$Fold_Enrichment)
pyrr$P_value <- convert_numeric(pyrr$P_value)

# Convert observed and background counts to integers
noca$Observed <- as.integer(noca$Observed)
noca$Background <- as.integer(noca$Background)
pyrr$Observed <- as.integer(pyrr$Observed)
pyrr$Background <- as.integer(pyrr$Background)

# Merge datasets on GO term
noca2 <- noca %>% janitor::remove_empty("cols")
pyrr2 <- pyrr %>% janitor::remove_empty("cols")
merged_data_temp <- full_join(noca2, pyrr2, by="GO_term", suffix=c("_noca", "_pyrr"))

# Extract unique GO terms
# Extract the GO terms from the column
noca.fstraisd.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", noca.fstraisd.fdr$PANTHER.GO.Slim.Biological.Process)
pyrr.fstraisd.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", pyrr.fstraisd.fdr$PANTHER.GO.Slim.Biological.Process)
noca.fst.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", noca.fst.fdr$PANTHER.GO.Slim.Biological.Process)
pyrr.fst.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", pyrr.fst.fdr$PANTHER.GO.Slim.Biological.Process)
noca.raisd.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", noca.raisd.fdr$PANTHER.GO.Slim.Biological.Process)
pyrr.raisd.fdr$GO_term <- gsub(".*\\((GO:[0-9]+\\.[0-9]+)\\)", "\\1", pyrr.raisd.fdr$PANTHER.GO.Slim.Biological.Process)

unique_go_terms <- unique(c(noca.fstraisd.fdr$GO_term, pyrr.fstraisd.fdr$GO_term, noca.fst.fdr$GO_term, pyrr.fst.fdr$GO_term, noca.raisd.fdr$GO_term, pyrr.raisd.fdr$GO_term))
merged_data <- merged_data_temp %>% filter(GO_term %in% unique_go_terms)

# Fisher's exact test function
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    fisher_p = ifelse(
      all(!is.na(c(Observed_noca, Background_noca, Observed_pyrr, Background_pyrr))), 
      fisher.test(matrix(c(Observed_noca, Background_noca - Observed_noca, 
                           Observed_pyrr, Background_pyrr - Observed_pyrr), 
                         nrow=2))$p.value, 
      NA_real_ # Assign NA if missing data
    )
  ) %>%
  ungroup()  # Ensure it's no longer rowwise

# Calculate log2 fold change in enrichment (handle zeros)
merged_data <- merged_data %>%
  mutate(log2_fold_change = log2((Fold_Enrichment_noca + 1e-6) / (Fold_Enrichment_pyrr + 1e-6))) # Avoid log(0)

# Adjust p-values using FDR correction
merged_data <- merged_data %>%
  mutate(fdr_p = p.adjust(fisher_p, method="fdr"))

# Check if fisher_p column is correctly generated
print(merged_data %>% select(GO_term, fisher_p, fdr_p, log2_fold_change) %>% head())

# Define significance threshold
p_threshold <- 0.2

# Subset significant GO terms
significant_results <- merged_data %>%
  filter(fdr_p < p_threshold)

# Save significant results to a CSV file
write.csv(significant_results, "GO_significant_results.csv", row.names=FALSE)

# Print first few rows to verify
print(head(significant_results))

# Plot enrichment comparison
ggplot(merged_data, aes(x=reorder(GO_term, log2_fold_change), y=log2_fold_change, fill=log2_fold_change > 0)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("FALSE" = "#2d6a4f", "TRUE" = "#7570B3")) +  # Custom color scheme
  labs(x="GO Term", y="Log2 Fold Change (NOCA vs. PYRR)", title="Functional Enrichment Comparison") +
  theme_minimal()

ggsave(
  "enrichmentcomparison.pdf",
  plot = last_plot(),
  width = 10,
  height = 15,
  units = "in",
  dpi = 300
)

# Save results
write.csv(merged_data, "GO_comparison_results.csv", row.names=FALSE)

