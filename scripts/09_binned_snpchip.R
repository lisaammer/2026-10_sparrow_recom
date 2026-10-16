#### Get binned SNP chip distribution for chromosomes
#### Updated 13/10/2026 by Lisa Ammer

rm(list = ls()) 
gc()

# Load necessary libraries
library(dplyr)
library(tidyr)

# Set directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
reference_dir <- "/mnt/loki/johnston/lisa/sparrow/reference"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"

# Import bin_lut from input directory
bin_lut <- read.table(file.path(input_dir, "bin_lut.txt"),
    header = TRUE,
    sep = "\t"
)

# Import snp data from input directory
snp_data <- read.table(file.path(reference_dir, "70K_200K_maf_geno_mind_v6_newpos_SUPER_sorted.bim"),
    header = FALSE,
    sep = "\t"
)

# Manually assign header names
colnames(snp_data) <- c("chromosome", "snp_id", "morgan_position", "base_pair_position", "allele1", "allele2")

# Chromosome as caracter for merging
snp_data$chromosome <- as.character(snp_data$chromosome)
bin_lut$chromosome <- as.character(bin_lut$chromosome)

# Merge snp_data to bin_lut to assign bins to each SNP
# base_pair_position in snp_data should be between bin_start and bin_end in bin_lut for the same chromosome
binned_snp_data <- bin_lut %>%
    left_join(snp_data, by = "chromosome") %>%
    filter(base_pair_position >= bin_start & base_pair_position <= bin_end) %>%
    select(chromosome, bin, snp_id, morgan_position, base_pair_position, allele1, allele2)

# Calculate SNP distribution per bin by counting the number of SNPs in each bin and dividing by the bin size
snp_distribution <- binned_snp_data %>%
    group_by(chromosome, bin) %>%
    summarise(snp_count = n(), .groups = "drop") %>%
    left_join(bin_lut, by = c("chromosome", "bin")) %>%
    mutate(
        bin_size = bin_end - bin_start + 1, # Calculate bin size
        snp_distribution = snp_count / bin_size
    ) %>% # Calculate SNP distribution
    select(chromosome, bin, snp_distribution)

# Join back to bin_lut to get binned_snp_distribution for all bins, filling missing values with 0
binned_snp_distribution <- bin_lut %>%
    left_join(snp_distribution, by = c("chromosome", "bin")) %>%
    mutate(snp_distribution = ifelse(is.na(snp_distribution), 0, snp_distribution)) %>%
    select(chromosome, bin, snp_distribution)

# Write to output directory
write.table(binned_snp_distribution, file.path(output_dir, "binned_snp_distribution.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(list = ls())
gc()

# End of script