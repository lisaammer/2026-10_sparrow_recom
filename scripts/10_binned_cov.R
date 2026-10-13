#### Get average sequence coverage for binned chromosomes
#### Updated 13/10/2026 by Lisa Ammer

rm(list = ls())
gc()

# Load necessary libraries
library(dplyr)
library(tidyr)
library(GenomicRanges)

input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
cov_dir <- "/mnt/loki/johnston/lisa/sparrow/reference/processed_yapp_coverage"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"

# Import bin_lut from input directory
bin_lut <- read.table(file.path(input_dir, "bin_lut.txt"),
    header = TRUE,
    sep = "\t"
) %>%
    # Chromosome as character
    mutate(chromosome = as.character(chromosome))

# Import coverage data from input directory
coverage_data <- read.table(
    file.path(cov_dir, "compiled_coverage_data.txt"),
    header = TRUE,
    sep = "\t"
) %>%
    # Chromosome as character
    mutate(chromosome = as.character(chromosome))

# Create pair_id to identify parent_offspring pairs
coverage_data <- coverage_data %>%
    mutate(pair_id = paste0(parent, "_", offspring))

# Filter bin_lut to only include chromosomes present in coverage_data
bin_lut <- bin_lut %>%
    filter(chromosome %in% coverage_data$chromosome)

# Calculate per-chromosome denominators
pairs_per_chrom <- coverage_data %>%
    group_by(chromosome) %>%
    summarise(n_pairs_chrom = n_distinct(pair_id), .groups = "drop")

# Create GRanges for coverage
cov_gr <- GRanges(
    seqnames = coverage_data$chromosome,
    ranges = IRanges(start = coverage_data$cov_start, end = coverage_data$cov_end),
    pair_id = coverage_data$pair_id
)

# Collapse segments within pairs and get absolute depth by counting unique pair_ids covering each base
# reduce() ensures we don't double-count if a single pair has overlapping segments
cov_gr_reduced <- unlist(GenomicRanges::reduce(split(cov_gr, ~pair_id)))
cvg_depth <- coverage(cov_gr_reduced)

# Define bins
bin_gr <- GRanges(
    seqnames = bin_lut$chromosome,
    ranges = IRanges(bin_lut$bin_start, bin_lut$bin_end)
)

# Restrict bin_gr to only the chromosomes present in cvg_depth
seqlevels(bin_gr, pruning.mode = "tidy") <- names(cvg_depth)

# Calculate average absolute depth for each bin
# This gives the average number of pairs covering each base in the bin
bin_res <- binnedAverage(bin_gr, cvg_depth, "avg_depth")

# Normalise average depth by the total number of pairs for that chromosome to get average coverage per pair
final_df <- as.data.frame(bin_res) %>%
    mutate(chromosome = as.character(seqnames)) %>%
    left_join(pairs_per_chrom, by = "chromosome") %>%
    mutate(
        # average depth / total possible pairs
        avg_cov = avg_depth / n_pairs_chrom,
        # Remove NAs for chromosomes without coverage data
        avg_cov = replace_na(avg_cov, 0)
    ) %>%
    left_join(bin_lut %>% select(chromosome, bin_start, bin_end, bin),
        by = c("chromosome", "start" = "bin_start", "end" = "bin_end")
    ) %>%
    select(chromosome, bin, avg_cov) %>%
    arrange(chromosome, bin)

# Write output
write.table(final_df, file.path(output_dir, "binned_avg_coverage.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

rm(list = ls())
gc()

# End of script