#### Get GC for binned chromosomes
#### Created 06/02/2026 by Lisa Ammer
#### Ran on ac3

rm(list = ls())
gc()

# Load necessary libraries
library(dplyr)
library(tidyr)
library(Biostrings)

# Get input directory from working directory/_in
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
# Get output direktory from working directory/_in
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Import bin_lut from input directory
bin_lut <- read.table(file.path(input_dir, "bin_lut.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))

# Load the sorted FASTA file
fasta_file <- file.path(input_dir, "filtered_fasta.fasta")
fasta <- Biostrings::readDNAStringSet(fasta_file, format = "FASTA")

# Initialize the container
gc_results_list <- list()

# Loop through each chromosome in the bin_lut and calculate GC content for each bin
for (chr in unique(bin_lut$chromosome)) {
    if (!chr %in% names(fasta)) {
        warning(paste("Chromosome", chr, "not found in fasta"))
        next
    }

    seq <- fasta[[chr]]
    seq_len <- length(seq)

    bins <- bin_lut[bin_lut$chromosome == chr, ]

    # Clip bin_end so it never exceeds chromosome length and filter out invalid starts
    bins <- bins %>%
        mutate(bin_end = pmin(bin_end, seq_len)) %>%
        filter(bin_start <= seq_len & bin_start > 0)

    if (nrow(bins) == 0) next

    # Calculate GC for all bins in this chromosome at once
    chr_views <- Biostrings::Views(seq, start = bins$bin_start, end = bins$bin_end)

    # Calculate frequency for all windows in one command
    freqs <- Biostrings::letterFrequency(chr_views, letters = c("G", "C"), as.prob = FALSE)

    # Calculate proportions (rowSums adds G + C)
    bin_gc_prop <- rowSums(freqs) / width(chr_views)

    # Store in the list correctly
    gc_results_list[[chr]] <- data.frame(
        chromosome = chr,
        bin = bins$bin,
        gc_prop = bin_gc_prop
    )
}

# Combine everything at the end
gc_prop_per_bin <- do.call(rbind, gc_results_list)

# Write the result to the output directory
write.table(gc_prop_per_bin, file.path(output_dir, "binned_gc_prop.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

message(">>> GC proportion per bin calculated successfully.")

rm(list = ls())
gc()

# End of script