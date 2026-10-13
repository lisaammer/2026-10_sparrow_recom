#### Get integrated repeats for binned chromosomes
#### Broadclass summary only
#### Updated 12/10/2026 by Lisa Ammer

rm(list = ls())
gc()

# Load necessary libraries
library(dplyr)
library(plyranges)
library(data.table)
library(GenomicRanges)

# Input and output directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"

# Import bin_lut from input directory
bin_lut <- read.table(file.path(input_dir, "bin_lut.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))

# Import repeat positions from input directory
repeats <- data.table::fread(file.path(input_dir, "master_repeat_list.txt"),
    sep = "\t",
    header = TRUE
) %>%
    as.data.frame() %>%
    mutate(
        chromosome = as.character(chromosome),
        broad_class = as.character(broad_class)
    ) %>%
    # Explicitly replace NAs and blanks, and drop any rows where broad_class is completely missing
    mutate(
        broad_class = tidyr::replace_na(broad_class, "Other/Unknown Repeat"),
        broad_class = if_else(broad_class == "" | broad_class == "NA", "Other/Unknown Repeat", broad_class)
    ) %>%
    filter(!is.na(broad_class))

# Only keep chromosomes present in both datasets
bin_lut <- bin_lut %>%
    filter(chromosome %in% repeats$chromosome)
repeats <- repeats %>%
    filter(chromosome %in% bin_lut$chromosome)

# Convert to GRanges
bins_gr <- GRanges(
    seqnames = bin_lut$chromosome,
    ranges = IRanges(start = bin_lut$bin_start, end = bin_lut$bin_end),
    bin = bin_lut$bin
)

repeats_gr <- GRanges(
    seqnames = repeats$chromosome,
    ranges = IRanges(start = repeats$start, end = repeats$end),
    broad_class = repeats$broad_class
)

message(">>> Reducing overlapping repeats per broad class...")

# Reduce ranges per broad class natively to prevent double-counting overlapping repeats
broad_reduced_list <- lapply(unique(repeats_gr$broad_class), function(bc) {
    sub_gr <- repeats_gr[repeats_gr$broad_class == bc]
    red_sub <- reduce(sub_gr)
    mcols(red_sub)$broad_class <- bc
    return(red_sub)
})
broad_reduced <- do.call(c, broad_reduced_list)

message(">>> Calculating overlaps with genomic bins...")

# Find overlaps between bins and reduced repeats
overlaps <- findOverlaps(bins_gr, broad_reduced)
intersect_gr <- pintersect(bins_gr[queryHits(overlaps)], broad_reduced[subjectHits(overlaps)])

# Build overlap summary table
overlap_df <- data.frame(
    chromosome = as.character(seqnames(bins_gr[queryHits(overlaps)])),
    bin = mcols(bins_gr)$bin[queryHits(overlaps)],
    broad_class = mcols(broad_reduced)$broad_class[subjectHits(overlaps)],
    overlap_width = width(intersect_gr)
)

# Summarise base pairs per broad class per bin
bin_repeat_summary <- overlap_df %>%
    group_by(chromosome, bin, broad_class) %>%
    summarise(total_bp = sum(overlap_width), .groups = "drop")

# Calculate proportions relative to bin size and pivot wider
broad_class_summary <- bin_lut %>%
    mutate(bin_width = bin_end - bin_start + 1) %>%
    left_join(bin_repeat_summary, by = c("chromosome", "bin")) %>%
    mutate(
        total_bp = tidyr::replace_na(total_bp, 0),
        prop = total_bp / bin_width
    ) %>%
   tidyr::pivot_wider(
       id_cols = c(chromosome, bin),
       names_from = broad_class,
       values_from = prop,
       values_fill = 0
   ) %>%
       # Explicitly remove the NA column if it managed to creep in
       dplyr::select(-any_of("NA"))

# ---------------------------------------------------------
# EXPORT RESULTS
# ---------------------------------------------------------

write.table(broad_class_summary, file.path(output_dir, "binned_repeat_prop.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

message(">>> Binned repeat proportions successfully calculated and saved!")

rm(list = ls())
gc()

# End of script