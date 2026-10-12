#### Get TSS for binned chromosomes
#### Differentiate between CGI overlap and no CGI overlap
#### Updated 12/10/2026 by Lisa Ammer

rm(list = ls())
gc()

# Load necessary libraries
library(dplyr)
library(tidyr)
library(GenomicRanges)

# Set Directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# PART 1: MAP GENES AND TSS TO BINS
# =============================================================================

message(">>> Mapping genes and TSS to genomic bins...")

# Import bin_lut from input directory
bin_lut <- read.table(file.path(input_dir, "bin_lut.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))

# Import filtered gff from input directory
gff_filtered_df <- read.table(file.path(input_dir, "filtered_gff.txt"), header = TRUE, sep = "\t", quote = "")

# Select columns and define all TSS based on strand orientation
genes <- gff_filtered_df %>%
    select(chromosome, start, end, genomic_length, exonic_length, gene, type, gene_biotype, strand) %>%
    mutate(
        chromosome = as.character(chromosome),
        genomic_length = as.numeric(genomic_length),
        exonic_length = as.numeric(exonic_length),
        # For negative strands, TSS is the 'end' (highest coordinate)
        # For positive strands, TSS is the 'start' (lowest coordinate)
        TSS = if_else(strand == "-", end, start)
    )

# Add bin information based on the TSS column using non-equi join
genes <- genes %>%
    left_join(bin_lut,
        by = join_by(
            chromosome == chromosome,
            TSS >= bin_start,
            TSS <= bin_end
        )
    ) %>%
    dplyr::select(chromosome, bin, gene, type, gene_biotype, strand, start, end, TSS, genomic_length, exonic_length)

# Ensure only one entry per unique "gene" and "chromosome" and "bin" combination, if multiple entries exist, keep the the one with gene_biotype information (non-NA), otherwise keep the first entry
genes_final <- genes %>%
    group_by(chromosome, bin, gene) %>%
    slice_min(order_by = is.na(gene_biotype), n = 1)

# Write the master gene list with bins to output directory
write.table(genes_final, file.path(input_dir, "all_genes_with_bins.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

message(">>> Gene-to-bin mapping complete and saved.")

# =============================================================================
# PART 2: INTEGRATE CGI OVERLAPS & SUMMARISE TSS
# =============================================================================

message(">>> Integrating CpG islands and calculating TSS summary...")

# Import CGI map (generated from your BED file)
cgi_map <- read.table(file.path(input_dir, "cpg_islands_map_from_bed.txt"), header = TRUE, sep = "\t") %>%
    filter(!is.na(chromosome)) %>%
    mutate(chromosome = as.character(chromosome))

# Convert to GRanges
# Note: use width = 1 to represent the exact TSS point
tss_gr <- GRanges(
    seqnames = genes_final$chromosome,
    ranges = IRanges(start = genes_final$TSS, width = 1),
    gene = genes_final$gene
)

cgi_gr <- GRanges(
    seqnames = cgi_map$chromosome,
    ranges = IRanges(start = cgi_map$start_bp, end = cgi_map$end_bp)
)

# Identify which TSS fall within a CGI
overlaps <- findOverlaps(tss_gr, cgi_gr)

# Create a logical vector for overlap status in the genes dataframe
genes_final$has_cgi <- FALSE
genes_final$has_cgi[queryHits(overlaps)] <- TRUE

# Summarise per bin (Total TSS, with CGI, without CGI - without biotype split)
tss_cgi_summary <- genes_final %>%
    group_by(chromosome, bin) %>%
    summarise(
        total_tss = n(),
        tss_with_cgi = sum(has_cgi == TRUE),
        tss_no_cgi = sum(has_cgi == FALSE),
        .groups = "drop"
    ) %>%
    # Join with bin_lut to fill in gaps for empty bins
    right_join(
        bin_lut %>% select(chromosome, bin),
        by = c("chromosome", "bin")
    ) %>%
    # Replace resulting NAs with 0
    mutate(across(c(total_tss, tss_with_cgi, tss_no_cgi), ~ replace_na(., 0))) %>%
    arrange(as.character(chromosome), as.integer(bin))

# Export results
write.table(genes_final,
    file.path(input_dir, "all_genes_with_bins_and_cgi_status.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(tss_cgi_summary,
    file.path(output_dir, "binned_tss_and_cgi.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

message(">>> Pipeline complete! All gene and TSS CGI summaries saved.")

rm(list = ls())
gc()

# End of script