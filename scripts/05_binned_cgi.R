#### Get CPG islands from BED and calculate background GC
#### Created 14/10/2026 by Lisa Ammer
#### Updated to use Pdom_cpgIslands_Unmasked.bed as source of CpG islands

rm(list = ls())
gc()

library(dplyr)
library(Biostrings)
library(GenomicRanges)

# Paths
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Import bin_lut
bin_lut <- read.table(file.path(input_dir, "changed/bin_lut.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))
bins_gr <- GRanges(
    seqnames = bin_lut$chromosome,
    ranges = IRanges(start = bin_lut$bin_start, end = bin_lut$bin_end)
)

# Import filtered sequence report
report <- read.table(file.path(input_dir, "changed/filtered_sequence_report.txt"), header = TRUE, sep = "\t") %>%
    mutate(
        chromosome = as.character(chrName),
        refseqAccession = as.character(refseqAccession)
    )

# Load FASTA
fasta_file <- file.path(input_dir, "changed/filtered_fasta.fasta")
fasta <- Biostrings::readDNAStringSet(fasta_file, format = "FASTA")
names(fasta) <- trimws(sub(" .*", "", names(fasta)))

# --- Load and Map CpG Islands from BED File ---
message("Loading and mapping CpG islands from BED file...")
bed_file <- file.path(input_dir, "Pdom_cpgIslands_Unmasked.bed")

bed_df <- read.table(bed_file, header = FALSE, sep = "\t")

# Map BED sequence accessions (V1) to bin_lut chromosome names using filtered_sequence_report
bed_mapped <- bed_df %>%
    left_join(filtered_sequence_report %>% select(refseqAccession, chromosome),
        by = c("V1" = "refseqAccession")
    ) %>%
    filter(!is.na(chromosome)) # Drop unplaced scaffolds if not in bin_lut

# BED format is 0-based start, 1-based end (half-open -> add 1 to start)
all_cgi_gr <- GRanges(
    seqnames = bed_mapped$chromosome,
    ranges = IRanges(start = bed_mapped$V2 + 1, end = bed_mapped$V3)
)

message("Loaded and mapped ", length(all_cgi_gr), " CpG islands from BED file.")

# --- Generate CpG Island Map for Output ---
cgi_map <- as.data.frame(all_cgi_gr) %>%
    dplyr::select(seqnames, start, end, width) %>%
    dplyr::rename(chromosome = seqnames, start_bp = start, end_bp = end)

cgi_starts_gr <- GRanges(
    seqnames = cgi_map$chromosome,
    ranges = IRanges(start = cgi_map$start_bp, width = 1)
)

bin_overlaps <- findOverlaps(cgi_starts_gr, bins_gr, select = "first")
cgi_map$bin_id <- bin_lut$bin[bin_overlaps]

cgi_map <- cgi_map %>%
    select(chromosome, bin_id, start_bp, end_bp, width)
    dplyr::rename(bin = bin_id)

write.table(cgi_map, file.path(input_dir, "changed/cpg_islands_map_from_bed.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# --- Calculate Stats per Bin ---

# 1. CGI proportion (coverage)
cgi_overlaps <- findOverlaps(bins_gr, all_cgi_gr)
cgi_intersect <- pintersect(bins_gr[queryHits(cgi_overlaps)], all_cgi_gr[subjectHits(cgi_overlaps)])

cgi_width_per_bin <- aggregate(width(cgi_intersect), list(queryHits(cgi_overlaps)), sum)

bin_lut$cgi_prop <- 0
bin_lut$cgi_prop[cgi_width_per_bin$Group.1] <- cgi_width_per_bin$x / width(bins_gr)[cgi_width_per_bin$Group.1]

# 2. Background GC content (excluding BED CpG island regions)
message("Calculating Background GC content...")

bin_lut$bg_gc_prop <- NA

for (chr in unique(bin_lut$chromosome)) {
    if (!chr %in% names(fasta)) next

    chr_seq <- fasta[[chr]]
    chr_bins <- which(bin_lut$chromosome == chr)

    # Mask the CpG islands from the BED file by replacing them with 'N'
    masked_seq <- chr_seq
    chr_cgis_in_chr <- all_cgi_gr[seqnames(all_cgi_gr) == chr]

    if (length(chr_cgis_in_chr) > 0) {
        masked_seq <- replaceAt(chr_seq, ranges(chr_cgis_in_chr), DNAStringSet(rep("N", length(chr_cgis_in_chr)), use.names = FALSE))
    }

    # Create views for the bins on the MASKED sequence
    chr_views <- Views(masked_seq, start = bin_lut$bin_start[chr_bins], end = bin_lut$bin_end[chr_bins])

    # Count G/C and total valid bases (excluding 'N')
    bg_freqs <- letterFrequency(chr_views, letters = c("G", "C", "A", "T"), as.prob = FALSE)

    gc_counts <- rowSums(bg_freqs[, c("G", "C")])
    atgc_counts <- rowSums(bg_freqs)

    bin_lut$bg_gc_prop[chr_bins] <- ifelse(atgc_counts > 0, gc_counts / atgc_counts, 0)
}

# Select final columns for export
final_output <- bin_lut %>%
    dplyr::select(chromosome, bin, cgi_prop, bg_gc_prop)

# Export
write.table(final_output, file.path(output_dir, "binned_background_gc_and_cgi.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

message("Processed and mapped using Pdom_cpgIslands_Unmasked.bed.")

rm(list = ls())
gc()

# End of script