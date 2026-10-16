#### Sparrows Chromosome Binning, TSS Extraction & FASTA Prep
#### Updated 09/10/2026 by Lisa Ammer

rm(list = ls()) # Clear workspace
gc()

# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)
library(fs)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(purrr)
library(tidyverse)
library(jsonlite)

# Set Directories
input_dir  <- "/mnt/loki/johnston/lisa/sparrow/reference/hap1"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed" 

# Target Chromosome Selection (1:37, Z, W)
target_chroms <- c(as.character(1:37), "Z", "W")

# ==================================================
# PART 1: Process Sequence Report & Create Bin Lookup
# ==================================================

# Read sequence report
sequence_report <- jsonlite::stream_in(
    file(file.path(input_dir, "sequence_report.jsonl")),
    verbose = FALSE
)

# Filter sequence report for target primary chromosomes
filtered_report <- sequence_report %>%
    rename_with(trimws) %>%
    filter(as.character(chrName) %in% target_chroms) %>%
    filter(!grepl("unloc", sequenceName, ignore.case = TRUE))

# Write filtered sequence report
write.table(filtered_report, file.path(output_dir, "filtered_sequence_report.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Generate 100kb Bin Lookup Table
bin_size <- 100000

bin_lut <- filtered_report %>%
  rowwise() %>%
  do({
    chr_length <- .$length
    chr_name <- .$chrName
    num_bins   <- ceiling(chr_length / bin_size)
    
    data.frame(
      chromosome        = chr_name,
      bin               = 1:num_bins,
      bin_start         = seq(1, by = bin_size, length.out = num_bins),
      bin_end           = pmin(seq(bin_size, by = bin_size, length.out = num_bins), chr_length),
      chromosome_length = chr_length
    ) %>%
      mutate(bin_mid = (bin_start + bin_end) / 2)
  }) %>%
  dplyr::select(chromosome, bin, bin_start, bin_mid, bin_end, chromosome_length) %>%
  bind_rows()

write.table(bin_lut, file.path(output_dir, "bin_lut.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==================================================
# PART 2: Process GFF, Extract TSS & Exonic Lengths
# ==================================================

# Read filtered sequence report
report <- read.table(file.path(output_dir, "filtered_sequence_report.txt"), header = TRUE, sep = "\t") %>%
  mutate(chrName = as.character(chrName))

# Import GFF3
gff <- import(file.path(input_dir, "genomic.gff"), format = "gff3")

# Convert GRanges to dataframe
gff_df <- as.data.frame(gff) %>%
  mutate(seqnames = as.character(seqnames)) %>%
    dplyr::select(-any_of("chromosome")) %>% # Remove conflicting old chromosome column
    mutate(across(where(is.list), ~ map_chr(., ~ paste(unique(.x), collapse = ","))))

lookup <- report[, c("refseqAccession", "chrName")]

# 4. Perform join, rename, and filter
gff_mapped <- gff_df %>%
  inner_join(lookup, by = c("seqnames" = "refseqAccession")) %>%
  dplyr::rename(chromosome = chrName) %>%
  filter(chromosome %in% target_chroms) %>%
  filter(chromosome != "Unknown")

# 5. Verify the mapped output
head(gff_mapped[, c("seqnames", "chromosome", "start", "end", "type")])
head(gff_mapped[gff_mapped$chromosome == "4", c("seqnames", "chromosome", "start", "end", "type")])
head(gff_mapped[gff_mapped$chromosome == "37", c("seqnames", "chromosome", "start", "end", "type")])

# --- Section A: Extract Annotated TSS (All Transcripts) ---
tss_df <- gff_mapped %>%
  filter(!type %in% c("gene", "pseudogene", "exon", "CDS", "region", "cDNA_match", "match")) %>%
  filter(!is.na(Parent) | !is.na(transcript_id) | str_detect(type, "RNA|transcript|gene_segment")) %>%
  mutate(
    tss_pos = ifelse(strand == "+", start, end),
    bed_start = tss_pos - 1,
    bed_end   = tss_pos,
    feature_label = ifelse(!is.na(gene) & gene != "", gene, 
                    ifelse(!is.na(transcript_id) & transcript_id != "", transcript_id, ID))
  ) %>%
  dplyr::select(chromosome, bed_start, bed_end, feature_label, score, strand, ID, type) %>%
  arrange(suppressWarnings(as.numeric(chromosome)), bed_start)

tss_df %>% filter(chromosome == "4") %>% head()
tss_df %>% filter(chromosome == "37") %>% head()

write.table(tss_df, file.path(output_dir, "extracted_tss_target_chroms.bed"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- Section B: Exonic Length Calculation ---
# Subset exons ONLY for target chromosomes
mapped_seqnames <- unique(gff_mapped$seqnames)
exons_target <- gff[gff$type == "exon" & as.character(seqnames(gff)) %in% mapped_seqnames]

exon_mcols <- mcols(exons_target)
group_ids  <- ifelse(!is.na(exon_mcols$gene) & exon_mcols$gene != "", exon_mcols$gene,
               ifelse(!is.na(exon_mcols$Parent), as.character(exon_mcols$Parent), as.character(exon_mcols$ID)))

exons_by_feature <- split(exons_target, group_ids)
reduced_exons    <- GenomicRanges::reduce(exons_by_feature)

exonic_lengths_df <- data.frame(
  feature_symbol = names(sum(width(reduced_exons))),
  exonic_length  = as.numeric(sum(width(reduced_exons)))
)

target_feature_types <- c("gene", "pseudogene", "mRNA", "lnc_RNA", "transcript", "rRNA", "snRNA", "snoRNA", "tRNA")

gff_filtered_df <- gff_mapped %>%
  filter(type %in% target_feature_types) %>%
  mutate(
    join_key = ifelse(!is.na(gene) & gene != "", gene,
               ifelse(!is.na(Parent), as.character(Parent), as.character(ID)))
  ) %>%
  left_join(exonic_lengths_df, by = c("join_key" = "feature_symbol")) %>%
  mutate(genomic_length = end - start + 1) %>%
  dplyr::select(chromosome, start, end, strand, gene, ID, type, gene_biotype, description, exonic_length, genomic_length)

write.table(gff_filtered_df, file.path(output_dir, "filtered_gff.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


# ==================================================
# PART 3: Process FASTA Sequences
# ==================================================
fasta_file <- "GCF_036417665.1_bPasDom1.hap1_genomic.fna"
fasta <- readDNAStringSet(file.path(input_dir, fasta_file), format = "fasta")

# Clean FASTA headers
names(fasta) <- sub(" .*", "", names(fasta))

# Map and filter target chromosomes
fasta_mapped <- fasta[names(fasta) %in% report$refseqAccession]
names(fasta_mapped) <- report$chrName[match(names(fasta_mapped), report$refseqAccession)]

sorted_fasta <- fasta_mapped[names(fasta_mapped) %in% target_chroms]
sorted_fasta <- sorted_fasta[match(intersect(target_chroms, names(sorted_fasta)), names(sorted_fasta))]

# Save output FASTA
writeXStringSet(sorted_fasta, file.path(output_dir, "filtered_fasta.fasta"), format = "fasta")

rm(list = ls())
gc()

# End of script