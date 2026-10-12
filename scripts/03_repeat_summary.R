#### Integrates TRASH and RepeatMasker TE annotations, generates summary statistics, and produces chromosome and Manhattan plots.
#### Updated 12/10/2026 by Lisa Ammer
#### Ran this on ac3

rm(list = ls())
gc()

# --- Load Necessary Libraries ---
library(tidyverse)
library(GenomicRanges)
library(data.table)
library(stringr)
library(ggplot2)
library(patchwork)

# --- Define Directories ---
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"

# Create specific subdirectories for outputs
plot_dir <- file.path("/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/03_repeat_chromosome_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# =============================================================================
# PART 1: INTEGRATE TRASH AND REPEATMASKER ANNOTATIONS
# =============================================================================

message(">>> Starting TRASH and RepeatMasker integration...")

# Import filtered_sequence_report
filtered_sequence_report <- read.table(file.path(input_dir, "filtered_sequence_report.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chrName))

# Import TRASH summary file
trash_repeats <- read.csv("/mnt/loki/johnston/lisa/sparrow/TRASH/TRASH_results/TRASH_hap1_win4000/all.repeats.from.GCF_036417665.1_bPasDom1.hap1_genomic.fna.csv") %>%
    left_join(filtered_sequence_report %>% select(refseqAccession, chromosome),
        by = c("seq.name" = "refseqAccession")
    ) %>%
    filter(!is.na(chromosome) & chromosome != "Un")

# Define standard headers for RepeatMasker output
rm_headers <- c(
    "SW_score", "perc_div", "perc_del", "perc_ins", "query_id",
    "query_start", "query_end", "query_left", "strand",
    "repeat_id", "class_family", "repeat_start", "repeat_end",
    "repeat_left", "ID", "extra"
)

# Import RepeatMasker summary file
te_summary <- read_tsv("/mnt/loki/johnston/lisa/sparrow/TEtools/Pas.dom1TE.tsv", col_names = rm_headers) %>%
    separate(class_family, into = c("class", "family"), sep = "/", fill = "right") %>%
    left_join(filtered_sequence_report %>% select(refseqAccession, chromosome),
        by = c("query_id" = "refseqAccession")
    ) %>%
    filter(!is.na(chromosome) & chromosome != "Un")

# TRASH GRanges
trash_gr <- GRanges(
    seqnames = trash_repeats$chromosome,
    ranges = IRanges(start = trash_repeats$start, end = trash_repeats$end),
    family = trash_repeats$seq
)

trash_gr_clean <- unlist(reduce(split(trash_gr, mcols(trash_gr)$family)))
mcols(trash_gr_clean)$family <- names(trash_gr_clean)
mcols(trash_gr_clean)$source <- "TRASH"
mcols(trash_gr_clean)$class <- "Satellite_Tandem"

# RepeatMasker GRanges
te_gr <- GRanges(
    seqnames = te_summary$chromosome,
    ranges = IRanges(start = te_summary$query_start, end = te_summary$query_end),
    class = te_summary$class,
    family = te_summary$family
)
mcols(te_gr)$source <- "RepeatMasker"

# Identify overlaps
logical_hits <- findOverlaps(trash_gr_clean, te_gr, select = "first")

# Hybrid regions
overlaps_gr <- trash_gr_clean[!is.na(logical_hits)]
matched_rm <- te_gr[logical_hits[!is.na(logical_hits)]]

mcols(overlaps_gr)$class <- mcols(matched_rm)$class
mcols(overlaps_gr)$family <- paste0(
    mcols(matched_rm)$family,
    " (TRASH: ", mcols(overlaps_gr)$family, ")"
)
mcols(overlaps_gr)$source <- "Hybrid"
mcols(overlaps_gr)$note <- "Confirmed_by_both"

# Unique TRASH
unique_trash_gr <- trash_gr_clean[is.na(logical_hits)]
mcols(unique_trash_gr)$note <- "TRASH_only"

# Unique RM
rm_hits_in_trash <- findOverlaps(te_gr, trash_gr_clean)
unique_rm_gr <- te_gr[-queryHits(rm_hits_in_trash)]
mcols(unique_rm_gr)$note <- "RM_only"

names(unique_trash_gr) <- NULL
names(unique_rm_gr) <- NULL
names(overlaps_gr) <- NULL

# Build master_repeats dataframe
master_repeats <- bind_rows(
    as.data.frame(unique_trash_gr),
    as.data.frame(unique_rm_gr),
    as.data.frame(overlaps_gr)
) %>%
    dplyr::select(chromosome = seqnames, start, end, class, family, source, note) %>%
    arrange(chromosome, start)

# Clean Hybrid strings and categories
master_repeats <- master_repeats %>%
    mutate(
        trash_sequence = str_extract(family, "(?<=TRASH: ).*?(?=\\))"),
        family_clean = str_remove(family, " \\(TRASH: .*?\\)"),
        family_clean = ifelse(family_clean == "NA" | is.na(family_clean), NA_character_, family_clean)
    ) %>%
    mutate(
        trash_sequence = ifelse(note == "TRASH_only", family, trash_sequence),
        family = ifelse(note == "TRASH_only", NA_character_, family_clean)
    ) %>%
    select(chromosome, start, end, class, family, trash_sequence, source, note) %>%
    mutate(class = str_c(class, "_Repeat")) %>%
    mutate(class = case_when(
        class == "Low_complexity_Repeat" ~ "Low_Complexity_Repeat",
        class == "Simple_repeat_Repeat" ~ "Simple_Repeat",
        TRUE ~ class
    )) %>%
    mutate(broad_class = case_when(
        class %in% c("Satellite_Repeat", "Satellite_Tandem_Repeat") ~ "Satellite/Tandem_Repeat",
        class %in% c("LINE_Repeat", "SINE_Repeat", "LTR_Repeat", "PLE_Repeat") ~ "Retrotransposon_Repeat",
        class %in% c("DNA_Repeat", "RC_Repeat") ~ "DNA_Transposon_Repeat",
        class %in% c("Simple_Repeat", "Low_Complexity_Repeat") ~ "Simple/Low_Complexity_Repeat",
        TRUE ~ "Other/Unknown Repeat"
    ))

# Write master list
fwrite(master_repeats, file.path(output_dir, "master_repeat_list.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
message(">>> Master repeat list saved.")


# =============================================================================
# PART 2: CHROMOSOME WINDOW PLOTS (100kb Windows)
# =============================================================================

message(">>> Generating per-chromosome landscape plots...")
all_chroms <- unique(master_repeats$chromosome)

for (chrom in all_chroms) {
    chrom_data <- master_repeats %>%
        filter(chromosome == chrom) %>%
        mutate(window = floor(start / 100000) * 100000) %>%
        group_by(window, broad_class) %>%
        summarise(coverage = sum(end - start + 1) / 100000 * 100, .groups = "drop")
    p <- ggplot(chrom_data, aes(x = window / 1e6, y = coverage, fill = broad_class)) +
        geom_area(alpha = 0.8, linewidth = 0.2) +
        facet_grid(broad_class ~ ., scales = "free_y") +
        scale_fill_brewer(palette = "Set1") +
        labs(
            title = paste("Repeat Landscape: Chromosome", chrom),
            subtitle = "100kb Windows",
            x = "Position (Mb)", y = "% Coverage", fill = "Category"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            strip.text.y = element_text(angle = 0, face = "bold"),
            legend.position = "none"
        )

    file_name <- paste0("03_Chr_", chrom, "_repeats.png")
    ggsave(file.path(plot_dir, file_name), plot = p, width = 10, height = 10, dpi = 150)
}

# =============================================================================
# PART 3: SUMMARY STATISTICS & TEXT REPORT EXPORT
# =============================================================================

message(">>> Calculating genome-wide summary statistics...")

bin_lut <- read.delim(file.path(input_dir, "bin_lut.txt"))

genome_size <- bin_lut %>%
    group_by(chromosome) %>%
    summarise(chrom_len = max(bin_end), .groups = "drop") %>%
    summarise(total_bp = sum(as.numeric(chrom_len))) %>%
    pull(total_bp)

repeats_clean <- master_repeats %>%
    mutate(
        length = end - start + 1,
        broad_class = case_when(
            class %in% c("Satellite_Repeat", "Satellite_Tandem_Repeat") ~ "Broad_Satellite/Tandem_Repeat",
            class %in% c("LINE_Repeat", "SINE_Repeat", "LTR_Repeat", "PLE_Repeat") ~ "Broad_Retrotransposon_Repeat",
            class %in% c("DNA_Repeat", "RC_Repeat") ~ "Broad_DNA_Transposon_Repeat",
            class %in% c("Simple_Repeat", "Low_Complexity_Repeat") ~ "Broad_Simple/Low_Complexity_Repeat",
            TRUE ~ "Broad_Other/Unknown_Repeat"
        )
    )

broad_summary <- repeats_clean %>%
    group_by(broad_class) %>%
    summarise(count = n(), bp = sum(as.numeric(length)), percent = (bp / genome_size) * 100, .groups = "drop") %>%
    arrange(desc(bp))

class_summary <- repeats_clean %>%
    group_by(broad_class, class) %>%
    summarise(count = n(), bp = sum(as.numeric(length)), percent = (bp / genome_size) * 100, .groups = "drop") %>%
    arrange(broad_class, desc(bp))

family_summary <- repeats_clean %>%
    group_by(broad_class, class, family) %>%
    summarise(count = n(), bp = sum(as.numeric(length)), percent = (bp / genome_size) * 100, .groups = "drop") %>%
    arrange(broad_class, class, desc(bp))

format_df_for_file <- function(df) {
    df_std <- as.data.frame(df)
    capture.output(write.table(df_std, row.names = FALSE, quote = FALSE, sep = "\t"))
}

summary_output_file <- file.path(output_dir, "master_repeat_summary.txt")
writeLines(c(
    "==================================================",
    "GENOME COMPOSITION SUMMARY",
    paste("Total Genome Size:", format(genome_size, big.mark = ","), "bp"),
    "==================================================",
    "",
    "1. BROAD CLASS SUMMARY",
    "--------------------------------------------------",
    format_df_for_file(broad_summary),
    "",
    "2. DETAILED CLASS SUMMARY",
    "--------------------------------------------------",
    format_df_for_file(class_summary),
    "",
    "3. FAMILY LEVEL SUMMARY",
    "--------------------------------------------------",
    format_df_for_file(family_summary)
), con = summary_output_file)

# =============================================================================
# PART 4: MANHATTAN PLOTS ACROSS GENOME WIDELY
# =============================================================================

message(">>> Preparing data for Manhattan plots...")

# Compute binned features dynamically from bin_lut and master_repeats if not pre-computed
bin_size <- 1e5
repeats_gr <- GRanges(
    seqnames = master_repeats$chromosome,
    ranges = IRanges(start = master_repeats$start, end = master_repeats$end),
    broad_class = master_repeats$broad_class
)

bins_gr <- makeGRangesFromDataFrame(
    bin_lut,
    keep.extra.columns = TRUE,
    seqnames.field = "chromosome",
    start.field = "bin_start",
    end.field = "bin_end"
)

# Calculate overlap proportion per bin
overlaps <- findOverlaps(bins_gr, repeats_gr)
hit_df <- data.frame(
    bin_idx = queryHits(overlaps),
    repeat_class = mcols(repeats_gr)$broad_class[subjectHits(overlaps)],
    overlap_len = width(pintersect(bins_gr[queryHits(overlaps)], repeats_gr[subjectHits(overlaps)]))
)

bin_summary <- hit_df %>%
    group_by(bin_idx, repeat_class) %>%
    summarise(overlap_bp = sum(overlap_len), .groups = "drop") %>%
    mutate(proportion = overlap_bp / bin_size)

# Merge back with bin_lut
binned_features <- bin_lut %>%
    mutate(bin_idx = row_number()) %>%
    left_join(bin_summary, by = "bin_idx") %>%
    mutate(proportion = replace_na(proportion, 0)) %>%
    pivot_wider(names_from = repeat_class, values_from = proportion, values_fill = 0)

# Ensure all Broad class columns exist
expected_cols <- c(
    "Broad_Satellite/Tandem_Repeat", "Broad_Retrotransposon_Repeat",
    "Broad_DNA_Transposon_Repeat", "Broad_Simple/Low_Complexity_Repeat", "Broad_Other/Unknown_Repeat"
)
for (col in expected_cols) {
    if (!col %in% names(binned_features)) binned_features[[col]] <- 0
}

repeats_long <- binned_features %>%
    pivot_longer(
        cols = starts_with("Broad_"),
        names_to = "repeat_class",
        values_to = "proportion"
    ) %>%
    mutate(repeat_class = gsub("Broad_|_Repeat", "", repeat_class))

# Chromosome ordering and offsets
all_chrs <- unique(as.character(bin_lut$chromosome))
numeric_chrs <- all_chrs[grepl("^[0-9]+$", all_chrs)]
numeric_chrs <- numeric_chrs[order(as.numeric(numeric_chrs))]
final_chr_order <- c(numeric_chrs, "W", "Z")

chr_offsets <- bin_lut %>%
    mutate(chromosome = factor(chromosome, levels = final_chr_order)) %>%
    group_by(chromosome) %>%
    summarise(chrom_len = max(bin_end), .groups = "drop") %>%
    arrange(chromosome) %>%
    mutate(tot = cumsum(as.numeric(lag(chrom_len, default = 0))))

repeats_long <- repeats_long %>%
    mutate(chromosome = factor(chromosome, levels = final_chr_order)) %>%
    left_join(chr_offsets, by = "chromosome") %>%
    mutate(pos_cum = (bin * bin_size) + tot)
    
axis_df <- repeats_long %>%
    group_by(chromosome) %>%
    summarize(center = mean(pos_cum), .groups = "drop") %>%
    arrange(chromosome)

chr_colors <- rep(c("grey40", "skyblue3"), length.out = length(final_chr_order))
names(chr_colors) <- final_chr_order
if ("W" %in% names(chr_colors)) chr_colors["W"] <- "grey40"
if ("Z" %in% names(chr_colors)) chr_colors["Z"] <- "skyblue3"

repeats_long <- repeats_long %>%
    mutate(repeat_class = case_when(
        repeat_class == "Retrotransposon" ~ "Retrotransposon",
        repeat_class == "DNA_Transposon" ~ "DNA Transposon",
        repeat_class == "Other.Unknown" ~ "Other / Unknown Repeat",
        repeat_class == "Satellite.Tandem" ~ "Satellite / Tandem Repeat",
        repeat_class == "Simple.Low_Complexity" ~ "Simple / Low Complexity Repeat",
        TRUE ~ repeat_class
    )) %>%
    mutate(repeat_class = factor(repeat_class, levels = c(
        "DNA Transposon",
        "Retrotransposon",
        "Other / Unknown Repeat",
        "Satellite / Tandem Repeat",
        "Simple / Low Complexity Repeat"
    )))

# --- Full Genome Manhattan Plot ---
message(">>> Generating full-genome Manhattan plot...")
ggplot(repeats_long, aes(x = pos_cum / 1e6, y = proportion)) +
    geom_point(aes(color = chromosome), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = chr_colors) +
    scale_x_continuous(labels = axis_df$chromosome, breaks = axis_df$center / 1e6) +
    facet_wrap(~repeat_class, ncol = 1, scales = "free_y") +
    labs(
        title = "Repeat Distribution Across the Sparrow Genome"
        x = "Chromosome", y = "Proportion"
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank()
    )

ggsave(file.path(plot_dir, "03_Repeat_Distribution_Manhattan_Plot.png"), width = 48, height = 8, dpi = 300)

message(">>> Pipeline execution complete successfully!")

rm(list = ls())
gc()

# End of script