#### Calculate the binned germline-specific tpm for protein-coding genes
#### Updated 12/10/2026 by Lisa Ammer

rm(list = setdiff(ls(), c("germ", "bin_lut", "genes_with_bins")))
gc()

# Load necessary libraries
library(dplyr)
library(Seurat)
library(GenomicRanges)
library(tidyr)

# Set global seed for reproducibility
set.seed(42)

# Set Directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/changed"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Seurat germ cell object path
seurat_path <- "/mnt/loki/johnston/lisa/sparrow/snRNAseq/integrated_downstream/5_germ_cell_subset_ovary.rds"

# If not loaded, load the germ cell data
if (!exists("germ")) {
    germ <- readRDS(seurat_path)
}

# If not loaded, load the bin lookup table
if (!exists("bin_lut")) {
    bin_lut <- read.table(
        file.path(input_dir, "bin_lut.txt"),
        header = TRUE,
        sep = "\t"
    ) %>% mutate(chromosome = as.character(chromosome))
}

# If not loaded, load genes with bins and CGI status
if (!exists("all_genes_with_bins")) {
    genes_with_bins <- read.table(
        file.path(input_dir, "all_genes_with_bins_and_cgi_status.txt"),
        header = TRUE,
        sep = "\t"
    ) %>% mutate(chromosome = as.character(chromosome))
}

# Extract raw counts
raw_counts <- germ[["RNA"]]$counts

# Find genes present in both Seurat object and your annotation table
common_genes <- intersect(rownames(raw_counts), genes_with_bins$gene)

# Subset counts and metadata to ensure 1-to-1 mapping
counts_subset <- raw_counts[common_genes, ]
gene_metadata <- genes_with_bins %>%
    filter(gene %in% common_genes) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(match(gene, common_genes))

# Filter for genes with non-zero exonic length to avoid division by zero
genes_with_lengths <- !is.na(gene_metadata$exonic_length)
counts_subset <- counts_subset[genes_with_lengths, ]
gene_metadata <- gene_metadata[genes_with_lengths, ]

# RPK = Counts / (exonic_length / 1000)
rpk_matrix <- as.matrix(counts_subset) / (gene_metadata$exonic_length / 1000)

# TPM = (RPK / Column Sum of RPKs) * 1e6
col_sums_rpk <- colSums(rpk_matrix)
tpm_matrix <- sweep(rpk_matrix, 2, col_sums_rpk, "/") * 1e6
tpm_matrix[is.na(tpm_matrix)] <- 0 # Handle cells with 0 total counts

# Get individual IDs from Seurat metadata
cell_info <- germ@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, sample_id)

unique_samples <- unique(cell_info$sample_id)

# Calculate Mean TPM per gene per Individual
sample_means_list <- lapply(unique_samples, function(s_id) {
    cells <- cell_info$cell_id[cell_info$sample_id == s_id]
    rowMeans(tpm_matrix[, cells, drop = FALSE])
})

names(sample_means_list) <- unique_samples
gene_tpm_by_ind <- as.data.frame(sample_means_list) %>%
    mutate(gene = rownames(tpm_matrix))

# Combine Individual means with Bin info
binned_tpm <- gene_tpm_by_ind %>%
    tidyr::pivot_longer(
        cols = -gene,
        names_to = "individual_id",
        values_to = "mean_tpm"
    ) %>%
    left_join(gene_metadata %>% dplyr::select(gene, chromosome, bin), by = "gene") %>%
    group_by(individual_id, chromosome, bin) %>%
    summarise(TPM_bin = sum(mean_tpm), .groups = "drop") %>%
    mutate(chromosome = as.character(chromosome))

# Ensure even bins with ZERO expression are represented
binned_tpm_per_id <- bin_lut %>%
    dplyr::select(chromosome, bin) %>%
    tidyr::expand_grid(individual_id = unique_samples) %>%
    left_join(binned_tpm, by = c("chromosome", "bin", "individual_id")) %>%
    mutate(TPM_bin = tidyr::replace_na(TPM_bin, 0))

# Save output
write.table(
    binned_tpm_per_id,
    file.path(output_dir, "binned_tpm_per_id.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# Get average TPM per bin across all individuals
binned_tpm_avg <- binned_tpm_per_id %>%
    group_by(chromosome, bin) %>%
    summarise(avg_tpm = mean(TPM_bin), .groups = "drop")

# Save output
write.table(
    binned_tpm_avg,
    file.path(output_dir, "binned_tpm_avg.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# Calculate mean TPM for each gene across all individuals
gene_tpm_by_ind$global_gene_mean <- rowMeans(gene_tpm_by_ind %>% dplyr::select(-gene))

# Activity threshold: global mean >= 1
gene_tpm_by_ind <- gene_tpm_by_ind %>%
    mutate(is_active = if_else(global_gene_mean >= 1, 1, 0))

# Calculate active gene counts per bin
binned_activity <- gene_tpm_by_ind %>%
    dplyr::select(gene, is_active) %>%
    left_join(gene_metadata %>% dplyr::select(gene, chromosome, bin), by = "gene") %>%
    group_by(chromosome, bin) %>%
    summarise(active_gene_count = sum(is_active), .groups = "drop") %>%
    mutate(chromosome = as.character(chromosome))

# Load the total TSS count to get the denominator
total_tss_per_bin <- gene_metadata %>%
    group_by(chromosome, bin) %>%
    summarise(total_tss = n(), .groups = "drop") %>%
    mutate(chromosome = as.character(chromosome))

# Combine and calculate proportion
binned_activity_final <- bin_lut %>%
    mutate(chromosome = as.character(chromosome)) %>%
    left_join(total_tss_per_bin, by = c("chromosome", "bin")) %>%
    left_join(binned_activity, by = c("chromosome", "bin")) %>%
    mutate(
        total_tss = tidyr::replace_na(total_tss, 0),
        active_gene_count = tidyr::replace_na(active_gene_count, 0),
        prop_active = if_else(total_tss > 0, active_gene_count / total_tss, 0)
    ) %>%
    dplyr::select(chromosome, bin, active_gene_count, active_gene_prop = prop_active)

# Save the new activity metrics
write.table(
    binned_activity_final,
    file.path(output_dir, "binned_prop_expr.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# Determine active genes per individual using the fixed TPM >= 1 threshold
gene_activity_per_id <- gene_tpm_by_ind %>%
    dplyr::select(gene, all_of(unique_samples)) %>%
    tidyr::pivot_longer(
        cols = -gene,
        names_to = "individual_id",
        values_to = "tpm"
    ) %>%
    mutate(
        is_active = if_else(tpm >= 1, 1, 0)
    )

# Join with bin metadata
binned_active_counts_per_id <- gene_activity_per_id %>%
    left_join(gene_metadata %>% dplyr::select(gene, chromosome, bin), by = "gene") %>%
    group_by(individual_id, chromosome, bin) %>%
    summarise(active_gene_count = sum(is_active), .groups = "drop") %>%
    mutate(chromosome = as.character(chromosome))

# Merge and calculate proportion per individual
binned_prop_per_id <- bin_lut %>%
    dplyr::select(chromosome, bin) %>%
    tidyr::expand_grid(individual_id = unique_samples) %>%
    left_join(binned_active_counts_per_id,
        by = c("individual_id", "chromosome", "bin")
    ) %>%
    left_join(total_tss_per_bin, by = c("chromosome", "bin")) %>%
    mutate(
        active_gene_count = tidyr::replace_na(active_gene_count, 0),
        total_tss = tidyr::replace_na(total_tss, 0),
        prop_active = if_else(total_tss > 0, active_gene_count / total_tss, 0)
    ) %>%
    dplyr::select(chromosome, bin, individual_id, active_gene_count, active_gene_prop = prop_active)

# Save the output
write.table(
    binned_prop_per_id,
    file.path(output_dir, "binned_prop_expr_per_id.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# =============================================================================
# CGI PROMOTER ACTIVITY SPLIT (Using pre-calculated has_cgi status)
# =============================================================================

# Join activity data with the pre-calculated has_cgi and gene_biotype from gene_metadata
gene_activity_cgi <- gene_tpm_by_ind %>%
    left_join(
        gene_metadata %>%
            dplyr::select(gene, chromosome, bin, has_cgi, gene_biotype),
        by = "gene"
    )

# Calculate bin-level metrics using has_cgi
binned_cgi_split <- gene_activity_cgi %>%
    group_by(chromosome, bin) %>%
    summarise(
        # ALL active genes
        all_active_with_cgi = sum(is_active == 1 & has_cgi == TRUE, na.rm = TRUE),
        all_active_no_cgi = sum(is_active == 1 & has_cgi == FALSE, na.rm = TRUE),

        # PROTEIN-CODING active genes
        pc_active_with_cgi = sum(is_active == 1 & has_cgi == TRUE & gene_biotype == "protein_coding", na.rm = TRUE),
        pc_active_no_cgi = sum(is_active == 1 & has_cgi == FALSE & gene_biotype == "protein_coding", na.rm = TRUE),

        # NON PROTEIN-CODING active genes (everything else)
        not_pc_active_with_cgi = sum(is_active == 1 & has_cgi == TRUE & gene_biotype != "protein_coding", na.rm = TRUE),
        not_pc_active_no_cgi = sum(is_active == 1 & has_cgi == FALSE & gene_biotype != "protein_coding", na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(chromosome = as.character(chromosome))

# Merge with bin_lut to ensure zero-counts are preserved
binned_cgi_final <- bin_lut %>%
    mutate(chromosome = as.character(chromosome)) %>%
    left_join(binned_cgi_split, by = c("chromosome", "bin")) %>%
    mutate(across(starts_with(c("all_", "pc_", "not_pc_")), ~ tidyr::replace_na(., 0))) %>%
    dplyr::select(
        chromosome, bin,
        all_active_with_cgi, all_active_no_cgi,
        pc_active_with_cgi, pc_active_no_cgi,
        not_pc_active_with_cgi, not_pc_active_no_cgi
    )

# Save the variables for your INLA model
write.table(
    binned_cgi_final,
    file.path(output_dir, "binned_active_genes_by_promoter_type.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

message(">>> Germline TPM and CGI activity pipeline complete!")

rm(list = ls())
gc()

# End of script