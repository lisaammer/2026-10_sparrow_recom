####  Investigate possible centromere positions for chromosomes
####  Updated 16/10/2026 by Lisa Ammer

# Clear workspace and run garbage collection
rm(list = ls())
gc()

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# Directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/centromere_candidates"

####
#### Part 1: Create equal-length panels for chromosomes to visualize potential centromere positions
####

# Read filtered sequence report
report <- read.table(file.path(input_dir, "changed/filtered_sequence_report.txt"), header = TRUE, sep = "\t") %>%
    mutate(chrName = as.character(chrName))

# Create a data frame with chromosome lengths
chrom_data <- report %>%
    select(chr = chrName, length) %>%
    mutate(chr = as.character(chr)) %>%
    rename(chromosome = chr)

# Bin-packing function (First-Fit Decreasing)
pack_chromosomes <- function(df) {
    df <- df %>% arrange(desc(length))
    max_len <- df$length[1]

    bins <- list()
    bin_sums <- c()

    for (i in 1:nrow(df)) {
        chrom <- df$chromosome[i]
        len <- df$length[i]
        placed <- FALSE

        if (length(bin_sums) > 0) {
            for (j in seq_along(bin_sums)) {
                if (bin_sums[j] + len <= max_len) {
                    bins[[j]] <- rbind(bins[[j]], df[i, ])
                    bin_sums[j] <- bin_sums[j] + len
                    placed <- TRUE
                    break
                }
            }
        }

        if (!placed) {
            bins[[length(bins) + 1]] <- df[i, , drop = FALSE]
            bin_sums <- c(bin_sums, len)
        }
    }

    return(list(max_len = max_len, bins = bins))
}

packed_result <- pack_chromosomes(chrom_data)
max_len <- packed_result$max_len

# Process each bin to create local plotting coordinates
plot_list <- lapply(seq_along(packed_result$bins), function(i) {
    bin_df <- packed_result$bins[[i]]

    # Order chromosome factor levels by length (descending)
    bin_df <- bin_df %>%
        mutate(
            chromosome = factor(chromosome, levels = chromosome[order(length, decreasing = TRUE)]),
            end_cum = cumsum(length),
            start_cum = end_cum - length
        )

    # Dynamically generate enough colors for any number of chromosomes (handles >12 safely)
    n_chroms <- nrow(bin_df)
    custom_colors <- if (n_chroms <= 12) {
        RColorBrewer::brewer.pal(max(n_chroms, 3), "Set3")[1:n_chroms]
    } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_chroms)
    }

    p <- ggplot(bin_df) +
        geom_rect(aes(xmin = start_cum, xmax = end_cum, ymin = 0, ymax = 1, fill = chromosome), alpha = 0.2) +
        # Add vertical lines for chromosome boundaries (skip the first boundary start)
        {
            if (nrow(bin_df) > 1) geom_vline(xintercept = bin_df$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        xlim(0, max_len) +
        labs(
            title = paste("Panel", i),
            x = "Genome Coordinates (bp)",
            y = "Value"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")

    return(p)
})

# Combine all panels vertically into a single view
final_plot <- wrap_plots(plot_list, ncol = 1)
ggsave(filename = file.path(output_dir, "chromosome_panels.png"), plot = final_plot, width = 30, height = 15, dpi = 300)

###
### Part 2: Identify potential centromere positions based on repeat density, recombination rates, GC patterns, TSS density, and SNP density
###

# Import binned datasets
binned_tss <- read.table(file.path(input_dir, "binned_features/binned_tss_and_cgi.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome)) %>%
    select(-c(tss_with_cgi, tss_no_cgi))
binned_gc <- read.table(file.path(input_dir, "binned_features/binned_gc_prop.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))
    binned_repeats <- read.table(file.path(input_dir, "binned_features/binned_repeat_prop.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))
binned_recomb <- read.table(file.path(input_dir, "binned_features/co_datasetA.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))
binned_snp <- read.table(file.path(input_dir, "binned_features/binned_snp_distribution.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))

# Print colnames of each binned dataset for verification
message(">>> Column names of binned datasets:")
message("Binned GC: ", paste(colnames(binned_gc), collapse = ", "))
message("Binned TSS: ", paste(colnames(binned_tss), collapse = ", "))
message("Binned Repeats: ", paste(colnames(binned_repeats), collapse = ", "))
message("Binned Recombination: ", paste(colnames(binned_recomb), collapse = ", "))
message("Binned SNP: ", paste(colnames(binned_snp), collapse = ", "))

# 1. Merge all binned datasets together by chromosome and bin
# (Using binned_recomb for bin_start/bin_end coordinates)
binned_all <- binned_gc %>%
    left_join(binned_tss, by = c("chromosome", "bin")) %>%
    left_join(binned_repeats, by = c("chromosome", "bin")) %>%
    left_join(binned_recomb %>% select(chromosome, bin, bin_start, bin_end, male_co, female_co), by = c("chromosome", "bin")) %>%
    left_join(binned_snp, by = c("chromosome", "bin"))

# 2. Loop through each bin panel to generate the 5-track plots
panel_feature_plots <- lapply(seq_along(packed_result$bins), function(i) {
    bin_df <- packed_result$bins[[i]]

    # Calculate cumulative offsets for chromosomes in this panel
    bin_df_offset <- bin_df %>%
        mutate(
            chromosome = factor(chromosome, levels = chromosome[order(length, decreasing = TRUE)]),
            end_cum = cumsum(length),
            start_cum = end_cum - length
        )

    # Filter binned data for chromosomes in this specific panel and add global coordinates
    panel_data <- binned_all %>%
        inner_join(bin_df_offset, by = "chromosome") %>%
        mutate(
            global_start = start_cum + bin_start,
            global_end = start_cum + bin_end,
            global_mid = (global_start + global_end) / 2
        )

    # Common theme adjustments for alignment
    track_theme <- theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = margin(2, 5, 2, 5)
        )

    # Track 1: Total TSS (Line)
    p_tss <- ggplot(panel_data, aes(x = global_mid, y = total_tss, color = chromosome)) +
        geom_line(alpha = 0.8) +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        xlim(0, max_len) +
        labs(y = "TSS", title = paste("Panel", i)) +
        track_theme +
        theme(legend.position = "none")

    # Track 2: GC Proportion (Line)
    p_gc <- ggplot(panel_data, aes(x = global_mid, y = gc_prop, color = chromosome)) +
        geom_line(alpha = 0.8)  +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        xlim(0, max_len) +
        labs(y = "GC Prop") +
        track_theme +
        theme(legend.position = "none")

    # Track 3: Different Repeat Classes (Multiple colored lines)
    repeat_data_long <- panel_data %>%
        pivot_longer(
            cols = c(
                Retrotransposon_Repeat, Satellite.Tandem_Repeat,
                Simple.Low_Complexity_Repeat, DNA_Transposon_Repeat, Other.Unknown.Repeat
            ),
            names_to = "repeat_class",
            values_to = "repeat_prop"
        )

    p_repeats <- ggplot(repeat_data_long, aes(x = global_mid, y = repeat_prop, color = repeat_class)) +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        geom_line(alpha = 0.8) +
        xlim(0, max_len) +
        labs(y = "Repeats", color = "Repeat Class") +
        track_theme +
        theme(legend.position = "right") # Or "none" if legend gets too cluttered

    # Track 4: Recombination (Male = Blue, Female = Red)
    p_recomb <- ggplot(panel_data, aes(x = global_mid)) +
        geom_line(aes(y = male_co, color = "Male CO"), alpha = 0.8) +
        geom_line(aes(y = female_co, color = "Female CO"), alpha = 0.8)  +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        scale_color_manual(values = c("Male CO" = "blue", "Female CO" = "red")) +
        xlim(0, max_len) +
        labs(y = "Crossovers", color = "Sex") +
        track_theme +
        theme(legend.position = "right")

    # Track 5: SNP Distribution (Line) + X-axis enabled here
    p_snp <- ggplot(panel_data, aes(x = global_mid, y = snp_distribution, color = chromosome)) +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        geom_line(alpha = 0.8) +
        xlim(0, max_len) +
        labs(x = "Genome Coordinates (bp)", y = "SNPs") +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            legend.position = "none",
            plot.margin = margin(2, 5, 2, 5)
        )

    # Track 6: Chromosome indicator bar at the bottom to serve as the chromosome legend reference
    p_chroms <- ggplot(bin_df_offset) +
        geom_rect(aes(xmin = start_cum, xmax = end_cum, ymin = 0, ymax = 1, fill = chromosome), alpha = 0.7) +
        # Add vertical dashed lines for chromosome boundaries
        {
            if (nrow(bin_df_offset) > 1) geom_vline(xintercept = bin_df_offset$start_cum[-1], linetype = "dashed", color = "grey")
        } +
        xlim(0, max_len) +
        labs(fill = "Chromosome") +
        theme_void() +
        theme(legend.position = "bottom", legend.box = "horizontal")

    # Combine vertically using patchwork
    combined_panel <- (p_tss / p_gc / p_repeats / p_recomb / p_snp / p_chroms) +
        plot_layout(heights = c(1, 1, 1, 1, 1, 0.3))

    return(combined_panel)
})

# 3. Save individual plots or wrap them together
for (i in seq_along(panel_feature_plots)) {
    ggsave(
        filename = file.path(output_dir, paste0("feature_distribution_panel_", i, ".png")),
        plot = panel_feature_plots[[i]],
        width = 45, height = 15, dpi = 300
    )
}

message(">>> All panel feature distribution plots successfully generated and saved!")

# Save processed data into a single file for local use
saveRDS(list(packed_result = packed_result, binned_all = binned_all),
    file = file.path(output_dir, "shiny_app_data.rds")
)

rm(list = ls())
gc()

# End of script