#### Assign estimated crossover positions to the crossover dataset and create binned population-level crossover datasets
#### Updated 13/10/2026 by Lisa Ammer

rm(list = ls())
gc()

# Load necessary libraries
library(dplyr)
library(purrr)

# Set global seed for reproducibility
set.seed(123)

# Set Directories (Unified input/output structure)
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/1_in_/binned_features"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# PART 1: ASSIGN ESTIMATED CROSSOVER POSITIONS & ACC SUMMARY
# =============================================================================

message(">>> Step 1: Processing raw crossover positions and generating acc_summary...")

# Load crossover data
co_file <- file.path(input_dir, "_2_rectab_Crossover_Positions_post_QC.txt")
co_data <- read.table(co_file, header = TRUE, sep = "\t") %>%
    select(parent, sex, meiosis, chrom, left, right, CO_window, Total_COs_on_Chr) %>%
    dplyr::rename(chromosome = chrom, width = CO_window) %>%
    mutate(chromosome = as.character(chromosome))

# Get Total_COs_on_Chr per meiosis_id x chromosome combination
acc_summary <- co_data %>%
    group_by(meiosis, chromosome, sex) %>%
    summarise(
        Total_COs_on_Chr = dplyr::n(),
        .groups = "drop"
    ) %>%
    select(
        parent_sex = sex,
        meiosis,
        chromosome,
        acc = Total_COs_on_Chr
    )

# Save acc_summary
write.table(acc_summary, file.path(input_dir, "changed/acc_summary_per_meio.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Function to assign a random position between 'left' and 'right'
generate_est_CO_pos <- function(data) {
    data$position <- mapply(function(current_window_start, current_window_end) {
        if (is.na(current_window_start) || is.na(current_window_end) || current_window_start > current_window_end) {
            return(NA)
        }
        sample(seq(current_window_start, current_window_end), 1)
    }, data$left, data$right)

    return(data)
}

# Generate estimated crossover positions
co_data <- generate_est_CO_pos(co_data)
co_data$position <- as.numeric(co_data$position)

message(">>> Step 1 complete.")

# =============================================================================
# PART 2: CREATE POPULATION CROSSOVER DATASETS
# =============================================================================

message(">>> Step 2: Building population crossover datasets...")

# Load bin_lut
bin_lut <- read.table(file.path(input_dir, "changed/bin_lut.txt"), header = TRUE, sep = "\t") %>%
    mutate(chromosome = as.character(chromosome))

# Helper functions for dataset info and summarizing
get_dataset_info <- function(data) {
    distinct_meioses <- data %>% distinct(meiosis, sex)
    distinct_parents <- data %>% distinct(parent, sex)

    res <- list(
        total_co               = nrow(data),
        co_ratio               = sum(data$sex == "M") / sum(data$sex == "F"),
        total_meioses          = nrow(distinct_meioses),
        meioses_ratio          = sum(distinct_meioses$sex == "M") / sum(distinct_meioses$sex == "F"),
        total_parents          = nrow(distinct_parents),
        parents_ratio          = sum(distinct_parents$sex == "M") / sum(distinct_parents$sex == "F"),
        avg_width_male         = mean(data$width[data$sex == "M"], na.rm = TRUE),
        avg_width_female       = mean(data$width[data$sex == "F"], na.rm = TRUE),
        range_width_male_min   = range(data$width[data$sex == "M"], na.rm = TRUE)[1],
        range_width_male_max   = range(data$width[data$sex == "M"], na.rm = TRUE)[2],
        range_width_female_min = range(data$width[data$sex == "F"], na.rm = TRUE)[1],
        range_width_female_max = range(data$width[data$sex == "F"], na.rm = TRUE)[2]
    )
    return(as.data.frame(res))
}

summarise_co_by_bin <- function(co_subset, bin_subset) {
    counts <- co_subset %>%
        group_by(chromosome, bin) %>%
        summarise(
            male_co = sum(sex == "M", na.rm = TRUE),
            female_co = sum(sex == "F", na.rm = TRUE),
            .groups = "drop"
        )

    co_summary <- bin_subset %>%
        left_join(counts, by = c("chromosome", "bin")) %>%
        mutate(
            male_co = coalesce(male_co, 0L),
            female_co = coalesce(female_co, 0L)
        )

    return(co_summary)
}

# Define filtering logic for Datasets A - E
co_list <- list()
bin_list <- list()

# Dataset A: Base dataset
co_list$A <- co_data %>%
    left_join(bin_lut, by = join_by(chromosome == chromosome, position >= bin_start, position <= bin_end))
bin_list$A <- bin_lut

# Dataset B: Filter for widths <= 95th percentile in Dataset A
co_list$B <- co_list$A %>% filter(width <= quantile(width, 0.95, na.rm = TRUE))
bin_list$B <- bin_list$A

# Update info table and Export
dataset_names <- c("A", "B")
info_list <- map(dataset_names, ~ get_dataset_info(co_list[[.x]]))

co_dataset_info <- bind_rows(info_list) %>%
    mutate(dataset = dataset_names, .before = 1)

write.table(co_dataset_info, file.path("/mnt/loki/johnston/lisa/sparrow/2026-10_sparrow_recom/2_out_/co_dataset_info.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# Run bin summaries and export files
walk(dataset_names, function(name) {
    summary_data <- summarise_co_by_bin(co_list[[name]], bin_list[[name]])
    write.table(summary_data,
        file.path(output_dir, paste0("co_dataset", name, ".txt")),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
})

message(">>> Integrated crossover pipeline complete! All datasets and summaries saved.")

rm(list = ls())
gc()

# End of script