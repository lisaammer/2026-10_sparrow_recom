#### Compile sequence coverage per chromsome for all samples
#### Created 06/02/2026 by Lisa Ammer
#### Ran this on ac3

rm(list = ls()) # Clear the environment to ensure a clean run
gc()

# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)
library(fs)

# Define input and output directories
input_dir <- "/mnt/loki/johnston/lisa/sparrow/2026-01_sparrow_recom/yapp_coverage"
output_dir <- "/mnt/loki/johnston/lisa/sparrow/reference/processed_yapp_coverage"

# Create the output directory if it doesn't exist
dir_create(output_dir)

# Get a list of all files in the input directory
file_list <- list.files(input_dir, full.names = TRUE)

# Check if any files were found in the input directory
if (length(file_list) == 0) {
    message("No files found in the input directory: ", input_dir)
} else {
    message(paste("Found", length(file_list), "files to process."))

    # Loop through each file and apply transformations
    for (file_path in file_list) {
        # Extract the filename from the full path
        file_name <- basename(file_path)
        message(paste("Processing file:", file_name))

        tryCatch(
            {
                # Read the whitespace-delimited file using readr::read_table
                df <- readr::read_table(
                    file_path,
                    col_names = c("parent", "sex", "offspring", "chrom", "left", "right"),
                    # Explicitly set column names
                    skip = 1 # Skip the original header row
                )

                df_transformed <- df %>%
                    dplyr::rename(
                        chromosome = chrom,
                        cov_start = left,
                        cov_end = right
                    ) %>%
                    mutate(
                        chromosome = as.character(chromosome),
                        # Ensure it's character for str_pad
                        chromosome = str_pad(
                            chromosome,
                            width = 2,
                            side = "left",
                            pad = "0"
                        ),
                        chromosome = as.factor(chromosome)
                    )

                # Define the output path for the transformed file
                output_file_path <- file.path(output_dir, file_name)

                # Save the transformed data frame to the new folder
                write.table(
                    df_transformed,
                    output_file_path,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE
                )

                message(paste("Successfully processed and saved:", output_file_path))
            },
            error = function(e) {
                message(paste("Error processing file", file_name, ":", e$message))
            }
        )
    }
    message("All files processed.")
}

# Compile information from processed files

message("\nCompiling information from processed files...")

# Get a list of all processed files in the output directory
processed_file_list <- list.files(output_dir, full.names = TRUE)

# Initialize an empty list to store data frames
all_data_list <- list()

# Check if any processed files were found
if (length(processed_file_list) == 0) {
    message("No processed files found in the output directory: ", output_dir)
    compiled_table <- data.frame() # Create an empty data frame if no files
} else {
    # Loop through each processed file and read it
    for (file_path in processed_file_list) {
        file_name <- basename(file_path)
        message(paste("Reading processed file:", file_name))

        tryCatch(
            {
                # Read the processed file (which is tab-delimited)
                # We expect 'chromosome' to be a factor and other columns as they were saved.
                df_processed <- read.delim(
                    file_path,
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE
                ) %>%
                    # Ensure chromosome is treated as factor as it was saved
                    mutate(chromosome = as.factor(chromosome))

                # Select only the required columns
                df_selected <- df_processed %>%
                    select(parent, offspring, chromosome, cov_start, cov_end)

                # Add the selected data frame to the list
                all_data_list[[file_name]] <- df_selected
            },
            error = function(e) {
                message(paste("Error reading processed file", file_name, ":", e$message))
            }
        )
    }

    # Combine all data frames in the list into a single data frame
    if (length(all_data_list) > 0) {
        compiled_table <- bind_rows(all_data_list)
        message(paste(
            "Successfully compiled data from",
            length(all_data_list),
            "files."
        ))
    } else {
        message("No data frames were successfully read from processed files.")
        compiled_table <- data.frame() # Create an empty data frame if no data was read
    }
}

# Save the compiled table to file
write.table(
    compiled_table,
    file.path(output_dir, "compiled_coverage_data.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

message(paste(
    "Compiled table saved to:",
    file.path(output_dir, "compiled_coverage_data.txt")
))


# Calculate minimum cov_start and maximum cov_end per chromosome
max_coverage_per_chromosome <- compiled_table %>%
    group_by(chromosome) %>%
    summarise(
        min_cov_start = min(cov_start, na.rm = TRUE),
        max_cov_end = max(cov_end, na.rm = TRUE)
    ) %>%
    ungroup() # Ungroup the data for further operations if needed

# Write to output directory
write.table(
    max_coverage_per_chromosome,
    file.path(output_dir, "max_coverage_per_chromosome.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
message(paste(
    "Max coverage per chromosome saved to:",
    file.path(output_dir, "max_coverage_per_chromosome.txt")
))

# End of script