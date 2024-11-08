# Load necessary libraries
library(rjson)
library(tidyverse)

# Set working directory
setwd("C:/Users/MSI-NB/Jupyter Code/TCGA_R")

# Load JSON metadata
json_data <- jsonlite::fromJSON("metadata.cart.2024-11-08.json")
View(json_data)

sample_ids <- sapply(json_data$associated_entities, function(x) x[, 1])
sample_ids[1:10]

file_sample <- data.frame(sample_id = sample_ids, file_name = json_data$file_name)
View(file_sample)

# List count files
count_files <- list.files('gdc_download_20241108_123324.776410/', pattern = '*.tsv', recursive = TRUE)
count_files[1:10]

count_file_names <- sapply(strsplit(count_files, split = '/'), function(x) x[2])
count_file_names[1:10]

# Initialize main matrix with 60660 rows (based on assumption)
matrix <- data.frame(matrix(nrow = 60660, ncol = 0))

# Function to read and process each file
process_file <- function(file_path, sample_id) {
  data <- read.delim(file_path, fill = TRUE, header = FALSE, row.names = 1, skip = 1)

  # Ensure it has at least 6 columns
  if (ncol(data) < 6) stop("File format error: Less than 6 columns.")

  # Extract 6th column, clean it, convert to numeric, handle NAs
  tpm_values <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", data[[6]])))
  tpm_values[is.na(tpm_values)] <- 0

  # Ensure exactly 60660 rows
  length_diff <- 60660 - length(tpm_values)
  if (length_diff > 0) {
    tpm_values <- c(tpm_values, rep(0, length_diff))
  } else if (length_diff < 0) {
    tpm_values <- tpm_values[1:60660]
  }

  # Return as named data frame column
  return(data.frame(sample_id = sample_id, tpm_values = tpm_values))
}

# Process each file and bind to main matrix
for (i in seq_along(count_files)) {
  path <- file.path("gdc_download_20241108_123324.776410", count_files[i])
  sample_id <- file_sample$sample_id[file_sample$file_name == count_file_names[i]]

  matrix <- cbind(matrix, setNames(process_file(path, sample_id), sample_id))
}

# Extract gene names and types from the first file
initial_data <- read.delim(file.path("gdc_download_20241108_123324.776410", count_files[1]), fill = TRUE, header = FALSE, row.names = 1)
gene_names <- initial_data[-c(1:6), 1]
gene_names[1:10]

gene_types <- initial_data[-c(1:6), 2]
gene_types[1:10]

# Confirm gene_names and gene_types lengths match
stopifnot(length(gene_names) == length(gene_types))

# Combine gene data with main matrix
matrix0 <- data.frame(gene_type = gene_types, gene_name = gene_names, matrix)

# Check matrix0 before aggregation
cat("Before aggregation and filtering:\n")
print(dim(matrix0))

# Perform aggregation by gene_name, taking the maximum across columns
if (nrow(matrix0) > 0) {
  matrix0 <- aggregate(. ~ gene_name, data = matrix0, FUN = max)
  cat("After aggregation:\n")
  print(dim(matrix0))

  # Filter for protein-coding genes
  matrix0 <- subset(matrix0, gene_type == "protein_coding")
  cat("After filtering for protein-coding genes:\n")
  print(dim(matrix0))

  # Check if matrix0 still has rows after filtering
  if (nrow(matrix0) > 0) {
    # Set row names by gene_name and drop unnecessary columns
    rownames(matrix0) <- matrix0$gene_name
    matrix0 <- as.data.frame(matrix0[, -c(1, 2), drop = FALSE])

    cat("After removing gene_name and gene_type columns:\n")
    print(dim(matrix0))

    # Construct matrix1 with ID column
    matrix1 <- data.frame(ID = rownames(matrix0), matrix0)
    cat("Final matrix1 structure:\n")
    print(dim(matrix1))

    # Replace dots in column names with hyphens and export the final matrix
    colnames(matrix1) <- gsub("[.]", "-", colnames(matrix1))
    write.table(matrix1, 'TCGA_LUSC_TPM.txt', sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    stop("matrix0 is empty after filtering for protein-coding genes.")
  }
} else {
  stop("matrix0 is empty before aggregation.")
}
