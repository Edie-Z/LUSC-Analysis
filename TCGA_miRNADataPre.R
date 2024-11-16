# Load necessary libraries
library(stringr)  # For string manipulation functions
library(tidyverse)  # Includes dplyr, ggplot2, and other packages for data manipulation and visualization
library(miRBaseVersions.db)  # For miRNA-related databases

# Set working directory
setwd("D:/R/Rcode/LUSC_Analysis")

# List all files that match the updated pattern (isoforms.quantification.txt)
count_files <- dir("TCGA_Dataset/miRNA/gdc_download_20241115_102256.015891/",
                   pattern = "*isoforms.quantification.txt$",  # Match files with this pattern
                   recursive = TRUE)  # Recursively search subdirectories

# Check if any files were found
if (length(count_files) == 0) {
  available_files <- list.files("TCGA_Dataset/miRNA/gdc_download_20241115_102256.015891/", recursive = TRUE)
  stop(paste("No files found with the specified pattern. Available files are:",
             paste(available_files, collapse = "\n")))
}

# Process each file using purrr::map for cleaner code
exp <- purrr::map(count_files, ~ {
  file_path <- paste0("TCGA_Dataset/miRNA/gdc_download_20241115_102256.015891/", .x)

  # Read and process the file
  read.table(file_path, sep = "\t", header = TRUE) %>%
    dplyr::select(miRNA_region = 6, reads_per_million_miRNA_mapped = 4) %>%
    dplyr::group_by(miRNA_region) %>%
    dplyr::summarise(reads_per_million_miRNA_mapped = sum(reads_per_million_miRNA_mapped, na.rm = TRUE))
})

# Merge all data frames in the list 'exp' by 'miRNA_region'
m <- Reduce(function(x, y) merge(x, y, by = 'miRNA_region', all = TRUE), exp)

# Replace any NA values with 0
m[is.na(m)] <- 0

# Set 'miRNA_region' as row names for the merged data
exp <- column_to_rownames(m, var = "miRNA_region")

# Remove the last 2 rows, assuming they contain invalid data
exp <- exp[-((nrow(exp)-2):nrow(exp)),]

# Check for "mature," in row names and remove the prefix
if (any(str_detect(rownames(exp), "mature,"))) {
  rownames(exp) <- str_remove(rownames(exp), "mature,")
}

# Retrieve miRNA information using MIMAT IDs from the miRBaseVersions.db database
mh <- select(miRBaseVersions.db,
             keys = rownames(exp),  # Use MIMAT ID as key
             keytype = "MIMAT",  # Specify the key type as MIMAT
             columns = c("ACCESSION", "NAME", "VERSION"))  # Retrieve ACCESSION, NAME, and VERSION columns

# Filter for miRNAs with version "21"
mh <- mh[mh$VERSION == "21",]

# Match the row names (MIMAT IDs) with the ACCESSION in the database and reorder
mh <- mh[match(rownames(exp), mh$ACCESSION),]

# Update row names with the miRNA names from the database
rownames(exp) <- mh$NAME

# Read metadata in JSON format, containing information about the samples
meta <- jsonlite::fromJSON("TCGA_Dataset/miRNA/metadata.cart.2024-11-15.json")

# Create a mapping data frame between file names and IDs
file2id <- data.frame(file_name = meta$file_id, ID = sapply(meta$associated_entities, function(x) x$entity_submitter_id))

# Extract file names from the paths and match them with the file2id data frame
count_files2 <- stringr::str_split(count_files, "/", simplify = TRUE)[, 1]
file2id <- file2id[match(count_files2, file2id$file_name),]

# Assign the IDs to the columns of 'exp' (expression data)
colnames(exp) <- file2id$ID

# Filter out rows where the sum of values > 0 is not greater than 100 (miRNAs that are not expressed in enough samples)
k <- apply(exp, 1, function(x) sum(x > 0) > 100)

# Keep only the rows that passed the filter
exp <- exp[k,]

# Convert the data frame to a matrix
exp <- as.matrix(exp)

# Create a new data frame with miRNA IDs and expression data
exp1 <- data.frame(ID = rownames(exp), exp)

# Replace '.' with '-' in column names for better readability
colnames(exp1) <- gsub('[.]', '-', colnames(exp1))

# Write the processed expression data to a text file
write.table(exp1, 'miRNA.RPM.txt', sep = "\t", quote = FALSE, row.names = FALSE)
