library("XML")
library("methods")

# Set working directory
setwd("D:/R/Rcode/LUSC_Analysis")

# Define the directory containing XML files
dir <- "D:/R/Rcode/LUSC_Analysis/TCGA_Dataset/clinical/gdc_download_20241113_195625.183908"

# List all XML files in the directory (including subdirectories)
all_files <- list.files(path = dir, pattern = '*.xml$', recursive = TRUE, full.names = TRUE)

# Parse XML files and convert them to data frames
cl <- lapply(all_files, function(x) {
  # Parse the XML, extract the root node, and convert the second child node to a data frame
  xmldataframe <- xmlToDataFrame(xmlRoot(xmlParse(file = x))[2])
  return(xmldataframe)
})

# Combine all data frames into one
clinical <- do.call(rbind, cl)

# Write the combined data frame to a text file
write.table(clinical, file = "clinical.txt", sep = "\t", quote = FALSE, row.names = FALSE)
