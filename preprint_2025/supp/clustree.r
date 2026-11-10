# Load required libraries
library(clustree)
library(readr)
library(ggplot2)

# Set the input and output file paths
input_file <- “/n/eddy_lab/Lab/mckinley/cagri_output/20250523-AmhrGsD_Final/Clustree/20250523-AmhrGsD_Final_DataForClustree.csv”
output_file <- “/n/eddy_lab/Lab/mckinley/cagri_output/20250523-AmhrGsD_Final/Clustree/20250523-AmhrGsD_Final-Clustree.png”

# Read the data
data <- read_csv(input_file)

# Generate the clustree plot
clustree_plot <- clustree(data, prefix = “leiden_scVI_“)

# Save the plot
ggsave(output_file, clustree_plot, width = 10, height = 8, dpi = 300)