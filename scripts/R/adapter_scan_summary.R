#!/usr/bin/env Rscript

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(optparse)
library(patchwork)  # New library for combining plots

verbose <- TRUE

# Function to load the data
load_data <- function(data_file, n_max=Inf){    
    # Read the data
    df <- read_tsv(
        data_file,
        na = c("", "NA", "None"),
        n_max=n_max
    )
    
    return(df)
}

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input TSV file path", metavar="FILE"),
  make_option(c("-s", "--summary"), type="character", default="nanopore_summary.csv", 
              help="Output file for summary CSV [default= %default]", metavar="FILE"),
  make_option(c("-p", "--plots"), type="character", default="nanopore_plots.pdf", 
              help="Output file for combined plots [default= %default]", metavar="FILE"),
  make_option(c("-d", "--device"), type="character", default="pdf", 
              help="Output device for plots (pdf, png, svg) [default= %default]", metavar="DEVICE")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Input file path must be provided. Use -h for help.")
}

# Read the TSV file
if(verbose){print("Reading in tsv...")}
data <- load_data(opt$input)

# Generate summary grouped by 'lab'
if(verbose){print("Computing stats...")}
summary <- data %>%
  group_by(lab) %>%
  summarise(
    count = n(),
    avg_readlen = mean(readlen),
    med_readlen = median(readlen),
    min_readlen = min(readlen),
    max_readlen = max(readlen),
    pct_full_length = sum(fl) / n() * 100,
    pct_stranded = sum(stranded) / n() * 100
  )
  
# Write summary to CSV file
if(verbose){print("Writing summary...")}
write_csv(summary, opt$summary)

# Generate plots
if(verbose){print("Generating plots...")}

# 1. Read count by lab
p1 <- ggplot(summary, aes(x = lab, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Read Count by Lab", x = "Lab", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Average read length by lab
p2 <- ggplot(summary, aes(x = lab, y = avg_readlen)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Average Read Length by Lab", x = "Lab", y = "Average Read Length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Percentage of full-length reads by lab
p3 <- ggplot(summary, aes(x = lab, y = pct_full_length)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_minimal() +
  labs(title = "Percentage of Full-Length Reads by Lab", x = "Lab", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Percentage of stranded reads by lab
p4 <- ggplot(summary, aes(x = lab, y = pct_stranded)) +
  geom_bar(stat = "identity", fill = "purple") +
  theme_minimal() +
  labs(title = "Percentage of Stranded Reads by Lab", x = "Lab", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Read length distribution by lab
p5 <- ggplot(data, aes(x = readlen, fill = lab)) +
  geom_histogram(bins = 30, position = "dodge") +
  theme_minimal() +
  labs(title = "Read Length Distribution by Lab", x = "Read Length", y = "Count") +
  theme(legend.position = "bottom")

# Combine all plots using patchwork
combined_plot <- (p1 + p2) / (p3 + p4) / p5 +
  plot_layout(heights = c(1, 1, 1.5)) +
  plot_annotation(title = "Nanopore Sequencing Analysis",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))

# Save combined plot
if(verbose){print("Saving combined plot...")}
ggsave(opt$plots, combined_plot, width = 15, height = 20, units = "in", device = opt$device)

cat("Analysis complete. Summary and combined plot have been saved.\n")