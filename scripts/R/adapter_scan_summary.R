#!/usr/bin/env Rscript

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(optparse)

verbose <- T

# Function to load the data
load_data <- function(data_file, n_max=Inf){    
    # Read the data
    df <- read_tsv(
        data_file,
        na = c("", "NA", "None"),
        n_max=n_max
    )

    # # Downsample
    # if(n_max < nrow(df)){
    #     df <- df[sample(nrow(df), n_max), ]
    # }
    
    return(df)
}

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input TSV file path", metavar="FILE"),
  make_option(c("-s", "--summary"), type="character", default="nanopore_summary.csv", 
              help="Output file for summary CSV [default= %default]", metavar="FILE"),
  make_option(c("-p", "--plots"), type="character", default="nanopore_plot", 
              help="Output prefix for plot files [default= %default]", metavar="PREFIX"),
  make_option(c("-d", "--device"), type="character", default="png", 
              help="Output device for plots (png, pdf, svg) [default= %default]", metavar="DEVICE")
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
# Generate summary grouped by 'lab'
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

# Function to save plot with specified device
save_plot <- function(plot, filename, device = opt$device, width = 10, height = 6) {
  if (device == "png") {
    ggsave(filename, plot, device = "png", width = width, height = height, dpi = 300)
  } else if (device == "pdf") {
    ggsave(filename, plot, device = "pdf", width = width, height = height)
  } else if (device == "svg") {
    ggsave(filename, plot, device = "svg", width = width, height = height)
  } else {
    stop("Unsupported device. Use png, pdf, or svg.")
  }
}

# Generate plots
if(verbose){print("Plotting summary...")}
# 1. Read count by lab
p1 <- ggplot(summary, aes(x = lab, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Read Count by Lab", x = "Lab", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p1, paste0(opt$plots, "_read_count.", opt$device))

# 2. Average read length by lab
p2 <- ggplot(summary, aes(x = lab, y = avg_readlen)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  theme_minimal() +
  labs(title = "Average Read Length by Lab", x = "Lab", y = "Average Read Length") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p2, paste0(opt$plots, "_avg_readlen.", opt$device))

# 3. Percentage of full-length reads by lab
p3 <- ggplot(summary, aes(x = lab, y = pct_full_length)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_minimal() +
  labs(title = "Percentage of Full-Length Reads by Lab", x = "Lab", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p3, paste0(opt$plots, "_pct_full_length.", opt$device))

# 4. Percentage of stranded reads by lab
p4 <- ggplot(summary, aes(x = lab, y = pct_stranded)) +
  geom_bar(stat = "identity", fill = "purple") +
  theme_minimal() +
  labs(title = "Percentage of Stranded Reads by Lab", x = "Lab", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p4, paste0(opt$plots, "_pct_stranded.", opt$device))

# 5. Read length distribution by lab
p5 <- ggplot(data, aes(x = readlen, fill = lab)) +
  geom_histogram(bins = 30, position = "dodge") +
  theme_minimal() +
  labs(title = "Read Length Distribution by Lab", x = "Read Length", y = "Count") +
  theme(legend.position = "bottom")

if(verbose){print("Saving plots...")}
save_plot(p5, paste0(opt$plots, "_readlen_dist.", opt$device))

cat("Analysis complete. Summary and plots have been saved.\n")