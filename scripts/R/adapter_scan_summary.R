#!/usr/bin/env Rscript

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(optparse)
library(patchwork) 


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
  make_option(c("-s", "--summary"), type="character", default="adapter_scan_summary.csv", 
              help="Output file for summary CSV [default= %default]", metavar="FILE"),
  make_option(c("-p", "--plots"), type="character", default="adapter_scan_summary.pdf", 
              help="Output file for combined plots [default= %default]", metavar="FILE"),
  make_option(c("-d", "--device"), type="character", default="pdf", 
              help="Output device for plots (pdf, png, svg) [default= %default]", metavar="DEVICE"),
  make_option(c("-v", "--verbose"), type="logical", default=TRUE, 
              help="Whether or not to print updates [default= %default]", metavar="VERBOSE"),
  make_option(c("-n", "--nrows"), type="integer", default=Inf, 
              help="Number of rows to read from the input file [default= %default]", metavar="NROWS")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

verbose <- opt$verbose

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Input file path must be provided. Use -h for help.")
}

# Read the TSV file
if(verbose){message("Reading in tsv...")}
data <- load_data(opt$input, n_max=opt$nrows)

data$lab <- outFactor <- factor(
  x = data$lab, 
  levels = c(
    "full_len",
    "single_adapter1",
    "double_adapter1",
    "single_adapter2",
    "double_adapter2",
    "other",
    "no_adapters"
  )
)

# Generate summary grouped by 'lab'
if(verbose){message("Computing stats...")}

# Full read length ("readlen" is just the insert length plus primers)
data <- data %>%
  mutate(readlength = end + start)

data <- data %>%
  mutate(normStart = start/readlength)

summary <- data %>%
  group_by(lab) %>%
  summarise(
    count = n(),
    avg_insertlen = mean(readlen),
    avg_readlen = mean(readlength),
    # med_insertlen = median(readlen),
    # med_readlen = median(readlength),
    # min_readlen = min(readlength),
    # max_readlen = max(readlength),
    avg_start = mean(start),
    avg_normStart = mean(normStart),
    pct_stranded = sum(stranded) / n() * 100
  )
  
# Write summary to CSV file
if(verbose){message("Writing summary...")}
write_csv(summary, opt$summary)

# Generate plots
if(verbose){message("Generating plots...")}

current_theme <- theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(size=14, angle = 45, hjust = 1)
  )

# 1. Read count by lab
p1 <- ggplot(summary, aes(x = lab, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Read Count by Lab", x = "Lab", y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 2. Average read length by lab
p2 <- ggplot(data, aes(x = lab, y = readlength)) +
  geom_violin(fill = "darkgreen") +
  labs(
    title = "Read Length Distribution by Lab", 
    x = "Read Type", 
    y = "Read Length"
  ) +
  theme_minimal() +
  current_theme

# 3. Average start position
p3 <- ggplot(data, aes(x = lab, y = start)) +
  geom_violin(fill = "orange") +
  labs(
    title = "Start Position of Adapter 1", 
    x = "Read Type", 
    y = "Start Position"
  ) +
  theme_minimal() +
  current_theme

# 4. Average start position, normalized to read length
p4 <- ggplot(data, aes(x = lab, y = normStart)) +
  geom_violin(fill = "purple") +
  labs(
    title = "Start Position of Adapter 1, relative to read length", 
    x = "Read Type", 
    y = "Normalized Start Position"
  ) +
  theme_minimal() +
  current_theme

# 5. Read length distribution by lab
p5 <- ggplot(data, aes(x = readlen, fill = lab)) +
  geom_histogram(bins = 30, position = "dodge") +
  labs(
    title = "Read Length Distribution by Lab", 
    x = "Read Length", 
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine all plots using patchwork
# combined_plot <- (p1 + p2) / (p3 + p4) / p5 +
#   plot_layout(heights = c(1, 1, 1.5)) +
#   plot_annotation(title = "Adapter Scan Summary",
#                   theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face="bold")))

combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Adapter Scan Summary",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face="bold")))

# Save combined plot
if(verbose){message("Saving combined plot...")}
ggsave(
  opt$plots, 
  combined_plot, 
  width = 12, 
  height = 12, 
  units = "in", 
  dpi = 300,
  device = opt$device
  # create.dir = TRUE
)

message("Analysis complete. Summary and combined plot have been saved.\n")