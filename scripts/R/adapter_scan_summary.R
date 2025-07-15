#!/usr/bin/env Rscript

# adapter_scan_summary.R

library(readr)
library(dplyr)
library(ggplot2)
library(optparse)
library(patchwork)
library(glue)

fill_color <- "black"

# Usage:
# Rscript scripts/R/adapter_scan_summary.R \
# -i out/E105_3C2_P5_maxima_txg_kapahifi/ont/adapter_scan.tsv \
# -s out/E105_3C2_P5_maxima_txg_kapahifi/ont/misc_logs/1a_adapter_scan_summary.csv \
# -p out/E105_3C2_P5_maxima_txg_kapahifi/ont/misc_logs/1a_adapter_scan_summary.pdf,out/E105_3C2_P5_maxima_txg_kapahifi/ont/misc_logs/1a_adapter_scan_summary.png \
# 1> out/E105_3C2_P5_maxima_txg_kapahifi/ont/misc_logs/1a_adapter_scan_summary.log \
# 2> out/E105_3C2_P5_maxima_txg_kapahifi/ont/misc_logs/1a_adapter_scan_summary.err

# Start time
start_time <- Sys.time()
cat(
  glue("Script started at: {format(start_time, '%Y-%m-%d %H:%M:%S')}"),
  "\n"
)
cat("\n")

# Function to load the data
load_data <- function(data_file, n_max = Inf) {
  # Read the data with specified column types
  df <- read_tsv(
    data_file,
    col_types = cols(
      readlen = col_double(),
      start = col_double(),
      end = col_double(),
      strand = col_character(),
      lab = col_character(),
      adapter_config = col_character(),
      read_id = col_character()
    ),
    na = c("", "NA", "None"),
    n_max = n_max
  )

  return(df)
}

# Define command line options
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Input TSV file path",
    metavar = "FILE"
  ),
  make_option(
    c("-s", "--summary"),
    type = "character",
    default = "adapter_scan_summary.csv",
    help = "Output file for summary CSV [default= %default]",
    metavar = "FILE"
  ),
  make_option(
    c("-p", "--plots"),
    type = "character",
    default = "adapter_scan_summary.pdf",
    help = "Output file(s) for combined plots, comma-separated [default= %default]",
    metavar = "FILES"
  ),
  make_option(
    c("-v", "--verbose"),
    type = "logical",
    default = TRUE,
    help = "Whether or not to print updates [default= %default]",
    metavar = "VERBOSE"
  ),
  make_option(
    c("-n", "--nrows"),
    type = "integer",
    default = Inf,
    help = "Number of rows to read from the input file [default= %default]",
    metavar = "NROWS"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

verbose <- opt$verbose

# Log the parameters passed in
cat("Parameters passed in:\n")
cat(glue("Input file:     {opt$input}"), "\n")
cat(glue("Summary file:   {opt$summary}"), "\n")
cat(glue("Plot files:     {opt$plots}"), "\n")
cat(glue("Verbose:        {opt$verbose}"), "\n")
cat(glue("Number of rows: {opt$nrows}"), "\n")
cat("\n")

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Input file path must be provided. Use -h for help.")
}

# Ensure n_max is properly handled
n_max <- ifelse(is.na(opt$nrows) || opt$nrows == 0, Inf, opt$nrows)

# Read the TSV file
if (verbose) {
  cat("Reading in tsv...\n")
}
data <- load_data(opt$input,
  n_max = n_max
)
cat(glue("Number of reads analyzed: {nrow(data)}"), "\n")
cat("\n")

data$lab <- outFactor <- factor(
  x = data$lab,
  levels = c(
    "full_len",
    "adapter1_single",
    "adapter1_double",
    "adapter2_single",
    "adapter2_double",
    "other",
    "other_ambiguous",
    "no_adapters"
  )
)

# Generate summary grouped by 'lab'
if (verbose) {
  cat("Computing stats...\n")
}

# Full read length ("readlen" is just the insert length plus primers)
data <- data %>%
  mutate(readlength = end + start)

data <- data %>%
  mutate(normStart = start / readlength)

summary <- data %>%
  group_by(lab) %>%
  summarise(
    count = n(),
    avg_insertlen = round(mean(readlen), 2),
    avg_readlen = round(mean(readlength), 2),
    med_insertlen = round(median(readlen), 2),
    med_readlen = round(median(readlength), 2),
    min_readlen = round(min(readlength), 2),
    max_readlen = round(max(readlength), 2),
    avg_start = round(mean(start), 2),
    avg_normStart = round(mean(normStart), 2)
  )
cat(
  glue("Summary statistics computed for {nrow(summary)} groups"),
  "\n"
)

# Write summary to CSV file
if (verbose) {
  cat("Writing summary...\n")
}
write_csv(summary, opt$summary)

# Generate plots
if (verbose) {
  cat("Generating plots...\n")
}

# Filter data for p2 to within 3 standard deviations
readlen_mean <- mean(data$readlength, na.rm = TRUE)
readlen_sd <- sd(data$readlength, na.rm = TRUE)
data_filtered <- data[abs(data$readlength - readlen_mean) <= 3 * readlen_sd, ]

# Filter data for p3 to within 3 standard deviations (excluding no_adapters)
labels_to_exclude <- c(
  "adapter2_single",
  "adapter2_double",
  "other",
  "other_ambiguous",
  "no_adapters"
)
data_no_adapters <- data[!data$lab %in% labels_to_exclude, ]
start_mean <- mean(data_no_adapters$start, na.rm = TRUE)
start_sd <- sd(data_no_adapters$start, na.rm = TRUE)
data_p3_filtered <- data_no_adapters[abs(data_no_adapters$start - start_mean) <= 3 * start_sd, ]

# Log number of points removed
cat(glue("Read length filtering: {nrow(data) - nrow(data_filtered)} points removed (kept {nrow(data_filtered)}/{nrow(data)})"), "\n")
cat(glue("Start position filtering: {nrow(data_no_adapters) - nrow(data_p3_filtered)} points removed (kept {nrow(data_p3_filtered)}/{nrow(data_no_adapters)})"), "\n")

current_theme <- theme(
  plot.title = element_text(hjust = 0.5, face = "bold"),
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.title.x = element_blank()
)

# 1. Read count by lab
p1 <- ggplot(summary, aes(x = lab, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Read Count by Type",
    x = "Read Type",
    y = "Count"
  ) +
  theme_minimal() +
  current_theme

# 2. Average read length by lab
p2 <- ggplot(
  data_filtered,
  aes(x = lab, y = readlength)
) +
  geom_violin(
    fill = fill_color,
    scale = "width"
  ) +
  labs(
    title = "Read Length Distribution by Type",
    x = "Read Type",
    y = "Read Length"
  ) +
  theme_minimal() +
  current_theme

# 3. Average start position for insert
p3 <- ggplot(
  data_p3_filtered,
  aes(x = lab, y = start)
) +
  geom_violin(
    fill = fill_color,
    scale = "width"
  ) +
  labs(
    title = "Start Position of Adapter 1",
    x = "Read Type",
    y = "Start Position"
  ) +
  theme_minimal() +
  current_theme

# 4. Average start position for insert, normalized to read length
p4 <- ggplot(
  data_p3_filtered,
  aes(x = lab, y = normStart)
) +
  geom_violin(
    fill = fill_color,
    scale = "width"
  ) +
  labs(
    title = "Start Position of Adapter 1, relative to read length",
    x = "Read Type",
    y = "Normalized Start Position"
  ) +
  theme_minimal() +
  current_theme

# 5. Read length distribution by lab
# p5 <- ggplot(
#   data,
#   aes(x = readlen, fill = lab)) +
#   geom_histogram(bins = 30, position = "dodge") +
#   labs(
#     title = "Read Length Distribution by Lab",
#     x = "Read Length",
#     y = "Count"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "bottom")


# combine into a single plot
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Adapter Scan Summary",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
  )

# Save combined plot
if (verbose) {
  cat("Saving combined plot...\n")
}
plot_files <- strsplit(opt$plots, ",")[[1]]
for (plot_file in plot_files) {
  device <- tools::file_ext(plot_file)
  if (device == "") {
    stop("File extension must be provided to infer the device.")
  }
  ggsave(
    plot_file,
    combined_plot,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    device = device
  )
}

cat("Analysis complete. Summary and combined plot have been saved.\n")

# End time and run time
end_time <- Sys.time()
cat("\n")
cat(
  glue("Script ended at: {format(end_time, '%Y-%m-%d %H:%M:%S')}"),
  "\n"
)
cat(
  glue(
    "Total run time: {difftime(end_time, start_time, units = 'secs')} seconds"
  ),
  "\n"
)
