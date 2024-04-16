#!/usr/bin/env Rscript

library(jsonlite)
library(ggplot2)
library(optparse)

# Function to plot a summary of the Cutadapt JSON report
plot_cutadapt_summary <- function(json_file, out_file) {
 # Read the JSON file
 data <- fromJSON(json_file)
  
 # Extract relevant data for plotting
 summary_data_r1 <- data.frame(
    Reads = c("Total", "With Adapters"),
    Count = c(data$total_reads, data$reads_with_adapters)
 )
for(name in names(data$adapters_read1))

 # Plot the summary
 p <- ggplot(summary_data, aes(x = Reads, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Cutadapt Summary", x = "Reads", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
 # Infer the device from the filename extension
 device <- tools::file_ext(out_file)
 device <- gsub("^\\.", "", device) # Remove the dot from the extension
  
 # Save the plot to the specified file
 ggsave(out_file, plot = p, device = device)
}

# Define the options for command-line arguments
option_list <- list(
 make_option(c("-f", "--file"), type="character", default=NULL,
              help="Path to the Cutadapt JSON report file", metavar="character"),
 make_option(c("-o", "--out"), type="character", default="cutadapt_summary.png",
              help="Output file name for the plot [default= %default]", metavar="character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if the input file is provided
if (is.null(opt$file)) {
 print_help(opt_parser)
 stop("Input file must be supplied.\n", call.=FALSE)
}

# Run the function with the provided arguments
plot_cutadapt_summary(opt$file, opt$out)
