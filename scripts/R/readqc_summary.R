#!/usr/bin/env Rscript

# Load necessary libraries
library(readr)
library(ggplot2)
library(patchwork)
library(optparse)
library(scales)

# Define the options for command-line arguments
option_list <- list(
    make_option(
        c("-f", "--file"),
        type = "character", default = NULL, help = "Input data file name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character", default = "combined_plot.png", help = "Output file name for the combined plot [default= %default]",
        metavar = "character"
    )
)

# Read_ID Read_Length     GC_Percent      First_Base      Last_Base       Longest_Homopolymer     Homopolymer_Base
# @67ae57cf-e0ec-438f-a07b-45f7dbe312e0_0 14      42.86   A       C       2       A
# @d918b6df-d6f4-4431-acf4-11fb944b9b40_0 155     57.42   G       T       4       T
# @f49ccee4-1d37-4191-9948-8ef8063190cf_0 205     52.2    G       T       5       A
# @d93965e7-51eb-4d22-a4c0-10a5e4c8bf99_0 156     58.33   G       T       4       C
# @d8bc4554-51e6-4192-9c40-6e5ae2b7820f_0 116     59.48   G       T       3       A
# @6b90a0af-a09e-413f-a4c9-230f1c8050f8_0 103     51.46   A       T       4       T
# @320e06cf-fbd1-40ea-9953-b2fcbb523e25_0 267     50.19   A       T       4       A
# @d95e1228-0ea6-4fc6-ba70-b38040715b30_0 92      50.0    G       C       4       C
# @8eec0f77-78d0-4a62-a64b-0b82914d86b1_0 155     57.42   G       T       4       T
# @e1a5289d-263d-448d-b34d-f22f73500e14_0 78      53.85   G       G       2       T
# @4f43c08a-8547-4798-8a55-1c1bfac4cadd_0 178     60.11   C       G       5       G
# @7a7086a9-a2b2-4636-b830-d5df76afdb78_0 210     30.95   G       T       5       T
# @180ad7da-a3a9-4f08-9ea8-ae51c3ba51c3_0 154     56.49   G       T       4       A
# @fa0d831f-6cea-4b1d-9d16-2f2751de8061_0 102     58.82   A       G       5       A
# @30aaf6c4-bd60-4932-a79d-4aff18565fbc_0 73      30.14   A       T       4       G


# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input file is provided
if (is.null(opt$file)) {
    print_help(opt_parser)
    stop("Input file must be supplied.\n", call. = FALSE)
}

custom_theme <- theme_minimal() +
 theme(
    # Centered titles
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5),
    
    # Black axes
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    
    # Scientific notation for y-axis
    axis.text.y = element_text(angle = 0, hjust = 1)
) 

# Function to read data and create plots
create_plots <- function(data_file, out_file) {
    # Read the data
    data <- read_tsv(
        data_file,
        na = c("", "NA", "None")
    )

    # Create summary plots
    plot.gc <- ggplot(data, aes(x = GC_Percent)) +
        geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5) +
        custom_theme +
        scale_y_continuous(labels = label_scientific()) +
        ggtitle("GC Percentage Distribution")

    plot.readLength <- ggplot(data, aes(x = Read_Length)) +
        geom_histogram(binwidth = 50, fill = "blue", alpha = 0.5) +
        custom_theme +
        scale_y_continuous(
            trans="log10",
            labels = label_scientific()
        ) +
        ggtitle("Read Length Distribution")

    plot.homoploymer <- ggplot(
            data, aes(x = Longest_Homopolymer)
        ) +
        geom_histogram(
            aes(
                fill=Homopolymer_Base
            ),
            binwidth = 1, 
            # fill = "red",
            alpha = 0.5
        ) +
        custom_theme +
        scale_x_continuous(trans='log2') +
        scale_y_continuous(labels = label_scientific()) +
        ggtitle("Longest Homopolymer Length")

    # Combine plots using patchwork
    combined_plot <- plot.gc + 
        plot.readLength + 
        plot.homoploymer + 
        plot_layout(ncol = 1)

    # Infer the device from the filename extension
    device <- tools::file_ext(out_file)
    device <- gsub("^\\.", "", device)  # Remove the dot from the extension

    # Save the combined plot
    ggsave(
        out_file, 
        plot = combined_plot, 
        device = device
    )
}

# Run the function with the provided arguments
create_plots(opt$file, opt$out)
