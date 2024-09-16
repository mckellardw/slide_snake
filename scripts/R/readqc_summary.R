#!/usr/bin/env Rscript

# Load necessary libraries
library(readr)
library(ggplot2)
library(patchwork)
library(optparse)
library(scales)
library(ggforce)

# Define the options for command-line arguments
option_list <- list(
    make_option(
        c("-f", "--file"),
        type = "character", default = NULL, help = "Input data file name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character", 
        default = "combined_plot.png",
        help = "Output file name for the combined plot [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-s", "--sample"),
        type = "integer", 
        default = 1000000, 
        help = "Number of reads to downsample [default= %default]",
        metavar = "integer"
    )
)

## Example input:
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

# Parse the command-line arguments --------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input file is provided --------------------------------------------
if (is.null(opt$file)) {
    print_help(opt_parser)
    stop("Input file must be supplied.\n", call. = FALSE)
}

big.text <- 10
small.text <- 8

custom_theme <- theme_minimal() +
 theme(
    # Centered titles
    plot.title = element_text(hjust = 0.5, face="bold", size=big.text),
    axis.title.x = element_text(hjust = 0.5, size=small.text, face = "bold"),
    axis.title.y = element_text(hjust = 0.5, size=small.text, face = "bold"),
    
    # Black axes
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black", size=small.text),
    
    # Scientific notation for y-axis
    axis.text.y = element_text(angle = 0, hjust = 1)
) 

# Function to load the data
load_data <- function(data_file, n_reads){    
    # Read the data
    df <- read_tsv(
        data_file,
        na = c("", "NA", "None"),
        n_max=n_reads
    )

    # Downsample
    if(n_reads < nrow(df)){
        df <- df[sample(nrow(df), n_reads), ]
    }
    
    return(df)
}


# Function to read data and create plots --------------------------------------------
create_plots <- function(df, data_file, out_file) {

    # Summary histograms (left column)
    plot.gc <- ggplot(df, aes(x = GC_Percent)) +
        geom_histogram(binwidth = 1, fill = "red", alpha = 0.5) +
        custom_theme +
        scale_y_continuous(labels = label_scientific()) +
        ggtitle("GC Percentage Distribution")

    plot.readLength <- ggplot(df, aes(x = Read_Length)) +
        geom_histogram(
            binwidth = 50, fill = "blue", alpha = 0.5
        ) +
        custom_theme +
        lims(
            x=c(0,10000)
        )+
        scale_y_continuous(
            # trans="log10",
            labels = label_scientific()
        ) +
        ggtitle("Read Length Distribution")

    plot.homoploymer <- ggplot(
            df, aes(x = Longest_Homopolymer/Read_Length)
        ) +
        geom_histogram(
            aes(
                fill=Homopolymer_Base
            ),
            binwidth = 0.01
            # alpha = 0.5
        ) +
        custom_theme +
        # scale_x_continuous(trans='log2') +
        xlim(0,0.25)+
        scale_y_continuous(labels = label_scientific()) +
        ggtitle("Longest Homopolymer Fraction")


    # Summary scatter plots (right column)
    scatter.gc <- ggplot(
        df, 
        aes(
            x = Read_Length,
            y = GC_Percent
            )
        ) +
        # geom_hex(binwidth = c(NA,1)) +
        geom_point(
            size=0.1,
            alpha=0.2
        )+
        custom_theme +
        # scale_y_continuous(labels = label_scientific()) +
        ggtitle("Read Length vs. GC Percentage")

    
    plot.readLength.zoom <- ggplot(
            df, 
            aes(x = Read_Length)
        ) +
        geom_histogram(
            binwidth = 1, fill = "blue", alpha = 0.5
        ) +
        custom_theme +
        lims(
            x=c(0,5000)
        )+
        scale_y_continuous(
            # trans="log10",
            labels = label_scientific()
        ) +
        ggtitle("Small RNAs")

    
    vln.homopolymer <- ggplot(
            df, 
            aes(
                x = Homopolymer_Base,
                y = Longest_Homopolymer,
                fill = Homopolymer_Base
            )
        ) +
        geom_violin() +
        custom_theme +
        lims(y=c(0,200))+
        # scale_y_continuous(trans="log10") +
        # scale_y_log10()+
        ggtitle("Longest Homopolymer Length")

    # Re-size Read1 plots..
    if(grepl("R1", data_file)){
        scatter.gc <- scatter.gc + xlim(0,1000)
        
        plot.readLength.zoom <- plot.readLength.zoom + xlim(0,100)
    }

    # Combine plots using patchwork
    combined_plot <- wrap_plots(
        wrap_plots(
            plot.gc,
            plot.readLength,
            plot.homoploymer,
            ncol=1
        ),
         wrap_plots(
            scatter.gc,
            plot.readLength.zoom,
            vln.homopolymer,
            ncol=1
        ),
        ncol=2
     )

    # Infer the device from the filename extension
    device <- tools::file_ext(out_file)
    device <- gsub("^\\.", "", device)  # Remove the dot from the extension

    # Save the combined plot
    ggsave(
        out_file, 
        plot = combined_plot, 
        device = device,
        width = 240,
        height = 240,
        units = "mm",#c("in", "cm", "mm", "px"),
        dpi = 300,
        create.dir = T
    )
}

# Run the function with the provided arguments --------------------------------------------
df <- load_data(
    opt$file, 
    opt$sample
)

create_plots(
    df, 
    opt$file,
    opt$out
)
