#!/usr/bin/env Rscript

# Suppress messages and warnings when loading libraries
suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(argparse)
})

# Define command-line arguments
parser <- ArgumentParser(description = "Generate QC plots for single-cell RNA-seq data")
parser$add_argument("--matrix", type = "character", required = TRUE, help = "Path to the matrix file")
parser$add_argument("--genes", type = "character", required = TRUE, help = "Path to the genes file")
parser$add_argument("--cells", type = "character", required = TRUE, help = "Path to the cells file")
parser$add_argument("--output", type = "character", required = TRUE, help = "Path to save the output file")

# Parse arguments
opt <- parser$parse_args()

# Infer the device from the output file extension
output_extension <- tools::file_ext(opt$output)
if (output_extension == "") {
    stop("Output file must have a valid extension (e.g., .png, .pdf).")
}

# Check if the file type is supported by ggsave
supported_devices <- c("png", "pdf", "jpeg", "tiff", "bmp", "svg")
if (!(output_extension %in% supported_devices)) {
    stop(paste("Unsupported file type:", output_extension, ". Supported types are:", paste(supported_devices, collapse = ", ")))
}

# Load the matrix
matrix <- ReadMtx(
    mtx = opt$matrix,
    features = opt$genes,
    cells = opt$cells
)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(
    counts = matrix,
    project = "SeuratProject",
    assay="RNA"
)

# Custom theme for plots
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

# Extract metadata for plotting
metadata <- seurat_obj@meta.data

# Generate QC plots using ggplot
# nFeat_vln <- ggplot(metadata, aes(x = "", y = nFeature_RNA)) +
#     geom_violin(
#         fill = "#CE4C51", 
#         color = "black"
#     ) +
#     geom_jitter(
#         width = 0.1, 
#         height = 0, 
#         size = 0.5, 
#         color = "black"
#     ) +
#     custom_theme +
#     ggtitle("Number of Features") +
#     theme(
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank()
#     )

# nCount_vln <- ggplot(
#     metadata, 
#     aes(x = "", y = nCount_RNA)
#     ) +
#     geom_violin(
#         fill = "#52A0DC", 
#         color = "black"
#     ) +
#     geom_jitter(
#         width = 0.1, 
#         height = 0, 
#         size = 0.5, 
#         color = "black"
#     ) +
#     custom_theme +
#     ggtitle("Number of Counts") +
#     theme(
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank()
#     )

# Add knee plot
sorted_counts <- sort(seurat_obj$nCount_RNA, decreasing = TRUE)
knee_plot <- ggplot(
    data.frame(
        Index = seq_along(sorted_counts), 
        Counts = sorted_counts
    ), 
    aes(x = Index, y = Counts)
    ) +
    geom_line() +
    custom_theme +
    ggtitle("Knee Plot") +
    xlab("Cell Index") +
    ylab("# UMIs") +
    scale_y_continuous(labels = scales::scientific) +
    scale_x_continuous(labels = scales::scientific)

# Add scatter plot
scatter_plot <- ggplot(
        metadata, 
        aes(x = nCount_RNA, y = nFeature_RNA)
    ) +
    geom_point(
        alpha = 0.4, 
        color = "black", 
        size = 0.2
    ) +
    geom_smooth(
        method = "gam",
        formula = y ~ s(x, bs = "cs"),
        color = "#CE4C51",
        linewidth = 0.5, 
        se = TRUE
    ) +
    scale_x_continuous(
        labels = scales::scientific
    ) +
    scale_y_continuous(
        labels = scales::scientific
    ) +
    custom_theme +
    ggtitle("# UMIs vs. # Features") +
    xlab("# UMIs") +
    ylab("# Features")

# Combine plots using patchwork
combined_plot <- (knee_plot | scatter_plot) + 
    plot_layout(guides = "collect")

# Save the plot
ggsave(
    filename = opt$output,
    plot = combined_plot, 
    width = 150, 
    height = 75, 
    units = "mm",
    dpi = 300,
    device = output_extension
)
