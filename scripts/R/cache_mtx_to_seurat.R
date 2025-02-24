library(Seurat)
library(argparse)

# Example usage:
# Rscript scripts/r/cache_mtx_to_seurat.R \
#   --mat_in input_matrix.mtx \
#   --feat_in input_features.tsv \
#   --bc_in input_barcodes.txt \
#   --bc_map input_spatial_map.tsv \
#   --seurat_out output_seurat.rds

parse_args <- function() {
  parser <- ArgumentParser(description = "Process spatial transcriptomics data with Seurat.")
  parser$add_argument("--mat_in", required = TRUE, help = "Input count matrix file (mtx format)")
  parser$add_argument("--feat_in", required = TRUE, help = "Input feature file (tsv format)")
  parser$add_argument("--bc_in", required = TRUE, help = "Input barcode file (txt format)")
  parser$add_argument("--bc_map", required = TRUE, help = "Input spatial map file (tsv format)")
  parser$add_argument("--seurat_out", required = TRUE, help = "Output Seurat object file (rds format)")
  parser$add_argument("--feat_col", type = "integer", default = 1, help = "Feature column index in the feature file (default: 1)")
  parser$add_argument("--transpose", type = "logical", default = FALSE, help = "Transpose count matrix? (default: FALSE)")
  return(parser$parse_args())
}

main <- function(mat_in, feat_in, bc_in, bc_map, seurat_out, feat_col, transpose) {
  # Load count matrix using Seurat
  counts <- Read10X(data.dir = dirname(mat_in), gene.column = feat_col, cell.column = 1)
  
  if (transpose) {
    counts <- t(counts)
  }
  
  # Load features
  features <- read.table(feat_in, sep = "\t", header = FALSE)
  rownames(counts) <- features$V1
  
  # Load barcodes
  barcodes <- read.table(bc_in, sep = "\t", header = FALSE)
  colnames(counts) <- barcodes$V1
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # Load spatial coordinates
  spatial_data <- read.table(bc_map, sep = "\t", header = FALSE, col.names = c("barcode", "x", "y"))
  rownames(spatial_data) <- spatial_data$barcode
  
  # Match cells between Seurat object and spatial data
  common_cells <- intersect(colnames(seurat_obj), rownames(spatial_data))
  if(length(common_cells) == 0) stop("No matching barcodes between data and spatial map")
  
  # Subset Seurat object and spatial data to the common cells
  seurat_obj <- subset(seurat_obj, cells = common_cells)
  embeddings <- as.matrix(spatial_data[common_cells, c("x", "y")])
  colnames(embeddings) <- paste0("spatial_", seq_len(ncol(embeddings)))
  seurat_obj[["spatial"]] <- CreateDimReducObject(embeddings = embeddings, key = "spatial_")
  
  # Save Seurat object
  saveRDS(seurat_obj, file = seurat_out)
  cat("Seurat object saved to", seurat_out, "\n")
}

args <- parse_args()
cat(
  "Matrix file:                  ", args$mat_in, "\n",
  "Features/genes file:          ", args$feat_in, "\n",
  "Barcodes file:                ", args$bc_in, "\n",
  "Barcode map file:             ", args$bc_map, "\n",
  "Output Seurat object file:    ", args$seurat_out, "\n",
  "Feature column index:         ", args$feat_col, "\n",
  "Transpose matrix:             ", args$transpose, "\n",
  sep=""
)
main(args$mat_in, args$feat_in, args$bc_in, args$bc_map, args$seurat_out, args$feat_col, args$transpose)
