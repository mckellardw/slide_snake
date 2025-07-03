# cache_tx_to_seurat.R

library(Seurat)
library(argparse)
library(glue)

# Example usage:
# Rscript scripts/r/cache_mtx_to_seurat.R \
#   --mat_in input_matrix.mtx \
#   --feat_in input_features.tsv \
#   --bc_in input_barcodes.txt \
#   --bc_map input_spatial_map.tsv \
#   --seurat_out output_seurat.rds

parse_args <- function() {
  parser <-
    ArgumentParser(description = "Process spatial transcriptomics data with Seurat.")
  parser$add_argument("--mat_in", required = TRUE, help = "Input count matrix file (mtx format)")
  parser$add_argument("--feat_in", required = TRUE, help = "Input feature file (tsv format)")
  parser$add_argument("--bc_in", required = TRUE, help = "Input barcode file (txt format)")
  parser$add_argument("--bc_map", required = TRUE, help = "Input spatial map file (tsv format)")
  parser$add_argument("--seurat_out", required = TRUE, help = "Output Seurat object file (rds format)")
  parser$add_argument("--feat_col",
                      type = "integer",
                      default = 1,
                      help = "Feature column index in the feature file (default: 1)")
  parser$add_argument("--transpose",
                      type = "logical",
                      default = FALSE,
                      help = "Transpose count matrix? (default: FALSE)")
  return(parser$parse_args())
}

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

main <-
  function(mat_in,
           feat_in,
           bc_in,
           bc_map,
           seurat_out,
           feat_col,
           transpose) {
    # Load count matrix using Seurat
    cat(timestamp(), " - Loading count matrix...\n")
    counts <-
      Read10X(
        data.dir = dirname(mat_in),
        gene.column = feat_col,
        cell.column = 1
      )
    cat(timestamp(), " - Done.\n\n")
    
    if (transpose) {
      counts <- t(counts)
    }
    
    # Load features
    features <- read.table(feat_in, sep = "\t", header = FALSE)
    rownames(counts) <- features$V1
    
    # Load barcodes
    barcodes <- read.table(bc_in, sep = "\t", header = FALSE)
    colnames(counts) <- barcodes$V1
    
    # Print the number of features and cells loaded
    cat(timestamp(), " - Number of features loaded: ", nrow(counts), "\n")
    cat(timestamp(), " - Number of cells loaded:    ", ncol(counts), "\n")
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = counts)
    
    # Load spatial coordinates
    cat(timestamp(), " - Loading spatial barcode map...\n")
    barcode_map <-
      read.table(
        bc_map,
        sep = "\t",
        header = FALSE,
        col.names = c("barcode", "x", "y")
      )
    cat(timestamp(), " - Done.\n")
    
    rownames(barcode_map) <- barcode_map$barcode
    
    # Match cells between Seurat object and spatial data
    common_cells <-
      intersect(colnames(seurat_obj), rownames(barcode_map))
    if (length(common_cells) == 0) {
      stop("No matching barcodes between data and spatial map")
    } else if (length(common_cells) == length(colnames(seurat_obj))) {
      cat(timestamp(),
          " - No barcodes missing between data and spatial map\n")
    } else {
      cat(timestamp(),
          " - Subsetting Seurat object and spatial data to the common cells\n")
      seurat_obj <- subset(seurat_obj, cells = common_cells)
    }
    
    # Warning if barcode_map has more rows than Seurat object has columns
    if (nrow(barcode_map) > ncol(seurat_obj)) {
      cat(
        timestamp(),
        glue(
          " - WARNING: barcode map has more rows [{nrow(barcode_map)}] than Seurat object has columns [{ncol(seurat_obj)}]\n"
        )
      )
    }
    
    # Subset Seurat object and spatial data to the common cells
    embeddings <- as.matrix(barcode_map[common_cells, c("x", "y")])
    colnames(embeddings) <-
      paste0("spatial_", seq_len(ncol(embeddings)))
    seurat_obj[["spatial"]] <-
      CreateDimReducObject(embeddings = embeddings, key = "spatial_")
    
    # Save Seurat object
    saveRDS(seurat_obj, file = seurat_out)
    cat(timestamp(), " - Seurat object saved to", seurat_out, "\n")
  }

args <- parse_args()
cat(
  "Matrix file:                  ",
  args$mat_in,
  "\n",
  "Features/genes file:          ",
  args$feat_in,
  "\n",
  "Barcodes file:                ",
  args$bc_in,
  "\n",
  "Barcode map file:             ",
  args$bc_map,
  "\n",
  "Output Seurat object file:    ",
  args$seurat_out,
  "\n",
  "Feature column index:         ",
  args$feat_col,
  "\n",
  "Transpose matrix:             ",
  args$transpose,
  "\n\n",
  sep = ""
)
main(
  args$mat_in,
  args$feat_in,
  args$bc_in,
  args$bc_map,
  args$seurat_out,
  args$feat_col,
  args$transpose
)
