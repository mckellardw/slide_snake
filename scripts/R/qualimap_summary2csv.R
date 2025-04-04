suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(optparse)
})

# Function to parse the Qualimap report
parse_qualimap_report <- function(qualimap_report_txt) {
  lines <- read_lines(qualimap_report_txt)
  data <- list()
  current_section <- NULL
  
  for (line in lines) {
    line <- str_trim(line)
    if (line == "") next  # Skip empty lines
    if (str_starts(line, ">>>>>>>")) {
      current_section <- str_remove(line, ">>>>>>> ") %>% str_trim()
      next
    }
    if (str_detect(line, "=")) {
      parts <- str_split(line, "=", n = 2)[[1]]
      key <- str_trim(parts[1])
      value <- str_trim(parts[2])
      if (!is.null(current_section)) {
        key <- paste(current_section, key, sep = "_")
      }
      value <- str_remove_all(value, ",")
      if (str_detect(value, "\\(")) {
        value <- str_split(value, "\\(", n = 2)[[1]][1] %>% str_trim()
      }
      if (str_detect(value, "%")) {
        value <- as.numeric(str_remove(value, "%"))
      } else if (str_detect(value, "^[0-9.]+$")) {
        value <- as.numeric(value)
      }
      data[[key]] <- value
    }
  }
  
  # Convert to data frame and clean column names
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  colnames(df) <- colnames(df) %>%
    str_replace_all("\\.\\.", "_") %>%    # Replace double periods with underscores
    str_replace_all("\\.", "_") %>%    # Replace periods with underscores
    str_replace_all("[^A-Za-z0-9_]", "")  %>% # Remove special characters
    str_replace_all(" ", "_")       # Replace spaces with underscores
  
  return(df)
}

# Function to parse command-line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(
      c("--input"), 
      type = "character", 
      help = "Path to the Qualimap report text file.", 
      metavar = "FILE"
    ),
    make_option(
      c("--output"), 
      type = "character", 
      help = "Path to the output CSV file.", 
      metavar = "FILE"
    ),
    make_option(
      c("--verbose"), 
      type = "logical", 
      default = TRUE, 
      help = "Print comments during execution. Default is TRUE."
    )
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(opt_parser)
  
  return(opts)
}

# Main function
main <- function() {
  args <- parse_arguments()
  
  if (args$verbose) message("Parsing Qualimap report...")
  data <- parse_qualimap_report(args$input)
  if (args$verbose) message(glue("   Parsed data: {ncol(data)} columns"))
  
  if (args$verbose) message("Writing parsed data to CSV...")
  write_csv(data, args$output)
  
  if (args$verbose) message("Execution completed.")
}

# Run the script
if (!interactive()) {
  main()
}
