#!/usr/bin/env python3
"""
Process adapter insertion for barcode sequences.

This script handles insertion of adapter sequences into barcodes
for different spatial RNA-seq technologies.

Example usage:
python scripts/py/process_adapter_insertion.py \
    --bc-map-file input/map.txt \
    --bc-adapter output/whitelist_adapter.txt \
    --bc-adapter-r1 output/whitelist_adapter_r1.txt \
    --bc-adapter-map output/map_adapter.txt \
    --bc-adapter-r1-map output/map_adapter_R1.txt \
    --adapter ACGTACGTACGT \
    --bc-primer ATCGATCG \
    --bc-lengths "8 6" \
    --recipes seeker \
    --recipes-to-split "seeker,miST,decoder"
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
from datetime import datetime


def currentTime():
    """Return current time as a formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def print_log(msg):
    """Print a log message with the current time."""
    print(f"{currentTime()} | {msg}")


def print_error(msg):
    """Print an error message with the current time."""
    print(f"{currentTime()} | ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def log_message(message, level="info", verbosity=2):
    """Print message to console based on verbosity level."""
    if level == "error" and verbosity >= 1:
        print_error(message)
    elif level == "info" and verbosity >= 2:
        print_log(message)
    elif level == "debug" and verbosity >= 3:
        print_log(f"DEBUG: {message}")


def get_recipe_type(recipes):
    """Determine the recipe type based on the recipes list."""
    recipe_str = "".join(recipes)
    
    # Check for specific recipe patterns
    if "seeker" in recipe_str:
        return "seeker"
    elif "decoder" in recipe_str:
        return "decoder"
    elif "miST" in recipe_str:
        return "miST"
    else:
        return "default"


def should_process_adapters(recipes, recipes_to_split):
    """Check if this recipe needs adapter insertion."""
    recipe_str = "".join(recipes)
    return any(recipe in recipe_str for recipe in recipes_to_split)


def get_barcode_split_position(bc_lengths, recipes):
    """Get the position to split barcodes based on recipe type and lengths."""
    recipe_type = get_recipe_type(recipes)
    
    # Get split position from recipe sheet or use defaults
    if bc_lengths and len(bc_lengths) > 0:
        # Parse lengths from recipe sheet (e.g., "8 6" -> split at position 8)
        if isinstance(bc_lengths[0], str):
            lengths = [int(x) for x in bc_lengths[0].split()]
            split_pos = lengths[0] if lengths else 8
        else:
            split_pos = bc_lengths[0]
    else:
        # Default split positions based on recipe type
        if recipe_type == "seeker":
            split_pos = 8
        elif recipe_type == "decoder":
            split_pos = 8
        elif recipe_type == "miST":
            split_pos = 10
        else:
            split_pos = 8
    
    return split_pos


def split_and_insert_adapter(bc_map, adapter, bc_primer, bc_lengths, recipes, verbosity=2):
    """Split barcodes and insert adapter sequences."""
    recipe_type = get_recipe_type(recipes)
    split_pos = get_barcode_split_position(bc_lengths, recipes)
    
    if recipe_type == "decoder":
        # TODO: Implement decoder-specific logic
        log_message(f"Adapter insertion for {recipe_type} not yet implemented", level="info", verbosity=verbosity)
        return None, None
    
    # Split barcodes at the specified position
    bc_1 = [bc[:split_pos] for bc in list(bc_map[0])]
    bc_2 = [bc[split_pos:] for bc in list(bc_map[0])]
    
    # Create adapter-inserted sequences
    bc_adapter = [f"{bc1}{adapter}{bc2}" for bc1, bc2 in zip(bc_1, bc_2)]
    bc_adapter_r1 = [f"{bc_primer}{bc1}{adapter}{bc2}" for bc1, bc2 in zip(bc_1, bc_2)]
    
    return bc_adapter, bc_adapter_r1


def process_adapters(bc_map, recipes, adapter, bc_primer, bc_lengths, recipes_to_split, output_files, verbosity=2):
    """Process adapter insertion and generate output files."""
    if not should_process_adapters(recipes, recipes_to_split):
        # Create empty files if adapter insertion is not needed
        for output_file in output_files.values():
            Path(output_file).touch()
        log_message("No adapter insertion needed for this recipe", level="info", verbosity=verbosity)
        return
    
    bc_adapter, bc_adapter_r1 = split_and_insert_adapter(bc_map, adapter, bc_primer, bc_lengths, recipes, verbosity)
    
    if bc_adapter is None or bc_adapter_r1 is None:
        # Create empty files if processing failed
        for output_file in output_files.values():
            Path(output_file).touch()
        log_message("Adapter insertion failed, created empty files", level="info", verbosity=verbosity)
        return
    
    # Create adapter maps with additional columns from original map
    bc_adapter_map = pd.DataFrame(list(zip(bc_adapter, bc_map[1], bc_map[2])))
    bc_adapter_r1_map = pd.DataFrame(list(zip(bc_adapter_r1, bc_map[1], bc_map[2])))
    
    # Save all files
    pd.Series(bc_adapter).to_csv(
        output_files["BC_ADAPTER"], sep="\n", header=False, index=False
    )
    pd.Series(bc_adapter_r1).to_csv(
        output_files["BC_ADAPTER_R1"], sep="\n", header=False, index=False
    )
    bc_adapter_map.to_csv(
        output_files["BC_ADAPTER_MAP"], sep="\t", header=False, index=False
    )
    bc_adapter_r1_map.to_csv(
        output_files["BC_ADAPTER_R1_MAP"], sep="\t", header=False, index=False
    )
    
    log_message(f"Successfully processed adapter insertion for {get_recipe_type(recipes)}", level="info", verbosity=verbosity)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process adapter insertion for barcode sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process Seeker barcodes with adapter insertion
  python process_adapter_insertion.py \\
    --bc-map-file input/map.txt \\
    --bc-adapter output/whitelist_adapter.txt \\
    --bc-adapter-r1 output/whitelist_adapter_r1.txt \\
    --bc-adapter-map output/map_adapter.txt \\
    --bc-adapter-r1-map output/map_adapter_R1.txt \\
    --adapter ACGTACGTACGT \\
    --bc-primer ATCGATCG \\
    --bc-lengths "8 6" \\
    --recipes seeker \\
    --recipes-to-split "seeker,miST,decoder"
        """
    )
    
    parser.add_argument(
        "--bc-map-file",
        required=True,
        help="Path to input barcode map file (TSV format)"
    )
    
    parser.add_argument(
        "--bc-adapter",
        required=True,
        help="Path to output adapter-inserted barcode whitelist file"
    )
    
    parser.add_argument(
        "--bc-adapter-r1",
        required=True,
        help="Path to output R1 adapter-inserted barcode whitelist file"
    )
    
    parser.add_argument(
        "--bc-adapter-map",
        required=True,
        help="Path to output adapter-inserted barcode map file"
    )
    
    parser.add_argument(
        "--bc-adapter-r1-map",
        required=True,
        help="Path to output R1 adapter-inserted barcode map file"
    )
    
    parser.add_argument(
        "--adapter",
        required=True,
        help="Internal adapter sequence to insert"
    )
    
    parser.add_argument(
        "--bc-primer",
        required=True,
        help="Forward primer sequence for R1 processing"
    )
    
    parser.add_argument(
        "--bc-lengths",
        default="8 6",
        help="Space-separated barcode lengths for splitting (e.g., '8 6')"
    )
    
    parser.add_argument(
        "--recipes",
        default="",
        help="Comma-separated list of recipes being processed"
    )
    
    parser.add_argument(
        "--recipes-to-split",
        default="seeker,miST,decoder",
        help="Comma-separated list of recipes that need adapter insertion"
    )
    
    parser.add_argument(
        "--verbosity",
        type=int,
        default=2,
        help="Verbosity level (0=quiet, 1=error, 2=info, 3=debug)"
    )
    
    return parser.parse_args()


def main():
    """Main function to process adapter insertion based on command line arguments."""
    args = parse_args()
    
    try:
        # Create output directories if they don't exist
        output_files = {
            "BC_ADAPTER": args.bc_adapter,
            "BC_ADAPTER_R1": args.bc_adapter_r1,
            "BC_ADAPTER_MAP": args.bc_adapter_map,
            "BC_ADAPTER_R1_MAP": args.bc_adapter_r1_map,
        }
        
        # Create directories for output files
        for output_file in output_files.values():
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        # Parse parameters
        bc_lengths = [args.bc_lengths.strip()] if args.bc_lengths.strip() else []
        
        # Parse recipes
        recipes = []
        if args.recipes.strip():
            recipes = [r.strip() for r in args.recipes.split(",")]
        
        # Parse recipes to split
        recipes_to_split = []
        if args.recipes_to_split.strip():
            recipes_to_split = [r.strip() for r in args.recipes_to_split.split(",")]
        
        # Load barcode map
        bc_map = pd.read_csv(args.bc_map_file, sep="\t", header=None)
        if args.verbosity >= 2:
            print_log(f"Loaded barcode map with {len(bc_map)} entries for adapter insertion")
        
        # Process adapters using functional approach
        process_adapters(
            bc_map=bc_map,
            recipes=recipes,
            adapter=args.adapter,
            bc_primer=args.bc_primer,
            bc_lengths=bc_lengths,
            recipes_to_split=recipes_to_split,
            output_files=output_files,
            verbosity=args.verbosity
        )
        
        if args.verbosity >= 2:
            print_log(f"Successfully processed adapter insertion for recipes: {', '.join(recipes) if recipes else 'default'}")
        
    except Exception as e:
        print_error(f"Error processing adapter insertion: {e}")
        if args.verbosity >= 3:
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()
