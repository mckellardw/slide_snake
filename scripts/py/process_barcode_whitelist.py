#!/usr/bin/env python3
"""
Process barcode whitelists for different spatial RNA-seq technologies.

This script handles barcode splitting and formatting for various recipes
based on configuration from the recipe sheet.

Example usage:
python scripts/py/process_barcode_whitelist.py \
    --bc-map-file resources/chromium_whitelists/3M-february-2018_map.txt \
    --bc-us-map output/sample/bc/map_underscore.txt \
    --bc-1 output/sample/bc/whitelist_1.txt \
    --bc-2 output/sample/bc/whitelist_2.txt \
    --bc-uniq-1 output/sample/bc/whitelist_uniq_1.txt \
    --bc-uniq-2 output/sample/bc/whitelist_uniq_2.txt \
    --bc-us output/sample/bc/whitelist_underscore.txt \
    --log-file output/sample/bc/info.log \
    --bc-lengths "8 6" \
    --bc-concat true \
    --recipes seeker
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


def get_barcode_lengths(bc_lengths, recipes):
    """Get barcode lengths for splitting."""
    recipe_type = get_recipe_type(recipes)

    # Get lengths from recipe sheet or use defaults
    if bc_lengths and len(bc_lengths) > 0:
        # Parse lengths from recipe sheet (e.g., "8 6" -> [8, 6])
        if isinstance(bc_lengths[0], str):
            lengths = [int(x) for x in bc_lengths[0].split()]
        else:
            lengths = bc_lengths
    else:
        # Default lengths based on recipe type
        if recipe_type == "seeker":
            lengths = [8, 6]
        elif recipe_type == "decoder":
            lengths = [8, 8]  # Default for decoder
        elif recipe_type == "miST":
            lengths = [10, 10]
        else:
            lengths = []

    return lengths


def should_concatenate_barcodes(bc_concat, recipes):
    """Check if barcodes should be concatenated."""
    if bc_concat and len(bc_concat) > 0:
        concat_val = bc_concat[0]
        if isinstance(concat_val, str):
            return concat_val.lower() == "true"
        else:
            return bool(concat_val)

    # Default behavior based on recipe type
    recipe_type = get_recipe_type(recipes)
    return recipe_type in ["seeker", "decoder", "miST"]


def split_barcodes(bc_map, bc_lengths, recipes, verbosity=2):
    """Split barcodes based on recipe configuration."""
    recipe_type = get_recipe_type(recipes)
    lengths = get_barcode_lengths(bc_lengths, recipes)

    if not lengths or len(lengths) < 2:
        log_message(
            f"No barcode splitting needed for {recipe_type}",
            level="debug",
            verbosity=verbosity,
        )
        return None, None, None

    bc_1_length = lengths[0]

    # Split barcodes
    bc_1 = [bc[:bc_1_length] for bc in list(bc_map[0])]
    bc_2 = [bc[bc_1_length:] for bc in list(bc_map[0])]
    bc_us = [f"{bc[:bc_1_length]}_{bc[bc_1_length:]}" for bc in list(bc_map[0])]

    log_message(
        f"Split barcodes for {recipe_type}: {bc_1_length}+{len(bc_2[0]) if bc_2 else 0}bp",
        level="info",
        verbosity=verbosity,
    )

    return bc_1, bc_2, bc_us


def create_default_files(bc_map, output_files, verbosity=2):
    """Create default files for recipes that don't need barcode splitting."""
    # Copy original map with underscore (no split, so just copy first column)
    bc_map.to_csv(output_files["BC_US_MAP"], header=False, index=False, sep="\t")

    # Always write unique barcode lists for every chemistry
    uniq_1 = pd.Series(bc_map[0]).drop_duplicates()
    uniq_1.to_csv(output_files["BC_UNIQ_1"], header=False, index=False)

    # For BC_UNIQ_2, if a second column exists, use it; otherwise, do not write the file
    if bc_map.shape[1] > 1:
        uniq_2 = pd.Series(bc_map[1]).drop_duplicates()
        uniq_2.to_csv(output_files["BC_UNIQ_2"], header=False, index=False)
    else:
        # Remove file if it exists and would be empty
        try:
            Path(output_files["BC_UNIQ_2"]).unlink(missing_ok=True)
        except Exception:
            pass

    # Only write BC_1 and BC_2 if they would have content (i.e., splitting is performed)
    # For default/no-split, do not write these files

    # Create simple whitelist (underscore version is just the first column)
    pd.Series(bc_map[0]).to_csv(output_files["BC_US"], header=False, index=False)

    log_message(
        "Processed default barcodes (no splitting, unique lists written)", level="info", verbosity=verbosity
    )


def process_barcodes(bc_map, recipes, bc_lengths, bc_concat, output_files, verbosity=2):
    """Process barcodes and generate output files."""
    recipe_type = get_recipe_type(recipes)
    should_split = (
        should_concatenate_barcodes(bc_concat, recipes)
        and len(get_barcode_lengths(bc_lengths, recipes)) >= 2
    )

    if should_split:
        bc_1, bc_2, bc_us = split_barcodes(bc_map, bc_lengths, recipes, verbosity)

        if bc_1 is not None and bc_2 is not None and bc_us is not None:
            # Create underscore map
            bc_us_map = bc_map.copy()
            bc_us_map.iloc[:, 0] = bc_us

            # Save all files
            bc_us_map.to_csv(
                output_files["BC_US_MAP"], header=False, index=False, sep="\t"
            )
            pd.Series(bc_1).to_csv(output_files["BC_1"], header=False, index=False)
            pd.Series(bc_2).to_csv(output_files["BC_2"], header=False, index=False)
            pd.Series(list(set(bc_1))).to_csv(
                output_files["BC_UNIQ_1"], header=False, index=False
            )
            pd.Series(list(set(bc_2))).to_csv(
                output_files["BC_UNIQ_2"], header=False, index=False
            )
            pd.Series(bc_us).to_csv(output_files["BC_US"], header=False, index=False)

            log_message(
                f"Processed {recipe_type} barcodes with splitting",
                level="info",
                verbosity=verbosity,
            )
        else:
            create_default_files(bc_map, output_files, verbosity)
    else:
        create_default_files(bc_map, output_files, verbosity)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process barcode whitelists for different spatial RNA-seq technologies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process Seeker barcodes with splitting
  python process_barcode_whitelist.py \\
    --bc-map-file input/map.txt \\
    --bc-us-map output/map_underscore.txt \\
    --bc-1 output/whitelist_1.txt \\
    --bc-2 output/whitelist_2.txt \\
    --bc-uniq-1 output/whitelist_uniq_1.txt \\
    --bc-uniq-2 output/whitelist_uniq_2.txt \\
    --bc-us output/whitelist_underscore.txt \\
    --log-file output/info.log \\
    --bc-lengths "8 6" \\
    --bc-concat true \\
    --recipes seeker

  # Process Visium barcodes (no splitting)
  python process_barcode_whitelist.py \\
    --bc-map-file input/map.txt \\
    --bc-us-map output/map_underscore.txt \\
    --bc-1 output/whitelist_1.txt \\
    --bc-2 output/whitelist_2.txt \\
    --bc-uniq-1 output/whitelist_uniq_1.txt \\
    --bc-uniq-2 output/whitelist_uniq_2.txt \\
    --bc-us output/whitelist_underscore.txt \\
    --log-file output/info.log \\
    --bc-lengths "" \\
    --bc-concat false \\
    --recipes visium
        """,
    )

    parser.add_argument(
        "--bc-map-file",
        required=True,
        help="Path to input barcode map file (TSV format)",
    )

    parser.add_argument(
        "--bc-us-map", required=True, help="Path to output underscore barcode map file"
    )

    parser.add_argument(
        "--bc-1", required=True, help="Path to output barcode 1 whitelist file"
    )

    parser.add_argument(
        "--bc-2", required=True, help="Path to output barcode 2 whitelist file"
    )

    parser.add_argument(
        "--bc-uniq-1",
        required=True,
        help="Path to output unique barcode 1 whitelist file",
    )

    parser.add_argument(
        "--bc-uniq-2",
        required=True,
        help="Path to output unique barcode 2 whitelist file",
    )

    parser.add_argument(
        "--bc-us",
        required=True,
        help="Path to output underscore barcode whitelist file",
    )

    parser.add_argument(
        "--log-file",
        help="Path to output log file (deprecated - now prints to console)",
    )

    parser.add_argument(
        "--bc-lengths",
        default="",
        help="Space-separated barcode lengths for splitting (e.g., '8 6')",
    )

    parser.add_argument(
        "--bc-concat",
        default="false",
        help="Whether to concatenate/split barcodes (true/false)",
    )

    parser.add_argument(
        "--recipes", default="", help="Comma-separated list of recipes being processed"
    )

    parser.add_argument(
        "--verbosity",
        type=int,
        default=2,
        help="Verbosity level (0=quiet, 1=error, 2=info, 3=debug)",
    )

    return parser.parse_args()


def main():
    """Main function to process barcodes based on command line arguments."""
    args = parse_args()

    try:
        # Create output directories if they don't exist
        output_files = {
            "BC_US_MAP": args.bc_us_map,
            "BC_1": args.bc_1,
            "BC_2": args.bc_2,
            "BC_UNIQ_1": args.bc_uniq_1,
            "BC_UNIQ_2": args.bc_uniq_2,
            "BC_US": args.bc_us,
        }

        # Create directories for output files
        for output_file in output_files.values():
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)

        # Create log file directory if specified (for backward compatibility)
        if args.log_file:
            Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)

        # Parse parameters
        bc_lengths = []
        if args.bc_lengths.strip():
            bc_lengths = [args.bc_lengths.strip()]

        bc_concat = [args.bc_concat.lower() == "true"]

        # Parse recipes
        recipes = []
        if args.recipes.strip():
            recipes = [r.strip() for r in args.recipes.split(",")]

        # Load barcode map
        bc_map = pd.read_csv(args.bc_map_file, sep="\t", header=None)
        if args.verbosity >= 2:
            print_log(f"Loaded barcode map with {len(bc_map)} entries")

        # Process barcodes using functional approach
        process_barcodes(
            bc_map=bc_map,
            recipes=recipes,
            bc_lengths=bc_lengths,
            bc_concat=bc_concat,
            output_files=output_files,
            verbosity=args.verbosity,
        )

        if args.verbosity >= 2:
            print_log(
                f"Successfully processed barcodes for recipes: {', '.join(recipes) if recipes else 'default'}"
            )

        # Write to log file if specified (for backward compatibility)
        if args.log_file:
            with open(args.log_file, "w") as f:
                f.write(
                    f"Processed barcodes for recipes: {', '.join(recipes) if recipes else 'default'}\n"
                )

    except Exception as e:
        print_error(f"Error processing barcodes: {e}")
        if args.verbosity >= 3:
            import traceback

            traceback.print_exc()


if __name__ == "__main__":
    main()
