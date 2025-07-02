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


def get_barcode_lengths(bc_lengths):
    """Get barcode lengths for splitting."""
    if bc_lengths and len(bc_lengths) > 0:
        # Parse lengths from recipe sheet (e.g., "8 6" -> [8, 6])
        if isinstance(bc_lengths[0], str):
            lengths = [int(x) for x in bc_lengths[0].split()]
        else:
            lengths = bc_lengths
    else:
        lengths = []
    return lengths


def split_barcodes(bc_map, bc_lengths):
    """Split barcodes based on provided lengths."""
    lengths = get_barcode_lengths(bc_lengths)

    if not lengths or len(lengths) < 2:
        print_log("No barcode splitting performed. Barcodes will be written as-is.")
        return None, None, None

    bc_1_length = lengths[0]
    bc_2_length = lengths[1] if len(lengths) > 1 else 0
    print_log(
        f"Splitting barcodes: first {bc_1_length} bases to BC_1, next {bc_2_length} bases to BC_2."
    )

    # Split barcodes
    bc_1 = [bc[:bc_1_length] for bc in list(bc_map[0])]
    bc_2 = [bc[bc_1_length:] for bc in list(bc_map[0])]
    bc_us = [f"{bc[:bc_1_length]}_{bc[bc_1_length:]}" for bc in list(bc_map[0])]

    print_log(f"Split barcodes: {bc_1_length}+{len(bc_2[0]) if bc_2 else 0}bp")

    return bc_1, bc_2, bc_us


def create_default_files(bc_map, output_files):
    """Create default files for recipes that don't need barcode splitting."""
    # Copy original map with underscore (no split, so just copy first column)
    bc_map.to_csv(output_files["BC_US_MAP"], header=False, index=False, sep="\t")
    print_log(
        f"Saved barcode map with underscores to: {output_files['BC_US_MAP']}\n  ({len(bc_map):,} barcodes)"
    )

    # Always write unique barcode lists for every chemistry
    uniq_1 = pd.Series(bc_map[0]).drop_duplicates()
    uniq_1.to_csv(output_files["BC_UNIQ_1"], header=False, index=False)
    percent_uniq_1 = 100.0 * len(uniq_1) / len(bc_map) if len(bc_map) > 0 else 0.0
    print_log(
        f"Saved unique barcode 1 whitelist to: {output_files['BC_UNIQ_1']}\n"
        f"  {len(uniq_1):,} unique barcodes\n"
        f"  {percent_uniq_1:.2f}% unique\n"
    )

    # For BC_UNIQ_2, if a second column exists, use it; otherwise, do not write the file
    if bc_map.shape[1] > 1:
        uniq_2 = pd.Series(bc_map[1]).drop_duplicates()
        uniq_2.to_csv(output_files["BC_UNIQ_2"], header=False, index=False)
        percent_uniq_2 = 100.0 * len(uniq_2) / len(bc_map) if len(bc_map) > 0 else 0.0
        print_log(
            f"Saved unique barcode 2 whitelist to: {output_files['BC_UNIQ_2']}\n"
            f"  {len(uniq_2):,} unique barcodes\n"
            f"  {percent_uniq_2:.2f}% unique\n"
        )
    else:
        # Remove file if it exists and would be empty
        try:
            Path(output_files["BC_UNIQ_2"]).unlink(missing_ok=True)
        except Exception:
            pass

    # For default/no-split, also write BC_1 and BC_2 as the first and second columns if present
    pd.Series(bc_map[0]).to_csv(output_files["BC_1"], header=False, index=False)
    print_log(
        f"Saved barcode 1 whitelist to: {output_files['BC_1']}\n"
        f"  {len(bc_map[0]):,} barcodes\n"
    )
    if bc_map.shape[1] > 1:
        pd.Series(bc_map[1]).to_csv(output_files["BC_2"], header=False, index=False)
        print_log(
            f"Saved barcode 2 whitelist to: {output_files['BC_2']}\n"
            f"  {len(bc_map[1]):,} barcodes\n"
        )

    # Create simple whitelist (underscore version is just the first column)
    pd.Series(bc_map[0]).to_csv(output_files["BC_US"], header=False, index=False)
    print_log(
        f"Saved underscore barcode whitelist to: {output_files['BC_US']}\n"
        f"  {len(bc_map):,} barcodes\n"
    )


def process_barcodes(bc_map, bc_lengths, output_files):
    """Process barcodes and generate output files."""
    lengths = get_barcode_lengths(bc_lengths)
    should_split = False
    if lengths and len(lengths) >= 2:
        should_split = True
    if should_split:
        bc_1, bc_2, bc_us = split_barcodes(bc_map, bc_lengths)
        if bc_1 is not None and bc_2 is not None and bc_us is not None:
            # Create underscore map
            bc_us_map = bc_map.copy()
            bc_us_map.iloc[:, 0] = bc_us

            # Save all files
            bc_us_map.to_csv(
                output_files["BC_US_MAP"], header=False, index=False, sep="\t"
            )
            print_log(
                f"Saved barcode map with underscores to: {output_files['BC_US_MAP']}\n"
                f"  ({len(bc_us_map):,} barcodes)"
            )

            pd.Series(bc_1).to_csv(output_files["BC_1"], header=False, index=False)
            print_log(
                f"Saved barcode 1 whitelist to: {output_files['BC_1']}\n"
                f"  {len(bc_1):,} barcodes\n"
            )

            pd.Series(bc_2).to_csv(output_files["BC_2"], header=False, index=False)
            print_log(
                f"Saved barcode 2 whitelist to: {output_files['BC_2']}\n"
                f"  {len(bc_2):,} barcodes\n"
            )

            uniq_1 = list(set(bc_1))
            pd.Series(uniq_1).to_csv(
                output_files["BC_UNIQ_1"], header=False, index=False
            )
            percent_uniq_1 = 100.0 * len(uniq_1) / len(bc_1) if len(bc_1) > 0 else 0.0
            print_log(
                f"Saved unique barcode 1 whitelist to: {output_files['BC_UNIQ_1']}\n"
                f"  {len(uniq_1):,} unique barcodes\n"
                f"  {percent_uniq_1:.2f}% unique\n"
            )

            uniq_2 = list(set(bc_2))
            pd.Series(uniq_2).to_csv(
                output_files["BC_UNIQ_2"], header=False, index=False
            )
            percent_uniq_2 = 100.0 * len(uniq_2) / len(bc_2) if len(bc_2) > 0 else 0.0
            print_log(
                f"Saved unique barcode 2 whitelist to: {output_files['BC_UNIQ_2']}\n"
                f"  {len(uniq_2):,} unique barcodes\n"
                f"  {percent_uniq_2:.2f}% unique\n"
            )

            pd.Series(bc_us).to_csv(output_files["BC_US"], header=False, index=False)
            print_log(
                f"Saved underscore barcode whitelist to: {output_files['BC_US']}\n"
                f"  {len(bc_us):,} barcodes\n"
            )

            print_log(f"Processed barcodes with splitting\n")
        else:
            create_default_files(bc_map, output_files)
    else:
        create_default_files(bc_map, output_files)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process barcode whitelists for spatial RNA-seq barcodes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process barcodes with splitting
  python process_barcode_whitelist.py \
    --bc-map-file input/map.txt \
    --bc-us-map output/map_underscore.txt \
    --bc-1 output/whitelist_1.txt \
    --bc-2 output/whitelist_2.txt \
    --bc-uniq-1 output/whitelist_uniq_1.txt \
    --bc-uniq-2 output/whitelist_uniq_2.txt \
    --bc-us output/whitelist_underscore.txt \
    --bc-lengths "8 6"

  # Process barcodes (no splitting)
  python process_barcode_whitelist.py \
    --bc-map-file input/map.txt \
    --bc-us-map output/map_underscore.txt \
    --bc-1 output/whitelist_1.txt \
    --bc-2 output/whitelist_2.txt \
    --bc-uniq-1 output/whitelist_uniq_1.txt \
    --bc-uniq-2 output/whitelist_uniq_2.txt \
    --bc-us output/whitelist_underscore.txt \
    --bc-lengths ""
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
        "--bc-lengths",
        default="",
        help="Space-separated barcode lengths for splitting (e.g., '8 6')",
    )

    return parser.parse_args()


def main():
    """Main function to process barcodes based on command line arguments."""
    args = parse_args()

    # Print all input and output parameters at the start
    print(
        f"Input barcode map:   {args.bc_map_file}\n"
        f"Output BC_US_MAP:    {args.bc_us_map}\n"
        f"Output BC_1:         {args.bc_1}\n"
        f"Output BC_2:         {args.bc_2}\n"
        f"Output BC_UNIQ_1:    {args.bc_uniq_1}\n"
        f"Output BC_UNIQ_2:    {args.bc_uniq_2}\n"
        f"Output BC_US:        {args.bc_us}\n"
        f"BC Lengths:          {args.bc_lengths}\n"
        f"\n"
    )

    # Parameter checks
    if not args.bc_map_file or not Path(args.bc_map_file).exists():
        print_error(f"Input barcode map file does not exist: {args.bc_map_file}")
    for out_file in [
        args.bc_us_map,
        args.bc_1,
        args.bc_2,
        args.bc_uniq_1,
        args.bc_uniq_2,
        args.bc_us,
    ]:
        if not out_file:
            print_error(f"Missing output file path for: {out_file}")
    if not args.bc_lengths:
        print_log("Warning: --bc-lengths is empty. No splitting will be performed.")
    else:
        try:
            _ = [int(x) for x in str(args.bc_lengths).split()]
        except Exception:
            print_error(f"Could not parse --bc-lengths: {args.bc_lengths}")

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

        # Parse parameters
        bc_lengths = []
        if args.bc_lengths.strip():
            bc_lengths = [args.bc_lengths.strip()]

        # Load barcode map
        bc_map = pd.read_csv(args.bc_map_file, sep="\t", header=None)
        print_log(f"Loaded barcode map with {len(bc_map)} entries")

        # Process barcodes using functional approach
        process_barcodes(
            bc_map=bc_map,
            bc_lengths=bc_lengths,
            output_files=output_files,
        )

    except Exception as e:
        print_error(f"Error processing barcodes: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
