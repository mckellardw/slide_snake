import argparse
from collections import defaultdict

# Usage
"""
python scripts/py/tsv_bc_correction_summary.py {input.TSV_FULL} > {output.SUMMARY}
"""


def infer_barcode_count(line):
    fields = line.strip().split("\t")
    return (len(fields) - 1) // 4


def calculate_stats(filename, max_lines=None):
    stats = defaultdict(lambda: defaultdict(int))
    levenshtein_distances = defaultdict(lambda: defaultdict(list))
    second_match_distances = defaultdict(lambda: defaultdict(list))
    levenshtein_tally = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    num_barcodes = None

    with open(filename, "r") as file:
        for line_num, line in enumerate(file, 1):
            if max_lines and line_num > max_lines:
                break

            if line_num == 1:
                num_barcodes = infer_barcode_count(line)

            fields = line.strip().split("\t")

            for i in range(num_barcodes):
                barcode_fields = fields[i * 4 : (i + 1) * 4 + 1]
                barcode_num = f"barcode_{i+1}"

                stats[barcode_num]["total_reads"] += 1

                if barcode_fields[2] != "-":  # Corrected barcode
                    stats[barcode_num]["corrected_barcodes"] += 1
                    levenshtein_distance = int(barcode_fields[3])
                    levenshtein_distances[barcode_num]["corrected"].append(
                        levenshtein_distance
                    )
                    levenshtein_tally[barcode_num]["corrected"][
                        levenshtein_distance
                    ] += 1
                    if barcode_fields[4] != "N":
                        second_match_distances[barcode_num]["corrected"].append(
                            int(barcode_fields[3]) + int(barcode_fields[4])
                        )
                else:  # Uncorrected barcode
                    levenshtein_distance = int(barcode_fields[3])
                    levenshtein_distances[barcode_num]["uncorrected"].append(
                        levenshtein_distance
                    )
                    levenshtein_tally[barcode_num]["uncorrected"][
                        levenshtein_distance
                    ] += 1
                    second_match_distances[barcode_num]["uncorrected"].append(
                        int(barcode_fields[3]) + int(barcode_fields[4])
                    )

    # Calculate means
    for barcode_num in stats:
        for correction_type in ["corrected", "uncorrected"]:
            if levenshtein_distances[barcode_num][correction_type]:
                stats[barcode_num][f"mean_levenshtein_{correction_type}"] = sum(
                    levenshtein_distances[barcode_num][correction_type]
                ) / len(levenshtein_distances[barcode_num][correction_type])
            else:
                stats[barcode_num][f"mean_levenshtein_{correction_type}"] = 0

            if second_match_distances[barcode_num][correction_type]:
                stats[barcode_num][f"mean_second_match_{correction_type}"] = sum(
                    second_match_distances[barcode_num][correction_type]
                ) / len(second_match_distances[barcode_num][correction_type])
            else:
                stats[barcode_num][f"mean_second_match_{correction_type}"] = 0

    return stats, levenshtein_tally, num_barcodes


def generate_stats_output(stats, levenshtein_tally, num_barcodes):
    output = []
    output.append(f"Number of barcodes per read: {num_barcodes}")
    for barcode_num in sorted(stats.keys()):
        output.append(f"\nStatistics for {barcode_num}:")
        output.append(
            f"Total reads processed:        {stats[barcode_num]['total_reads']:,}"
        )
        output.append(
            f"Number of corrected barcodes: {stats[barcode_num]['corrected_barcodes']:,}"
        )
        output.append(
            f"Correction efficiency:        {stats[barcode_num]['corrected_barcodes']/stats[barcode_num]['total_reads']:.3f}"
        )
        output.append(f"")
        output.append(
            f"Mean levenshtein distance of corrected barcodes:   {stats[barcode_num]['mean_levenshtein_corrected']:.2f}"
        )
        output.append(
            f"Mean levenshtein distance of uncorrected barcodes: {stats[barcode_num]['mean_levenshtein_uncorrected']:.2f}"
        )
        output.append(f"")
        output.append(
            f"Mean 2nd match levenshtein distance of corrected barcodes:   {stats[barcode_num]['mean_second_match_corrected']:.2f}"
        )
        output.append(
            f"Mean 2nd match levenshtein distance of uncorrected barcodes: {stats[barcode_num]['mean_second_match_uncorrected']:.2f}"
        )

        output.append("\nLevenshtein distance tally for corrected barcodes:")
        total_corrected = sum(levenshtein_tally[barcode_num]["corrected"].values())
        for distance, count in sorted(
            levenshtein_tally[barcode_num]["corrected"].items()
        ):
            percentage = (count / total_corrected) * 100 if total_corrected > 0 else 0
            output.append(f"  Distance {distance}: {count:,} ({percentage:.2f}%)")

        output.append("\nLevenshtein distance tally for uncorrected barcodes:")
        total_uncorrected = sum(levenshtein_tally[barcode_num]["uncorrected"].values())
        for distance, count in sorted(
            levenshtein_tally[barcode_num]["uncorrected"].items()
        ):
            percentage = (
                (count / total_uncorrected) * 100 if total_uncorrected > 0 else 0
            )
            output.append(f"  Distance {distance}: {count:,} ({percentage:.2f}%)")

    return output


def write_stats_to_file(stats, levenshtein_tally, num_barcodes, output_file):
    output = generate_stats_output(stats, levenshtein_tally, num_barcodes)
    with open(output_file, "w") as f:
        for line in output:
            f.write(line + "\n")


def print_stats_to_stdout(stats, levenshtein_tally, num_barcodes):
    output = generate_stats_output(stats, levenshtein_tally, num_barcodes)
    for line in output:
        print(line)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate statistics for barcode correction results."
    )
    parser.add_argument(
        "input_file", help="Input TSV file containing barcode correction results"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file to write results (if not specified, prints to stdout)",
    )
    parser.add_argument(
        "-m",
        "--max-lines",
        type=int,
        help="Maximum number of lines to process (default: process all lines)",
    )

    args = parser.parse_args()

    stats, levenshtein_tally, num_barcodes = calculate_stats(
        args.input_file, args.max_lines
    )

    if args.output:
        write_stats_to_file(stats, levenshtein_tally, num_barcodes, args.output)
    else:
        print_stats_to_stdout(stats, levenshtein_tally, num_barcodes)


if __name__ == "__main__":
    main()
