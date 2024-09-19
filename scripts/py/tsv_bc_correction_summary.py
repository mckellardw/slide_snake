import argparse
from collections import defaultdict


def infer_barcode_count(line):
    fields = line.strip().split("\t")
    return (len(fields) - 1) // 4


def calculate_stats(filename, max_lines=None):
    stats = defaultdict(lambda: defaultdict(int))
    hamming_distances = defaultdict(lambda: defaultdict(list))
    second_match_distances = defaultdict(lambda: defaultdict(list))
    hamming_tally = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
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
                    hamming_distance = int(barcode_fields[3])
                    hamming_distances[barcode_num]["corrected"].append(hamming_distance)
                    hamming_tally[barcode_num]["corrected"][hamming_distance] += 1
                    if barcode_fields[4] != "N":
                        second_match_distances[barcode_num]["corrected"].append(
                            int(barcode_fields[3]) + int(barcode_fields[4])
                        )
                else:  # Uncorrected barcode
                    hamming_distance = int(barcode_fields[3])
                    hamming_distances[barcode_num]["uncorrected"].append(
                        hamming_distance
                    )
                    hamming_tally[barcode_num]["uncorrected"][hamming_distance] += 1
                    second_match_distances[barcode_num]["uncorrected"].append(
                        int(barcode_fields[3]) + int(barcode_fields[4])
                    )

    # Calculate means
    for barcode_num in stats:
        for correction_type in ["corrected", "uncorrected"]:
            if hamming_distances[barcode_num][correction_type]:
                stats[barcode_num][f"mean_hamming_{correction_type}"] = sum(
                    hamming_distances[barcode_num][correction_type]
                ) / len(hamming_distances[barcode_num][correction_type])
            else:
                stats[barcode_num][f"mean_hamming_{correction_type}"] = 0

            if second_match_distances[barcode_num][correction_type]:
                stats[barcode_num][f"mean_second_match_{correction_type}"] = sum(
                    second_match_distances[barcode_num][correction_type]
                ) / len(second_match_distances[barcode_num][correction_type])
            else:
                stats[barcode_num][f"mean_second_match_{correction_type}"] = 0

    return stats, hamming_tally, num_barcodes


def generate_stats_output(stats, hamming_tally, num_barcodes):
    output = []
    output.append(f"Number of barcodes per read: {num_barcodes}")
    for barcode_num in sorted(stats.keys()):
        output.append(f"\nStatistics for {barcode_num}:")
        output.append(f"Total reads processed: {stats[barcode_num]['total_reads']}")
        output.append(
            f"Number of corrected barcodes: {stats[barcode_num]['corrected_barcodes']}"
        )
        output.append(
            f"Mean hamming distance of corrected barcodes: {stats[barcode_num]['mean_hamming_corrected']:.2f}"
        )
        output.append(
            f"Mean hamming distance of uncorrected barcodes: {stats[barcode_num]['mean_hamming_uncorrected']:.2f}"
        )
        output.append(
            f"Mean second match hamming distance of corrected barcodes: {stats[barcode_num]['mean_second_match_corrected']:.2f}"
        )
        output.append(
            f"Mean second match hamming distance of uncorrected barcodes: {stats[barcode_num]['mean_second_match_uncorrected']:.2f}"
        )

        output.append("\nHamming distance tally for corrected barcodes:")
        for distance, count in sorted(hamming_tally[barcode_num]["corrected"].items()):
            output.append(f"  Distance {distance}: {count}")

        output.append("\nHamming distance tally for uncorrected barcodes:")
        for distance, count in sorted(
            hamming_tally[barcode_num]["uncorrected"].items()
        ):
            output.append(f"  Distance {distance}: {count}")

    return output


def write_stats_to_file(stats, hamming_tally, num_barcodes, output_file):
    output = generate_stats_output(stats, hamming_tally, num_barcodes)
    with open(output_file, "w") as f:
        for line in output:
            f.write(line + "\n")


def print_stats_to_stdout(stats, hamming_tally, num_barcodes):
    output = generate_stats_output(stats, hamming_tally, num_barcodes)
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

    stats, hamming_tally, num_barcodes = calculate_stats(
        args.input_file, args.max_lines
    )

    if args.output:
        write_stats_to_file(stats, hamming_tally, num_barcodes, args.output)
    else:
        print_stats_to_stdout(stats, hamming_tally, num_barcodes)


if __name__ == "__main__":
    main()
