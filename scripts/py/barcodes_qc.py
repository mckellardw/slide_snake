import argparse
import pysam
import csv

def hamming(k1, k2):
    """Calculate the Hamming distance between two strings."""
    if len(k1) != len(k2):
        raise "kmers not the same length"
    first = np.array(list(k1))
    second = np.array(list(k2))
    dist = (first != second).sum()
    return dist

def parse_bam_file(bam_file, raw_barcode_tag, corrected_barcode_tag, output_file):
    """Parse the BAM file and compute Hamming distances."""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    no_corr_barcode_count = 0
    total_reads = 0
    results = []

    for read in samfile:
        total_reads += 1
        raw_barcode = read.get_tag(raw_barcode_tag)
        corrected_barcode = read.get_tag(corrected_barcode_tag)

        if corrected_barcode:
            distance = hamming(raw_barcode, corrected_barcode)
            results.append([raw_barcode, corrected_barcode, distance])
        else:
            no_corr_barcode_count += 1

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['Raw Barcode', 'Corrected Barcode', 'Hamming Distance'])
        writer.writerows(results)

    print(f"Total reads: {total_reads}, Reads without corrected barcode: {no_corr_barcode_count}")

def main():
    parser = argparse.ArgumentParser(description='Parse BAM file and compute Hamming distances between raw and corrected cell barcodes.')
    parser.add_argument('bam_file', help='Path to the BAM file.')
    parser.add_argument('raw_barcode_tag', help='Tag for raw cell barcodes.')
    parser.add_argument('corrected_barcode_tag', help='Tag for corrected cell barcodes.')
    parser.add_argument('output_file', help='Path to the output CSV or TSV file.')
    args = parser.parse_args()

    parse_bam_file(args.bam_file, args.raw_barcode_tag, args.corrected_barcode_tag, args.output_file)

if __name__ == "__main__":
    main()
