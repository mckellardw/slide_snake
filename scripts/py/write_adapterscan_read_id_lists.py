#!/usr/bin/env python

import argparse
import csv
import os

def process_tsv_file(tsv_file_path, output_directory):
    """
    Reads a TSV file line by line, finds all unique values in the "lab" column,   
    and writes a text file for each unique value containing the corresponding   
    "read_id" values.

    Parameters
    ----------
    tsv_file_path : str
        The path to the TSV file.
    output_directory : str
        The directory where the text files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Prepare a dictionary to hold the mapping of "lab" values to "read_id" values
    lab_to_read_ids = {}

    # Open the TSV file and read it line by line
    with open(tsv_file_path, 'r') as tsv_file:
        tsv_reader = csv.DictReader(tsv_file, delimiter='\t')
        for row in tsv_reader:
            # Extract the "lab" and "read_id" values from the current row
            lab = row['lab']
            read_id = row['read_id']

            # If the "lab" value is not yet in the dictionary, add it
            if lab not in lab_to_read_ids:
                lab_to_read_ids[lab] = []

            # Append the "read_id" to the list of "read_id" values for this "lab"
            lab_to_read_ids[lab].append(read_id)

    # Write the "read_id" values to a text file for each "lab" value
    for lab, read_ids in lab_to_read_ids.items():
        with open(os.path.join(output_directory, f'{lab}.txt'), 'w') as output_file:
            for read_id in read_ids:
                output_file.write(f'{read_id}\n')

    print("Text files have been written successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a TSV file to create text files for each unique "lab" value.')
    parser.add_argument('--tsv_file_path', type=str, help='The path to the TSV file.')
    parser.add_argument('--output_directory', type=str, help='The directory where the text files will be saved.')
    args = parser.parse_args()

    process_tsv_file(args.tsv_file_path, args.output_directory)
