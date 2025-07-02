import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import gzip
import argparse

# Usage:
## python long2mtx.py input.tsv.gz output_matrix.npz


def save_npz_matrix(matrix, var_names, obs_names, output_dir):
    sp.save_npz(output_dir, matrix)
    np.savez(
        output_dir.replace(".npz", "_metadata.npz"),
        var_names=var_names,
        obs_names=obs_names,
    )


# Save mtx file & features/barcodes, in the naming convention of STARsolo
def save_sparse_matrix(matrix, var_names, obs_names, output_dir):
    # Save the sparse matrix in Matrix Market (.mtx) format
    mmwrite(f"{output_dir}/matrix.mtx", matrix, field="real")

    # Save var_names (features) and obs_names (barcodes) as text files
    with open(f"{output_dir}/features.tsv", "w") as var_file:
        var_file.write("\n".join(var_names))

    with open(f"{output_dir}/barcodes.tsv", "w") as obs_file:
        obs_file.write("\n".join(obs_names))

    # Gzip the .mtx file and associated .txt files
    with open(f"{output_dir}/matrix.mtx", "rb") as file_in, gzip.open(
        f"{output_dir}/matrix.mtx" + ".gz", "wb"
    ) as file_out:
        file_out.writelines(file_in)

    with open(f"{output_dir}/features.tsv", "rb") as file_in, gzip.open(
        f"{output_dir}/features.tsv.gz", "wb"
    ) as file_out:
        file_out.writelines(file_in)

    with open(f"{output_dir}/barcodes.tsv", "rb") as file_in, gzip.open(
        f"{output_dir}/barcodes.tsv.gz", "wb"
    ) as file_out:
        file_out.writelines(file_in)

    # Remove the original uncompressed files
    import os

    os.remove(f"{output_dir}/matrix.mtx")
    os.remove(f"{output_dir}/features.tsv")
    os.remove(f"{output_dir}/barcodes.tsv")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert a tsv.gz file to a sparse matrix in mtx format."
    )
    parser.add_argument("tsv_gz_path", type=str, help="Path to the input tsv.gz file")
    parser.add_argument("output_dir", type=str, help="Path to the output directory")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose output (default: False)",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.verbose:
        print(f"Processing file: {args.tsv_gz_path}")
        print(f"Output directory: {args.output_dir}")

    # Open the tsv.gz file
    with gzip.open(args.tsv_gz_path, "rt") as f:
        # Skip the first line
        next(f)

        if args.verbose:
            print("Reading data from tsv.gz file...")

        # Initialize lists to store data
        data_vars = []
        data_obs = []
        data = []

        # Process each line in the file
        for line in f:
            var, obs, val = line.strip().split("\t")
            # print(f"var: {var}, obs: {obs}, val: {val}")
            if float(val) != 0:  # Only add non-zero entries to the sparse matrix
                data_vars.append(var)
                data_obs.append(obs)
                data.append(float(val))

    # Get unique var_names
    var_names = list(set(data_vars))
    obs_names = list(set(data_obs))

    if args.verbose:
        print(f"Found {len(var_names)} unique features")
        print(f"Found {len(obs_names)} unique observations")
        print(f"Non-zero entries: {len(data)}")
    # print(f"vars - {len(var_names)}")
    # print(f"obs - {len(obs_names)}")

    # Create mappings between string names and integer indices
    var_name_to_idx = {var_name: idx for idx, var_name in enumerate(var_names)}
    row_name_to_idx = {row_name: idx for idx, row_name in enumerate(obs_names)}

    # Convert string names to integer indices
    data_vars = [var_name_to_idx[var_name] for var_name in data_vars]
    data_obs = [row_name_to_idx[row_name] for row_name in data_obs]

    if args.verbose:
        print("Creating sparse matrix...")

    # Initialize sparse matrix
    matrix = sp.csr_matrix(
        (data, (data_obs, data_vars)), shape=(len(obs_names), len(var_names))
    )

    if args.verbose:
        print(f"Matrix shape: {matrix.shape}")
        print("Saving matrix to output directory...")

    # Save sparse matrix to npz file
    # save_npz_matrix(
    #     matrix,
    #     var_names,
    #     obs_names,
    #     args.output_dir
    # )

    save_sparse_matrix(matrix, var_names, obs_names, args.output_dir)

    print("Sparse matrix saved to", args.output_dir)


if __name__ == "__main__":
    main()
