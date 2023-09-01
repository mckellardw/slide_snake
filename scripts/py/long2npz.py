import numpy as np
import scipy.sparse as sp
import gzip
import argparse

# Usage:
## python long2npz.py input.tsv.gz output_matrix.npz

def save_sparse_matrix(matrix, var_names, row_names, output_path):
    sp.save_npz(output_path, matrix)
    np.savez(output_path.replace(".npz", "_metadata.npz"), var_names=var_names, row_names=row_names)
    
def main():
    parser = argparse.ArgumentParser(description="Convert a tsv.gz file to a sparse matrix in npz format.")
    parser.add_argument("tsv_gz_path", type=str, help="Path to the input tsv.gz file")
    parser.add_argument("output_path", type=str, help="Path to the output npz file")
    args = parser.parse_args()

    # Open the tsv.gz file
    with gzip.open(args.tsv_gz_path, 'rt') as f:
        # Skip the first line
        next(f)

        # Initialize lists to store data
        data_vars = []
        data_obs = []
        data = []

        # Process each line in the file
        for line in f:
            var, obs, val = line.strip().split('\t')
            # print(f"var: {var}, obs: {obs}, val: {val}")
            if float(val) != 0:  # Only add non-zero entries to the sparse matrix
                data_vars.append(var)
                data_obs.append(obs)
                data.append(float(val))

    # Get unique var_names
    var_names = list(set(data_vars))
    obs_names = list(set(data_obs))
    # print(f"vars - {len(var_names)}")
    # print(f"obs - {len(obs_names)}")

    # Create mappings between string names and integer indices
    var_name_to_idx = {var_name: idx for idx, var_name in enumerate(var_names)}
    row_name_to_idx = {row_name: idx for idx, row_name in enumerate(obs_names)}

    # Convert string names to integer indices
    data_vars = [var_name_to_idx[var_name] for var_name in data_vars]
    data_obs = [row_name_to_idx[row_name] for row_name in data_obs]

    # Initialize sparse matrix
    matrix = sp.csr_matrix(
        (data, (data_obs, data_vars)), 
        shape=(len(obs_names), len(var_names))
    )

    # Save sparse matrix to npz file
    save_sparse_matrix(
        matrix, 
        var_names,
        obs_names,
        args.output_path
    )
    print("Sparse matrix saved to", args.output_path)

if __name__ == "__main__":
    main()
