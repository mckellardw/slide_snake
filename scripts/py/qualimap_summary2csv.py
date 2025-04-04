import pandas as pd
import argparse


def parse_qualimap_report(qualimapReport_txt):
    """Parse the Qualimap report text file into a dictionary."""
    with open(qualimapReport_txt, "r") as file:
        lines = file.readlines()

    out_dict = {}
    current_section = None

    for line in lines:
        line = line.strip()
        if not line:
            continue  # Skip empty lines
        if line.startswith(">>>>>>>"):
            current_section = line.strip(">>>>>>> ").strip()
            continue  # Skip section headers
        if "=" in line:
            key, value = line.split("=", 1)
            key = (
                f"{current_section} - {key.strip()}" if current_section else key.strip()
            )
            value = value.strip()
            if "(" in value:
                value, _ = value.split("(")
            if "%" in value:
                value = float(value.rstrip("%").strip().replace(",", ""))
            elif "," in value and value.replace(",", "").replace(".", "").isdigit():
                value = float(value.strip().replace(",", ""))
            else:
                value = value.strip()
            out_dict[key] = value

    return out_dict


def write_csv(data_dict, output_csv):
    """Write the parsed data dictionary to a CSV file. Replace spaces with underscores in column names."""
    out_df = pd.DataFrame([data_dict])
    out_df.columns = [col.replace(" - ", "_") for col in out_df.columns]
    out_df.columns = [col.replace("'", "") for col in out_df.columns]
    out_df.columns = [col.replace(" ", "_") for col in out_df.columns]
    # Remove commas from numeric values
    for col in out_df.columns:
        if out_df[col].dtype == object:
            out_df[col] = out_df[col].str.replace(",", "", regex=False)
    out_df.to_csv(output_csv, index=False)


def validate_file_types(input_file, output_file):
    """Ensure the input file is a .txt file and the output file is a .csv file."""
    if not input_file.endswith(".txt"):
        raise ValueError(f"Input file must be a .txt file. Provided: {input_file}")
    if not output_file.endswith(".csv"):
        raise ValueError(f"Output file must be a .csv file. Provided: {output_file}")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert Qualimap report text file to CSV format."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the Qualimap report text file.",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to the output CSV file."
    )
    return parser.parse_args()


def main():
    """Main function to parse Qualimap report and write to CSV."""
    args = parse_args()
    validate_file_types(args.input, args.output)
    data_dict = parse_qualimap_report(args.input)
    write_csv(data_dict, args.output)


if __name__ == "__main__":
    main()
