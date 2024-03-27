import pandas as pd


def qualimap_summary2csv_STAR(qualimapReport_txt, qualimapReport_csv):
    with open(qualimapReport_txt, "r") as file:
        lines = file.readlines()

    out_dict = {}
    for line in lines:
        if "=" in line:
            key, value = line.split("=")
            if "(" in value:
                value, _ = value.split("(")
                value = float(value.strip().replace(",", ""))
            elif "%" in value:
                value = float(value.rstrip("%").strip().replace(",", ""))
            elif "," in value:
                value = float(value.strip().replace(",", ""))

            out_dict[key.strip()] = value

    out_df = pd.DataFrame.from_dict(out_dict, orient="index")
    out_df.T.to_csv(qualimapReport_csv, index=False)


if __name__ == "__main__":
    import sys

    qualimapReport_txt = sys.argv[1]
    qualimapReport_csv = sys.argv[2]
    qualimap_summary2csv_STAR(qualimapReport_txt, qualimapReport_csv)
