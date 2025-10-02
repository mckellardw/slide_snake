import os
import re
import glob
import pandas as pd


#### Utility rules ############################
# Index .bam file
rule utils_index_BAM:
    input:
        BAM="{BAM}",
    output:
        BAI="{BAM}.bai",
    # wildcard_constraints:
    #     BAM=".*\.(bam)$"
    # resources:
    threads: config["CORES"]
    conda:
        f"{workflow.basedir}/envs/sambamba.yml"
    shell:
        """
        sambamba index -q -t {threads} {input.BAM}
        """


#### Util functions ###########################
def check_sample_sheet(SAMPLE_SHEET):
    short_read_requirements = {
        "fastq_R1": [".fq.gz", ".fastq.gz"],
        "fastq_R2": [".fq.gz", ".fastq.gz"],
        "bwa_rRNA_ref": [".fa.gz"],
        "rRNA_gtf": [".gtf"],
        "STAR_ref": [""],
        "kb_idx": [".idx"],
        "kb_t2g": [".txt"],
        "genes_gtf": [".gtf"],
        "cdna_fa": [".fa", ".fa.gz"],  # Added for rRNA filtering
    }

    ONT_requirements = {
        "ONT": [".fq.gz", ".fastq.gz", ".bam", ".cram"],
        "genome_fa": [".fa"],
        "genes_gtf": [".gtf"],
        "cdna_fa": [".fa", ".fa.gz"],
    }

    # Recipe-specific requirements
    recipe_specific_requirements = {
        "rRNA-bwa": {
            "short_read": ["cdna_fa"],  # Need cDNA FASTA for BWA rRNA filtering
            "ont": [],
        },
        "ribodetector": {"short_read": [], "ont": []},
        "visium": {
            "short_read": ["STAR_ref", "genes_gtf"],
            "ont": ["genome_fa", "genes_gtf", "cdna_fa"],
        },
        "seeker": {
            "short_read": ["STAR_ref", "genes_gtf", "kb_idx", "kb_t2g"],
            "ont": ["genome_fa", "genes_gtf", "cdna_fa"],
        },
        "stomics": {
            "short_read": ["STAR_ref", "genes_gtf"],
            "ont": ["genome_fa", "genes_gtf", "cdna_fa"],
        },
        "decoder": {
            "short_read": ["STAR_ref", "genes_gtf"],
            "ont": ["genome_fa", "genes_gtf", "cdna_fa"],
        },
        "miST": {
            "short_read": ["STAR_ref", "genes_gtf"],
            "ont": ["genome_fa", "genes_gtf", "cdna_fa"],
        },
    }

    missing_files = []
    empty_fields = []
    recipe_errors = []
    duplicate_samples = SAMPLE_SHEET["sampleID"].duplicated()
    duplicate_sample_names = SAMPLE_SHEET["sampleID"][duplicate_samples].tolist()

    for index, row in SAMPLE_SHEET.iterrows():
        sample_id = row["sampleID"]

        # Check short read recipes and requirements
        if row["recipe"]:
            recipes = row["recipe"].split()

            # Check recipe-specific requirements
            for recipe in recipes:
                if recipe in recipe_specific_requirements:
                    required_fields = recipe_specific_requirements[recipe]["short_read"]
                    for field in required_fields:
                        if field not in row or not row[field]:
                            recipe_errors.append(
                                (
                                    sample_id,
                                    f"Recipe '{recipe}' requires field '{field}' but it's missing or empty",
                                )
                            )
                        elif field in short_read_requirements:
                            # Check file existence and extension
                            file_paths = str(row[field]).split()
                            for file_path in file_paths:
                                if file_path and file_path != "nan":
                                    if not any(
                                        file_path.endswith(ext)
                                        for ext in short_read_requirements[field]
                                    ):
                                        recipe_errors.append(
                                            (
                                                sample_id,
                                                f"File {file_path} for recipe '{recipe}' field '{field}' does not have correct extension {short_read_requirements[field]}",
                                            )
                                        )
                                    elif not glob.glob(file_path) and not re.search(
                                        r"\{.*\}", file_path
                                    ):
                                        missing_files.append((sample_id, file_path))

            # Check general short read requirements
            for column, exts in short_read_requirements.items():
                if column not in row:
                    continue
                file_paths = str(row[column]).split() if row[column] else []
                if not file_paths or (len(file_paths) == 1 and file_paths[0] == ""):
                    # Only mark as empty if it's required for the recipes being used
                    is_required = any(
                        column
                        in recipe_specific_requirements.get(recipe, {}).get(
                            "short_read", []
                        )
                        for recipe in recipes
                    )
                    if is_required:
                        empty_fields.append((sample_id, column))
                else:
                    for file_path in file_paths:
                        if file_path and file_path != "nan":
                            if column == "STAR_ref" and not os.path.isdir(file_path):
                                print(
                                    f"WARNING: {file_path} for sample {sample_id} was not found, or is not a directory."
                                )
                                missing_files.append((sample_id, file_path))
                            elif not any(
                                re.search(f"{ext}$", file_path) for ext in exts
                            ):
                                print(
                                    f"WARNING: File {file_path} for sample {sample_id} does not have the correct extension {exts}"
                                )
                            elif not glob.glob(file_path) and not re.search(
                                r"\{.*\}", file_path
                            ):
                                missing_files.append((sample_id, file_path))

        # Check ONT recipes and requirements
        if row["recipe_ONT"]:
            ont_recipes = row["recipe_ONT"].split()

            # Check recipe-specific requirements for ONT
            for recipe in ont_recipes:
                if recipe in recipe_specific_requirements:
                    required_fields = recipe_specific_requirements[recipe]["ont"]
                    for field in required_fields:
                        if field not in row or not row[field]:
                            recipe_errors.append(
                                (
                                    sample_id,
                                    f"ONT recipe '{recipe}' requires field '{field}' but it's missing or empty",
                                )
                            )
                        elif field in ONT_requirements:
                            # Check file existence and extension
                            file_paths = str(row[field]).split()
                            for file_path in file_paths:
                                if file_path and file_path != "nan":
                                    if not any(
                                        file_path.endswith(ext)
                                        for ext in ONT_requirements[field]
                                    ):
                                        recipe_errors.append(
                                            (
                                                sample_id,
                                                f"File {file_path} for ONT recipe '{recipe}' field '{field}' does not have correct extension {ONT_requirements[field]}",
                                            )
                                        )
                                    elif not glob.glob(file_path) and not re.search(
                                        r"\{.*\}", file_path
                                    ):
                                        missing_files.append((sample_id, file_path))

            # Check general ONT requirements
            for column, exts in ONT_requirements.items():
                if column not in row:
                    continue
                file_paths = str(row[column]).split() if row[column] else []
                if not file_paths or (len(file_paths) == 1 and file_paths[0] == ""):
                    # Only mark as empty if it's required for the recipes being used
                    is_required = any(
                        column
                        in recipe_specific_requirements.get(recipe, {}).get("ont", [])
                        for recipe in ont_recipes
                    )
                    if is_required:
                        empty_fields.append((sample_id, column))
                else:
                    for file_path in file_paths:
                        if (
                            file_path
                            and file_path != "nan"
                            and not any(re.search(f"{ext}$", file_path) for ext in exts)
                        ):
                            print(
                                f"WARNING: File {file_path} for sample {sample_id} does not have the correct extension {exts}"
                            )

    # Report issues
    if empty_fields:
        print("WARNING: The following required fields are empty:")
        for sample_id, column in empty_fields:
            print(f"Sample {sample_id}: {column}")

    if recipe_errors:
        print("WARNING: Recipe-specific requirement errors:")
        for sample_id, error in recipe_errors:
            print(f"Sample {sample_id}: {error}")

    if duplicate_samples.any():
        print("WARNING: Duplicate sample IDs found in the sample sheet:")
        for sample_id in duplicate_sample_names:
            print(f"Sample {sample_id}")

    if missing_files:
        print("WARNING: The following files are missing:")
        for sample_id, file_path in missing_files:
            print(f"Sample {sample_id}: {file_path}")
        raise FileNotFoundError(
            "Some required files are missing. Please check the log for details."
        )


def check_recipe_sheet(RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT):
    required_columns = [
        "R1_finalLength",
        "fwd_primer",
        "rev_primer",
        "BC_adapter",
        "BC_length",
        "BC_offset",
        "BC_position",
        "BC_max_ED",
        "BC_min_ED_diff",
        "BC_concat",
        "UMI_adapter",
        "UMI_length",
        "UMI_offset",
        "UMI_position",
        "internal_adapter",
        "STAR_soloType",
        "STAR_soloCBmatchWLtype",
        "STAR_soloCB",
        "STAR_soloUMI",
        "STAR_soloAdapter",
        "STAR_extra",
        "kb_x",
        "kb_extra",
        "featureCounts_extra",
        "mm2_extra",
    ]

    # Recipe-specific required columns
    recipe_specific_columns = {
        "rRNA-bwa": ["fwd_primer", "rev_primer"],
        "ribodetector": ["fwd_primer", "rev_primer"],
        "visium": [
            "STAR_soloType",
            "STAR_soloCB",
            "STAR_soloUMI",
            "BC_length",
            "UMI_length",
        ],
        "seeker": ["BC_length", "UMI_length", "STAR_soloType", "kb_x"],
        "stomics": [
            "STAR_soloType",
            "STAR_soloCB",
            "STAR_soloUMI",
            "BC_length",
            "UMI_length",
        ],
        "decoder": ["BC_length", "UMI_length", "STAR_soloType"],
        "miST": ["BC_length", "UMI_length", "STAR_soloType"],
    }

    missing_columns = [
        col for col in required_columns if col not in RECIPE_SHEET.columns
    ]
    duplicate_recipes = RECIPE_SHEET.index.duplicated()
    duplicate_recipe_names = RECIPE_SHEET.index[duplicate_recipes].tolist()
    recipe_config_errors = []

    if missing_columns:
        print(
            "WARNING: The following required columns are missing from the recipe sheet:"
        )
        for col in missing_columns:
            print(f"- {col}")

    if duplicate_recipes.any():
        print("WARNING: Duplicate recipes found in the recipe sheet:")
        for recipe in duplicate_recipe_names:
            print(f"Recipe {recipe}")

    # Get all recipes used in sample sheets
    RECIPES = unlist(RECIPE_DICT, unique=True) + unlist(RECIPE_ONT_DICT, unique=True)

    for RECIPE in RECIPES:
        if RECIPE not in list(RECIPE_SHEET.index):
            print(f"ERROR: Recipe '{RECIPE}' not found in recipe sheet!")
            print(f"Fix RECIPE_SHEET [{config['RECIPE_SHEET']}] and try again!")
            recipe_config_errors.append(f"Recipe '{RECIPE}' missing from recipe sheet")
        else:
            # Check recipe-specific requirements
            if RECIPE in recipe_specific_columns:
                required_cols = recipe_specific_columns[RECIPE]
                for col in required_cols:
                    if col in RECIPE_SHEET.columns:
                        value = RECIPE_SHEET.loc[RECIPE, col]
                        if pd.isna(value) or str(value).strip() == "":
                            recipe_config_errors.append(
                                f"Recipe '{RECIPE}' missing required parameter '{col}'"
                            )
                    else:
                        recipe_config_errors.append(
                            f"Recipe '{RECIPE}' requires column '{col}' which is missing from recipe sheet"
                        )

            # Validate specific parameter formats
            try:
                # Check BC_length format
                if "BC_length" in RECIPE_SHEET.columns and not pd.isna(
                    RECIPE_SHEET.loc[RECIPE, "BC_length"]
                ):
                    bc_length = str(RECIPE_SHEET.loc[RECIPE, "BC_length"])
                    if bc_length.strip():
                        # Should be either a single number or space-separated numbers
                        try:
                            [int(x) for x in bc_length.split()]
                        except ValueError:
                            recipe_config_errors.append(
                                f"Recipe '{RECIPE}' has invalid BC_length format: '{bc_length}'"
                            )

                # Check UMI_length format
                if "UMI_length" in RECIPE_SHEET.columns and not pd.isna(
                    RECIPE_SHEET.loc[RECIPE, "UMI_length"]
                ):
                    umi_length = str(RECIPE_SHEET.loc[RECIPE, "UMI_length"])
                    if umi_length.strip():
                        try:
                            int(umi_length)
                        except ValueError:
                            recipe_config_errors.append(
                                f"Recipe '{RECIPE}' has invalid UMI_length format: '{umi_length}'"
                            )

                # Check BC_concat format
                if "BC_concat" in RECIPE_SHEET.columns and not pd.isna(
                    RECIPE_SHEET.loc[RECIPE, "BC_concat"]
                ):
                    bc_concat = str(RECIPE_SHEET.loc[RECIPE, "BC_concat"]).lower()
                    if bc_concat not in ["true", "false", "1", "0"]:
                        recipe_config_errors.append(
                            f"Recipe '{RECIPE}' has invalid BC_concat value: '{bc_concat}' (should be True/False)"
                        )

            except KeyError as e:
                recipe_config_errors.append(
                    f"Recipe '{RECIPE}' not found in recipe sheet: {e}"
                )

    if recipe_config_errors:
        print("ERROR: Recipe configuration errors found:")
        for error in recipe_config_errors:
            print(f"- {error}")
        raise ValueError(
            "Recipe configuration errors found. Please fix your recipe sheet."
        )


# Select input reads based on alignment recipe
def get_fqs(w, return_type=["list", "dict"], mode=["ONT", "ILMN"]):
    # param defaults
    if type(return_type) is list and len(return_type) > 1:  # default
        return_type = "list"
    if type(mode) is list:  # default
        mode = "ONT"

    # get file paths
    try:
        if mode == "ILMN":
            if "rRNA-bwa" in w.RECIPE:  # Use trimmed & bwa-rRNA-filtered .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/bwa/no_rRNA_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/bwa/no_rRNA_R2.fq.gz"

                # additional R1 trimming options to remove bridge adapter
                if "internalTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/bwa/no_rRNA_internalTrim_R1.fq.gz"
                if "hardTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/bwa/no_rRNA_hardTrim_R1.fq.gz"

            elif "ribodetector" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_R2.fq.gz"

                # additional R1 trimming options to remove bridge adapter
                if "internalTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_internalTrim_R1.fq.gz"
                if "hardTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_hardTrim_R1.fq.gz"
            else:  # just trimmed .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz"

                if "internalTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_internalTrim_R1.fq.gz"
                if "hardTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_hardTrim_R1.fq.gz"
        elif mode == "ONT":
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_R2.fq.gz"
            if "internalTrim" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_internalTrim_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/twiceCut_internalTrim_R2.fq.gz"
            elif "hardTrim" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_hardTrim_R1.fq.gz"
            else:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_R1.fq.gz"
        else:
            print("get_fqs(): `mode` not found")
    except Exception:
        if mode == "ILMN":
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/tmp/twiceCut_R2.fq.gz"
        elif mode == "ONT":
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/ont/tmp/cut_R2.fq.gz"
        else:
            print("get_fqs(): `mode` not found")

    # return fq path(s)
    if return_type == "list":
        return [R1, R2]
    elif return_type == "dict":
        return {"R1": R1, "R2": R2}


def get_fwd_primer(w, mode=["ONT", "ILMN"]):
    """
    Get the forward primer sequence from the recipe sheet.
    w: wildcards
    mode: which workflow this is for
    """
    try:
        if mode == "ILMN":
            fwd_primer = list(set(get_recipe_info(w, "fwd_primer")))
        elif mode == "ONT":
            fwd_primer = list(set(get_recipe_info(w, "fwd_primer")))
        else:
            print("get_fwd_primer(): `mode` not found")
            return "ERROR"
    except Exception:
        if mode == "ILMN":
            print("TODO")
        elif mode == "ONT":
            print("TODO")
        else:
            print("get_fwd_primer(): `mode` not found")

    return fwd_primer


def get_rev_primer(w, mode=["ONT", "ILMN"]):
    """
    Get the reverse primer sequence from the recipe sheet.
    w: wildcards
    mode: which workflow this is for
    """
    try:
        if mode == "ILMN":
            rev_primer = list(set(get_recipe_info(w, "rev_primer")))
        elif mode == "ONT":
            rev_primer = list(set(get_recipe_info(w, "rev_primer")))
        else:
            print("get_rev_primer(): `mode` not found")
            return "ERROR"
    except Exception:
        if mode == "ILMN":
            print("TODO")
        elif mode == "ONT":
            print("TODO")
        else:
            print("get_rev_primer(): `mode` not found")

    return rev_primer


def get_bc_adapter(w, mode=["ONT", "ILMN"]):
    try:
        if mode == "ILMN":
            if "Trim" in w.RECIPE or "std" in w.RECIPE:
                return get_recipe_info(
                    w, "BC_start"
                )  # return integer start positions for each BC sequence
            elif "matchLinker" in w.RECIPE or "seeker" in w.RECIPE:
                return get_recipe_info(w, "BC_adapter")
        elif mode == "ONT":
            return get_recipe_info(w, "BC_adapter")
        else:
            print("get_bc_adapter(): `mode` not found")
            return "ERROR"
    except Exception:
        if mode == "ILMN":
            print("TODO")
        elif mode == "ONT":
            print("TODO")
        else:
            print("get_bc_adapter(): `mode` not found")


def get_umi_adapter(w, mode=["ONT", "ILMN"]):
    try:
        if mode == "ILMN":
            if "Trim" in w.RECIPE or "std" in w.RECIPE:
                return get_recipe_info(
                    w, "UMI_start"
                )  # return integer start positions for each BC sequence
            elif "matchLinker" in w.RECIPE or "seeker" in w.RECIPE:
                return get_recipe_info(w, "UMI_adapter")
        elif mode == "ONT":
            return get_recipe_info(w, "UMI_adapter")
        else:
            print("get_umi_adapter(): `mode` not found")
            return "ERROR"
    except Exception:
        if mode == "ILMN":
            print("TODO")
        elif mode == "ONT":
            print("TODO")
        else:
            print("get_umi_adapter(): `mode` not found")


# whitelist param handling for different recipes/technologies/chemistries/etc
def get_whitelist(w, return_type=None, mode="STAR"):
    try:
        if "internalTrim" in w.RECIPE:
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"
        elif "adapterInsert" in w.RECIPE:
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_adapter.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_adapter.txt"
        elif "hardTrim" in w.RECIPE:
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"
        elif "matchLinker" in w.RECIPE or "seeker_std" in w.RECIPE:
            # Use BC_concat=False for combinatorial barcode constructs (DBIT, microST, etc)
            if return_type == "list":
                if RECIPE_SHEET["BC_concat"][w.RECIPE] and mode in ["STAR", "ILMN"]:
                    # Barcode constructs where positional barcodes are NOT independent (must be concatenated)
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt",
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt",
                    ]
                elif RECIPE_SHEET["BC_concat"][w.RECIPE] and mode == "ONT":
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt",
                    ]
                else:
                    # Barcode constructs where positional barcodes ARE independent
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_uniq_1.txt",
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_uniq_2.txt",
                    ]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt {w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt"
        else:
            # default
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"
    except Exception:
        if return_type == "list":
            whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
        else:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"

    # return whitelist path(s)
    return whitelist


def get_default_whitelist(w, return_type=None):
    try:
        if any("seeker" in recipe for recipe in get_recipes(w, mode="ILMN")):
            if return_type == "list":
                whitelist = [
                    f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt",
                    f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt",
                ]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt {w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt"
        else:
            # visium, stomics, microST
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"
    except Exception:
        if return_type == "list":
            whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
        else:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"

    # return whitelist path(s)
    return whitelist


def get_all_used_whitelists(w, return_type=None):
    try:
        recipes = get_recipes(w.SAMPLE)
        whitelist = []
        for recipe in recipes:
            whitelist.append(get_whitelist(w, mode="recipe"))
    except Exception:
        if return_type == "list":
            whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
        else:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"

    # return whitelist path(s)
    return whitelist


# Barcode map param handling for different recipes/technologies/chemistries/etc
def get_bc_map(w, mode=["ONT", "ILMN"]):
    try:
        if "matchLinker" in w.RECIPE:
            if mode == "ILMN":
                bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map_underscore.txt"
            elif mode == "ONT":
                bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map.txt"
        elif "internalTrim" in w.RECIPE:
            bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map.txt"
        elif "adapterInsert" in w.RECIPE:
            bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map_adapter.txt"
        else:
            # visium, stomics, microST
            bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map.txt"
    except Exception:
        bc_map = f"{w.OUTDIR}/{w.SAMPLE}/bc/map.txt"

    # return whitelist path(s)
    return bc_map


# Get number of barcodes (does not include UMIs!)
def get_n_cells(w):
    """
    Get the number of cells/spots/beads/whatever

    w: wildcards
    """
    whitelist_path = get_whitelist(w, return_type="list")[0]
    n_bcs = count_lines_in_file(whitelist_path)

    return n_bcs


# Get number of barcodes (does not include UMIs!)
def get_n_bcs(w):
    """
    Get the number of barcode positions

    w: wildcards
    """
    bc_lengths = get_recipe_info(w, info_col="BC_length").split()
    n_bcs = len(bc_lengths)
    return n_bcs


def get_bwa_ref(w, mode=["genome", "rRNA"]):
    """
    Pull BWA info from sample sheet

    w: wildcards
    mode: which reference to pull
    """
    try:
        if type(mode) is list and len(mode) > 1:  # default
            mode = "rRNA"

        if mode == "genome":
            out = "No genome reference supported for BWA!"
        elif mode == "rRNA":
            out = SAMPLE_SHEET["bwa_rRNA_ref"][w.SAMPLE]
    except Exception:
        out = "No reference given! Check your sample sheet!"

    return out


def get_STAR_ref(w, mode="genome"):
    """
    Pull STAR info from sample sheet

    w: wildcards
    mode: which reference to pull
    """
    try:
        if type(mode) is list and len(mode) > 1:  # default
            mode = "genome"

        if mode == "genome":
            out = SAMPLE_SHEET["STAR_ref"][w.SAMPLE]
        else:
            print(f"Don't know what to do for mode [{mode}]")
    except Exception:
        out = "No reference given! Check your sample sheet!"

    return out


def get_kallisto_ref(w, mode=["idx", "t2g", "idx_velo", "t2g_velo"]):
    """
    Pull kalllisto info from sample sheet

    w: wildcards
    mode: which reference to pull
    """
    try:
        if mode == "idx":
            out = SAMPLE_SHEET["kb_idx"][w.SAMPLE]
        elif mode == "t2g":
            out = SAMPLE_SHEET["kb_t2g"][w.SAMPLE]
    except Exception:
        out = "No reference given! Check your sample sheet!"

    return out


# TODO- add STAR_rRNA functinoality (no w.RECIPE accessibility)
def get_STAR_extra_params(w):
    star_info = {
        "STAR_soloType": "",
        "STAR_soloUMI": "",
        "STAR_soloCB": "",
        "STAR_soloAdapter": "",
        "STAR_extra": "--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.05  --outFilterMatchNmin 12  --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0",
    }

    # Iterate over each key
    for key in star_info.keys():
        try:
            star_info[key] = RECIPE_SHEET[key][w.RECIPE]
        except Exception:
            # values for default
            recipe = get_recipes(w, mode=f"concatenate")

            if "stomics" in recipe:
                star_info[key] = RECIPE_SHEET[key]["stomics_total"]
            elif "visium" in recipe:
                star_info[key] = RECIPE_SHEET[key]["visium_total"]
            elif "seeker" in recipe and "matchLinker" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_matchLinker_total"]
            elif "seeker" in recipe and "internalTrim" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_internalTrim_total"]
            elif "seeker" in recipe and "adapterInsert" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_adapterInsert_total"]
            elif "seeker" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_total"]
            elif "miST" in recipe:
                star_info[key] = RECIPE_SHEET[key]["miST_ligation_total"]
            elif "decoder" in recipe:
                star_info[key] = RECIPE_SHEET[key]["decoderseq_total"]
            else:
                # star_info[key] = RECIPE_SHEET[key]["visium_total"]
                star_info[key] = "Couldn't find the recipe!"

    # remove empty values
    # star_info = [k for k, v in star_info.items() if v]
    return star_info


def get_recipes(w, mode=["ONT", "ILMN", "list"]):
    """
    Retrieves the recipe for a given sample based on the mode specified.

    Parameters:
    - w (object): A wildcards object containing sample information, including a 'SAMPLE' attribute.
    - mode (list, optional): A list specifying the mode of operation. Defaults to ["ONT","ILMN","list"].

    Returns:
    - str or list: The recipe for the sample if found, otherwise an empty string or list.

    Raises:
    - KeyError: If the sample is not found in the specified recipe dictionary.

    Note:
    - The function first checks if the mode contains "ONT" or "ILMN" to determine the recipe dictionary to use.
    - If "list" is in the mode, it returns the recipe as a list; otherwise, it returns the recipe as a space-separated string.
    - If no recipe is found, it prints a warning message and returns an empty string.
    """
    # default option
    if isinstance(mode, list) and len(mode) > 1:
        mode = "ONT"

    if "ONT" in mode:
        try:
            return RECIPE_ONT_DICT[w.SAMPLE]
        except Exception:
            print(
                f"No ONT recipe given for sample '{w.SAMPLE}'! Check your sample sheet!"
            )
            return ""
    elif "ILMN" in mode:
        try:
            return RECIPE_DICT[w.SAMPLE]
        except Exception:
            print(
                f"No ILMN recipe given for sample '{w.SAMPLE}'! Check your sample sheet!"
            )
            return ""
    else:
        all_recipes = []
        if "RECIPE_ONT_DICT" in globals():
            if w.SAMPLE in RECIPE_ONT_DICT:
                all_recipes.extend(RECIPE_ONT_DICT[w.SAMPLE])
        if "RECIPE_DICT" in globals():
            if w.SAMPLE in RECIPE_DICT:
                all_recipes.extend(RECIPE_DICT[w.SAMPLE])
        all_recipes = unlist(all_recipes)
        if len(all_recipes) > 0:
            if "list" in mode:
                return all_recipes
            else:
                return " ".join(all_recipes)
        else:
            print(f"No recipe found for {w.SAMPLE}! Check your sample sheet!")
            return ""


# Pull info from recipe sheet
def get_recipe_info(w, info_col, mode=["ONT", "ILMN", "list"]):
    try:
        return RECIPE_SHEET[info_col][w.RECIPE]
    except:
        # print(f"Couldn't find `{info_col}`")
        recipe = unlist(get_recipes(w, mode=mode))

        return [RECIPE_SHEET[info_col][x] for x in recipe]


def count_lines_in_file(file_path):
    with open(file_path, "r") as file:
        return sum(1 for _ in file)


# Convert a megabyte value (str) to bytes [int]
def megabytes2bytes(mb):
    bytes_out = mb * 1024 * 1024
    return bytes_out


# Build full path for one or more files; return delimited string for multiple files
def build_abs_path(files, abs_path, sep=" "):
    if isinstance(files, str):
        file_abs_path = f"{abs_path}/{files}"
    elif isinstance(whitelist, list):
        file_abs_path = sep.join([f"{abs_path}/{f}" for f in files])


def unlist(lst, unique=False):
    if not isinstance(lst, list) and not isinstance(lst, dict):
        # not a list or dict; return as a list
        return [lst]
    elif isinstance(lst, list) and len(lst) == 1:
        # list w/ 1 entry; return input as is
        return lst
    elif isinstance(lst, list) and len(lst) > 1:
        # nested lists- untangle
        result = []
        for item in lst:
            if isinstance(item, list):
                result.extend(unlist(item))
            else:
                result.append(item)
    elif isinstance(lst, dict):
        # dictionary; untangle values
        result = []
        for value in lst.values():
            if isinstance(value, (list, dict)):
                result.extend(unlist(value))
            else:
                result.append(value)

    if unique:
        result = list(set(result))

    return result


# USed for finding max BC length
def max_sum_of_entries(lst):
    max_sum = 0
    current_sum = 0

    for item in lst:
        try:
            # Try converting the item to int; if successful, it's a single number
            num = int(item)
            current_sum = num
        except ValueError:
            # If conversion fails, it's a string representing a list of numbers
            numbers = list(map(int, item.split()))
            current_sum = sum(numbers)

        if current_sum > max_sum:
            max_sum = current_sum

    return max_sum


# Helper functions for barcode processing
def get_barcode_config(w):
    """
    Get barcode configuration for a given recipe.

    Returns a dictionary with barcode processing parameters.
    """
    try:
        recipe_info = {}

        # Get basic barcode info
        recipe_info["bc_lengths"] = get_recipe_info(w, "BC_length", mode="list")
        recipe_info["bc_concat"] = get_recipe_info(w, "BC_concat", mode="list")

        # Parse barcode lengths
        if recipe_info["bc_lengths"] and isinstance(recipe_info["bc_lengths"][0], str):
            recipe_info["bc_lengths_parsed"] = [
                int(x) for x in recipe_info["bc_lengths"][0].split()
            ]
        else:
            recipe_info["bc_lengths_parsed"] = recipe_info["bc_lengths"]

        # Determine if barcodes should be split
        recipe_info["needs_splitting"] = (
            recipe_info["bc_lengths_parsed"]
            and len(recipe_info["bc_lengths_parsed"]) >= 2
            and recipe_info["bc_concat"]
            and str(recipe_info["bc_concat"][0]).lower() == "true"
        )

        return recipe_info

    except Exception as e:
        print(f"Error getting barcode config: {e}")
        return {}


def get_recipe_barcode_strategy(w):
    """
    Determine the barcode processing strategy for a recipe.

    Returns one of: 'split', 'single', 'none'
    """
    try:
        recipes = get_recipes(w, mode="list")
        recipe_str = "".join(recipes) if recipes else ""

        # Check if this is a recipe that needs barcode splitting
        split_recipes = ["seeker", "decoder", "miST"]
        if any(recipe in recipe_str for recipe in split_recipes):
            config = get_barcode_config(w)
            if config.get("needs_splitting", False):
                return "split"

        return "single"

    except Exception as e:
        print(f"Error determining barcode strategy: {e}")
        return "single"


def validate_recipe_compatibility(SAMPLE_SHEET, RECIPE_SHEET):
    """
    Validate that recipes specified in sample sheet are compatible and have all required parameters.
    """
    compatibility_errors = []

    for index, row in SAMPLE_SHEET.iterrows():
        sample_id = row["sampleID"]

        # Check short read recipes
        if row["recipe"]:
            recipes = row["recipe"].split()

            # Check for incompatible recipe combinations
            incompatible_combinations = [
                (
                    ["rRNA-bwa", "ribodetector"],
                    "Cannot use both BWA and ribodetector for rRNA filtering",
                ),
                (
                    ["visium", "seeker"],
                    "Visium and Seeker are incompatible spatial technologies",
                ),
                (
                    ["stomics", "seeker"],
                    "Stomics and Seeker are incompatible spatial technologies",
                ),
                (
                    ["visium", "stomics"],
                    "Visium and Stomics are incompatible spatial technologies",
                ),
            ]

            for incompatible_recipes, error_msg in incompatible_combinations:
                if all(recipe in recipes for recipe in incompatible_recipes):
                    compatibility_errors.append(
                        (sample_id, f"Short read recipes: {error_msg}")
                    )

            # Check if recipes exist in recipe sheet
            for recipe in recipes:
                if recipe not in RECIPE_SHEET.index:
                    compatibility_errors.append(
                        (
                            sample_id,
                            f"Short read recipe '{recipe}' not found in recipe sheet",
                        )
                    )

        # Check ONT recipes
        if row["recipe_ONT"]:
            ont_recipes = row["recipe_ONT"].split()

            # Check for incompatible ONT recipe combinations
            ont_incompatible_combinations = [
                (
                    ["visium", "seeker"],
                    "Visium and Seeker are incompatible spatial technologies",
                ),
                (
                    ["stomics", "seeker"],
                    "Stomics and Seeker are incompatible spatial technologies",
                ),
                (
                    ["visium", "stomics"],
                    "Visium and Stomics are incompatible spatial technologies",
                ),
            ]

            for incompatible_recipes, error_msg in ont_incompatible_combinations:
                if all(recipe in ont_recipes for recipe in incompatible_recipes):
                    compatibility_errors.append(
                        (sample_id, f"ONT recipes: {error_msg}")
                    )

            # Check if ONT recipes exist in recipe sheet
            for recipe in ont_recipes:
                if recipe not in RECIPE_SHEET.index:
                    compatibility_errors.append(
                        (sample_id, f"ONT recipe '{recipe}' not found in recipe sheet")
                    )

        # Check cross-platform compatibility
        if row["recipe"] and row["recipe_ONT"]:
            short_recipes = row["recipe"].split()
            ont_recipes = row["recipe_ONT"].split()

            # Warn about potential spatial technology mismatches
            spatial_techs = ["visium", "stomics", "seeker", "decoder", "miST"]
            short_spatial = [r for r in short_recipes if r in spatial_techs]
            ont_spatial = [r for r in ont_recipes if r in spatial_techs]

            if short_spatial and ont_spatial and short_spatial != ont_spatial:
                compatibility_errors.append(
                    (
                        sample_id,
                        f"Spatial technology mismatch: short read uses {short_spatial}, ONT uses {ont_spatial}",
                    )
                )

    if compatibility_errors:
        print("ERROR: Recipe compatibility issues found:")
        for sample_id, error in compatibility_errors:
            print(f"Sample {sample_id}: {error}")
        raise ValueError(
            "Recipe compatibility errors found. Please fix your sample sheet."
        )

    return True


def comprehensive_sample_sheet_validation(
    SAMPLE_SHEET, RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT
):
    """
    Comprehensive validation of sample sheet, recipe sheet, and their cross-references.
    """
    print("=== Starting comprehensive sample sheet validation ===")

    # Step 1: Basic sample sheet validation
    print("1. Checking sample sheet structure and file requirements...")
    check_sample_sheet(SAMPLE_SHEET)

    # Step 2: Recipe sheet validation
    print("2. Checking recipe sheet structure and parameters...")
    check_recipe_sheet(RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT)

    # Step 3: Recipe compatibility validation
    print("3. Checking recipe compatibility and cross-references...")
    validate_recipe_compatibility(SAMPLE_SHEET, RECIPE_SHEET)

    # Step 4: Advanced cross-validation
    print("4. Performing advanced cross-validation...")
    advanced_validation_errors = []

    for index, row in SAMPLE_SHEET.iterrows():
        sample_id = row["sampleID"]

        # Check if rRNA filtering is specified but cdna_fa is missing
        if row["recipe"]:
            recipes = row["recipe"].split()
            if "rRNA-bwa" in recipes:
                if (
                    not row.get("cdna_fa")
                    or str(row["cdna_fa"]).strip() == ""
                    or str(row["cdna_fa"]) == "nan"
                ):
                    advanced_validation_errors.append(
                        (
                            sample_id,
                            "Recipe 'rRNA-bwa' requires 'cdna_fa' field but it's missing or empty",
                        )
                    )

        if row["recipe_ONT"]:
            ont_recipes = row["recipe_ONT"].split()
            # Most ONT workflows need cdna_fa for gene annotation
            ont_needs_cdna = ["visium", "seeker", "stomics", "decoder", "miST"]
            if any(recipe in ont_needs_cdna for recipe in ont_recipes):
                if (
                    not row.get("cdna_fa")
                    or str(row["cdna_fa"]).strip() == ""
                    or str(row["cdna_fa"]) == "nan"
                ):
                    advanced_validation_errors.append(
                        (
                            sample_id,
                            f"ONT recipes {ont_needs_cdna} require 'cdna_fa' field but it's missing or empty",
                        )
                    )

        # Check barcode/UMI consistency
        if row["recipe"]:
            recipes = row["recipe"].split()
            spatial_recipes = ["visium", "seeker", "stomics", "decoder", "miST"]
            if any(recipe in spatial_recipes for recipe in recipes):
                # These recipes should have consistent barcode configurations
                for recipe in recipes:
                    if recipe in spatial_recipes and recipe in RECIPE_SHEET.index:
                        bc_length = (
                            RECIPE_SHEET.loc[recipe, "BC_length"]
                            if "BC_length" in RECIPE_SHEET.columns
                            else None
                        )
                        umi_length = (
                            RECIPE_SHEET.loc[recipe, "UMI_length"]
                            if "UMI_length" in RECIPE_SHEET.columns
                            else None
                        )

                        if pd.isna(bc_length) or str(bc_length).strip() == "":
                            advanced_validation_errors.append(
                                (
                                    sample_id,
                                    f"Spatial recipe '{recipe}' has empty BC_length",
                                )
                            )
                        if pd.isna(umi_length) or str(umi_length).strip() == "":
                            advanced_validation_errors.append(
                                (
                                    sample_id,
                                    f"Spatial recipe '{recipe}' has empty UMI_length",
                                )
                            )

    if advanced_validation_errors:
        print("ERROR: Advanced validation errors found:")
        for sample_id, error in advanced_validation_errors:
            print(f"Sample {sample_id}: {error}")
        raise ValueError(
            "Advanced validation errors found. Please fix your configuration."
        )

    print("=== Sample sheet validation completed successfully ===")
    return True
