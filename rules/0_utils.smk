import os
import re
import glob

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
    shell:
        """
        samtools index -@ {threads} {input.BAM}
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
        "genes_gtf": [".gtf"]
    }

    ONT_requirements = {
        "ONT": [".fq.gz", ".fastq.gz", ".bam", ".cram"],
        "genome_fa": [".fa"],
        "genes_gtf": [".gtf"],
        "cdna_fa": [".fa"]
    }

    missing_files = []
    empty_fields = []
    duplicate_samples = SAMPLE_SHEET["sampleID"].duplicated()
    duplicate_sample_names = SAMPLE_SHEET["sampleID"][duplicate_samples].tolist()

    for index, row in SAMPLE_SHEET.iterrows():
        sample_id = row["sampleID"]

        # Check additional requirements based on recipe
        if row["recipe"]:
            for column, exts in short_read_requirements.items():
                if column not in row:
                    continue
                file_paths = row[column].split()
                if not file_paths:
                    empty_fields.append((sample_id, column))
                for file_path in file_paths:
                    if file_path:
                        if column == "STAR_ref" and not os.path.isdir(file_path):
                            print(f"WARNING: {file_path} for sample {sample_id} was not found, or is not a directory.")
                            missing_files.append((sample_id, file_path))
                        elif not any(re.search(f"{ext}$", file_path) for ext in exts):
                            print(f"WARNING: File {file_path} for sample {sample_id} does not have the correct extension {exts}")
                        elif not glob.glob(file_path) and not re.search(r"\{.*\}", file_path):
                            missing_files.append((sample_id, file_path))

        if row["recipe_ONT"]:
            for column, exts in ONT_requirements.items():
                if column not in row:
                    continue
                file_paths = row[column].split()
                if not file_paths:
                    empty_fields.append((sample_id, column))
                for file_path in file_paths:
                    if file_path and not any(re.search(f"{ext}$", file_path) for ext in exts):
                        print(f"WARNING: File {file_path} for sample {sample_id} does not have the correct extension {exts}")

    if empty_fields:
        print("WARNING: The following required fields are empty:")
        for sample_id, column in empty_fields:
            print(f"Sample {sample_id}: {column}")

    if duplicate_samples.any():
        print("WARNING: Duplicate sample IDs found in the sample sheet:")
        for sample_id in duplicate_sample_names:
            print(f"Sample {sample_id}")

    if missing_files:
        print("WARNING: The following files are missing:")
        for sample_id, file_path in missing_files:
            print(f"Sample {sample_id}: {file_path}")
        raise FileNotFoundError("Some required files are missing. Please check the log for details.")


def check_recipe_sheet(RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT):
    required_columns = ["whitelist", "R1_finalLength", "fwd_primer", "rev_primer", "BC_adapter", "BC_length", "BC_offset", "BC_position", "BC_max_ED", "BC_min_ED_diff", "BC_concat", "UMI_adapter", "UMI_length", "UMI_offset", "UMI_position", "internal_adapter", "STAR_soloType", "STAR_soloCBmatchWLtype", "STAR_soloCB", "STAR_soloUMI", "STAR_soloAdapter", "STAR_extra", "kb_x", "kb_extra", "featureCounts_extra", "mm2_extra"]
    missing_columns = [col for col in required_columns if col not in RECIPE_SHEET.columns]
    duplicate_recipes = RECIPE_SHEET.index.duplicated()
    duplicate_recipe_names = RECIPE_SHEET.index[duplicate_recipes].tolist()

    if missing_columns:
        print("WARNING: The following required columns are missing from the recipe sheet:")
        for col in missing_columns:
            print(f"- {col}")

    if duplicate_recipes.any():
        print("WARNING: Duplicate recipes found in the recipe sheet:")
        for recipe in duplicate_recipe_names:
            print(f"Recipe {recipe}")

    RECIPES = unlist(RECIPE_DICT, unique=True) + unlist(RECIPE_ONT_DICT, unique=True)
    for RECIPE in RECIPES:
        if RECIPE not in list(RECIPE_SHEET.index):
            print(f"WARNING: `{RECIPE}` not found in recipe sheet!")
            print(f"Fix RECIPE_SHEET [{config['RECIPE_SHEET']}] and try again!")


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

                # TODO - update to match ribodetector style

            elif "ribodetector" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/short_read/rRNA/ribodetector/no_rRNA_R2.fq.gz"

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


def get_bc_adapter(w, mode=["ONT", "ILMN"]):
    try:
        if mode == "ILMN":
            if "Trim" in w.RECIPE or "std" in w.RECIPE:
                return get_recipe_info(w, "BC_start") # return integer start positions for each BC sequence
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
                return get_recipe_info(w, "UMI_start") # return integer start positions for each BC sequence
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
        elif "matchLinker" in w.RECIPE or "seeker" in w.RECIPE:
            # Use BC_concat=False for combinatorial barcode constructs (DBIT, microST, etc)
            if return_type == "list":
                if RECIPE_SHEET["BC_concat"][w.RECIPE] and mode in ["STAR", "ILMN"]:
                    # Barcode constructs where positional barcodes are NOT independent (must be concatenated)
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt",
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt",
                    ]
                elif RECIPE_SHEET["BC_concat"][w.RECIPE] and mode=="ONT":
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
            elif "microST" in recipe:
                star_info[key] = RECIPE_SHEET[key]["microST_ligation_total"]
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
    with open(file_path, 'r') as file:
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

