#### Utility rules ############################2


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


# Check to see if recipes are in the recipe sheet; does not kill run, only prints warnings.
def check_recipe_sheet(RECIPE_SHEET, RECIPE_DICT, RECIPE_ONT_DICT):
    RECIPES = unlist(RECIPE_DICT, unique=True) + unlist(RECIPE_ONT_DICT, unique=True)
    for RECIPE in RECIPES:
        if RECIPE not in list(RECIPE_SHEET.index):
            print(f"WARNING: `{RECIPE}` not found in recipe sheet!")
            print(f"Fix RECIPE_SHEET [{config['RECIPE_SHEET']}] and try again!")


def check_sample_sheet(SAMPLE_SHEET):
    print("TODO")


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
            if "rRNA.STAR" in w.RECIPE:  # Use trimmed & STAR-rRNA-filtered .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/noRibo_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/noRibo_R2.fq.gz"

                # TODO - update to match ribodetector style

            elif "rRNA.bwa" in w.RECIPE:  # Use trimmed & bwa-rRNA-filtered .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/noRibo_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/noRibo_R2.fq.gz"

                # TODO - update to match ribodetector style

            elif "ribodetector" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/ribodetector/noRibo_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/ribodetector/noRibo_R2.fq.gz"

                if "internalTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/ribodetector/noRibo_internalTrim_R1.fq.gz"
                if "hardTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/ribodetector/noRibo_hardTrim_R1.fq.gz"
            else:  # just trimmed .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R2.fq.gz"

                if "internalTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_internalTrim_R1.fq.gz"
                if "hardTrim" in w.RECIPE:
                    R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_hardTrim_R1.fq.gz"
        elif mode == "ONT":
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_R2.fq.gz"
            if "internalTrim" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_internalTrim_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/twiceCut_internalTrim_R2.fq.gz"
            elif "hardTrim" in w.RECIPE:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_hardTrim_R1.fq.gz"
            else:
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_R1.fq.gz"
        else:
            print("get_fqs(): `mode` not found")
    except Exception:
        if mode == "ILMN":
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R2.fq.gz"
        elif mode == "ONT":
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/ont/cut_R2.fq.gz"
        else:
            print("get_fqs(): `mode` not found")

    # return fq path(s)
    if return_type == "list":
        return [R1, R2]
    elif return_type == "dict":
        return {"R1": R1, "R2": R2}


# whitelist param handling for different recipes/technologies/chemistries/etc
def get_whitelist(w, return_type=None, mode=["recipe", "all_used", "all"]):
    try:
        if "matchLinker" in w.RECIPE:
            if return_type == "list":
                if RECIPE_SHEET["BC_concat"][w.RECIPE]:
                    # Barcode constructs where positional barcodes are NOT independent (must be concatenated)
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt",
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt",
                    ]
                else:
                    # Barcode constructs where positional barcodes ARE independent (shorter white list)
                    whitelist = [
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_uniq_1.txt",
                        f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_uniq_2.txt",
                    ]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_1.txt {w.OUTDIR}/{w.SAMPLE}/bc/whitelist_2.txt"
        elif "internalTrim" in w.RECIPE:
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist.txt"
        elif "adapterInsert" in w.RECIPE:
            if return_type == "list":
                whitelist = [f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_adapter.txt"]
            else:
                whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bc/whitelist_adapter.txt"
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
def get_n_bcs(w):
    bc_lengths = get_recipe_info(w, info_col="BC_length").split()
    n_bcs = len(bc_lengths)
    return n_bcs


# TODO move these values to recipe_sheet - also write better code than this...
def get_ont_barcode_pattern(w):
    ## SlideSeq/Seeker: R1="C"*22 | BC1="C"*8 | Linker="C"*18 | BC2="C"*6 | UMI="N"* 7
    try:
        if "stomics" in w.RECIPE:
            BC_PATTERN = "C" * 25 + "N" * 10
        elif "visium" in w.RECIPE:
            BC_PATTERN = "C" * 16 + "N" * 12
        elif "microST_ligation" in w.RECIPE:
            BC_PATTERN = "C" * 34
        elif "microST_klenow_v1" in w.RECIPE:
            BC_PATTERN = "C" * 34 + "N" * 12
        elif "seeker" in w.RECIPE and "matchLinker" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 6 + "N" * 7
        elif "seeker" in w.RECIPE and "internalTrim" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 6 + "N" * 7
        elif "seeker" in w.RECIPE and "adapterInsert" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 18 + "C" * 6 + "N" * 7
        else:
            BC_PATTERN = "C" * 16 + "N" * 12
    except Exception:
        BC_PATTERN = "C" * 16 + "N" * 12

    # return barcode pattern
    return BC_PATTERN

    # BC_PATTERN="(?P<discard_1>CTACACGACGCTCTTCCGATCT)"+ \
    #     "(?P<cell_1>.{{8}})"+ \
    #     "(?P<discard_2>TCTTCAGCGTTCCCGAGA)"+ \
    #     "(?P<cell_2>.{{6}})"+ \
    #     "(?P<umi_1>.{{7}})"

    # BC_PATTERN="(?P<discard_1>XXXXXXXXXXXXXXXXXXXXXXX)"+ \
    #     "(?P<cell_1>.{{8}})"+ \
    #     "(?P<discard_2>XXXXXXXXXXXXXXXXXX)"+ \
    #     "(?P<cell_2>.{{6}})"+ \
    #     "(?P<umi_1>.{{7}})"


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


def get_STAR_ref(w, mode=["genome", "rRNA"]):
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
        elif mode == "rRNA":
            out = SAMPLE_SHEET["STAR_rRNA_ref"][w.SAMPLE]
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
        "STAR_soloCBmatchWLtype": "",
        "STAR_soloAdapter": "",
        "STAR_extra": "--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.05  --outFilterMatchNmin 12  --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0",
    }

    # Iterate over each key
    for key in star_info.keys():
        try:
            star_info[key] = RECIPE_SHEET[key][w.RECIPE]
        except Exception:
            # values for rRNA/default
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
    if isinstance(mode, list) and len(mode) > 1:  # default
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


# Convert a megabyte value (str) to bytes [int]
def megabytes2bytes(mb):
    # mb_int = int(''.join(filter(str.isdigit, mb)))
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


# def unlist(*args, unique=False):
#     result = []
#     for arg in args:
#         if isinstance(arg, (list, dict)):
#             if isinstance(arg, list):
#                 result.extend(unlist(arg))
#             elif isinstance(arg, dict):
#                 result.extend(unlist(list(arg.values())))
#         else:
#             result.append(arg)

#     if unique:
#         result = list(set(result))  # Convert to set to remove duplicates, then back to list

#     return result

# def unlist(*args, unique=False):
#     def process_item(item):
#         if isinstance(item, (list, dict)):
#             if isinstance(item, list):
#                 yield from unlist(item, unique=unique)
#             elif isinstance(item, dict):
#                 yield from unlist(list(item.values()), unique=unique)
#         else:
#             yield item

#     result = []
#     for arg in args:
#         result.extend(process_item(arg))

#     if unique:
#         result = list(set(result))  # Convert to set to remove duplicates, then back to list

#     return result


########## TO BE DEPRECATED ########################################################
def get_barcode_length(w):
    """
    Get barcode length based on the recipe(s) passed.
    """
    bc_lengths = [RECIPE_SHEET["BC_length"][R] for R in RECIPE_DICT[w.SAMPLE]]
    if len(bc_lengths) == 1:
        out = bc_lengths[0]
    else:
        # TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(bc_lengths), key=bc_lengths.count)  # return mode
    return out


def get_umi_length(w):
    """
    Get UMI length based on the recipe(s) passed.
    """
    umi_lengths = [RECIPE_SHEET["UMI_length"][R] for R in RECIPE_DICT[w.SAMPLE]]
    if len(umi_lengths) == 1:
        out = umi_lengths[0]
    else:
        # TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(umi_lengths), key=umi_lengths.count)  # return mode
    return out


def get_split_ont_align_mem_gb(wildcards, threads):
    return config["RESOURCES_MM2_MEM_GB"] / threads
