#### Utility rules ############################2


# Index .bam file
rule index_BAM:
    input:
        BAM="{BAM}",
    output:
        BAI="{BAM}.bai",
    # wildcard_constraints:
    #     BAM=".*\.(bam)$"
    threads: config["CORES"]
    run:
        shell(
            f"""
            {EXEC['SAMTOOLS']} index -@ {threads} {input.BAM}
            """
        )


#### Util functions ###########################


# Select input reads based on alignment recipe
def get_fqs(w, return_type=["list", "dict"], mode=["ONT", "ILMN"]):
    # param defaults
    if len(return_type) > 1:  # default
        return_type = "list"
    if type(mode) is list:  # default
        mode = "ONT"

    # get file paths
    try:
        if mode == "ILMN":
            if "rRNA.STAR" in w.RECIPE:  # Use trimmed & STAR-rRNA-filtered .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz"
            elif "rRNA.bwa" in w.RECIPE:  # Use trimmed & bwa-rRNA-filtered .fq's
                R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz"
                R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz"
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


# whitelist param handling for different recipes
def get_whitelist(w):
    try:
        if "noTrim" in w.RECIPE or "matchLinker" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist_1.txt {w.OUTDIR}/{w.SAMPLE}/bb/whitelist_2.txt"
        elif "internalTrim" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"
        elif "adapterInsert" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist_adapter.txt"
        else:
            # visium, stomics, microST
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"
    except Exception:
        whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"

    # return whitelist path(s)
    return whitelist


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


def get_STAR_ref(w, mode=["genome", "rRNA"]):
    try:
        if len(mode) > 1:  # default
            mode = "genome"

        if mode == "genome":
            star_ref = REF_DICT[w.SAMPLE]
        elif mode == "rRNA":
            star_ref = rRNA_STAR_DICT[w.SAMPLE]
    except Exception:
        star_ref = "No reference given! Check your sample sheet!"

    return star_ref


# TODO- add STAR_rRNA functinoality (no w.RECIPE accessibility)
def get_STAR_extra_params(w):
    star_info = {
        "STAR.soloType": "",
        "STAR.soloUMI": "",
        "STAR.soloCB": "",
        "STAR.soloCBmatchWLtype": "",
        "STAR.soloAdapter": "",
        "STAR.extra": "--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.05  --outFilterMatchNmin 12  --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0",
    }

    # Iterate over each key
    for key in star_info.keys():
        try:
            star_info[key] = RECIPE_SHEET[key][w.RECIPE]
        except Exception:
            # values for rRNA/default?
            recipe = get_recipes(w, mode=f"concatenate")

            if "stomics" in recipe:
                star_info[key] = RECIPE_SHEET[key]["stomics_total"]
            elif "visium" in recipe:
                star_info[key] = RECIPE_SHEET[key]["visium_total"]
            elif "seeker" in recipe and "noTrim" in recipe or "matchLinker" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_matchLinker_total"]
            elif "seeker" in recipe and "internalTrim" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_internalTrim_total"]
            elif "seeker" in recipe and "adapterInsert" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_adapterInsert_total"]
            elif "seeker" in recipe:
                star_info[key] = RECIPE_SHEET[key]["seeker_total"]
            else:
                star_info[key] = RECIPE_SHEET[key]["visium_total"]
        # except Exception:
        #     pass

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
    if len(mode) > 1:  # default
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
        print("elseesese")
        all_recipes = unlist(all_recipes)

        if len(all_recipes) > 0:
            if "list" in mode:
                return all_recipes
            else:
                print(42424242)
                return " ".join(all_recipes)
        else:
            print(f"No recipe found for {w.SAMPLE}! Check your sample sheet!")
            return ""


# Pull info from recipe sheet
def get_recipe_info(w, info_col, mode=["ONT", "ILMN"]):
    recipe = get_recipes(w, mode=mode)

    # if len(recipe) > 1:
    #     return [RECIPE_SHEET[info_col][x] for x in recipe]
    # else:
    #     return RECIPE_SHEET[info_col][recipe]
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


def unlist(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(unlist(item))
        else:
            result.append(item)
    return result


########## TO BE DEPRECATED ########################################################
def get_barcode_length(w):
    """
    Get barcode length based on the recipe(s) passed.
    """
    bc_lengths = [RECIPE_SHEET["BC.length"][R] for R in RECIPE_DICT[w.SAMPLE]]
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
    umi_lengths = [RECIPE_SHEET["UMI.length"][R] for R in RECIPE_DICT[w.SAMPLE]]
    if len(umi_lengths) == 1:
        out = umi_lengths[0]
    else:
        # TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(umi_lengths), key=umi_lengths.count)  # return mode
    return out


def get_split_ont_align_mem_gb(wildcards, threads):
    return config["RESOURCES_MM2_MEM_GB"] / threads
