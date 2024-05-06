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
# whitelist param handling for different recipes
def get_whitelist(w):
    try:
        if "stomics" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"
        elif "noTrim" in w.RECIPE or "matchLinker" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist_1.txt {w.OUTDIR}/{w.SAMPLE}/bb/whitelist_2.txt"
        elif "internalTrim" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"
        elif "adapterInsert" in w.RECIPE:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist_adapter.txt"
        else:
            whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"
    except:
        whitelist = f"{w.OUTDIR}/{w.SAMPLE}/bb/whitelist.txt"

    # return whitelist path(s)
    return whitelist


# Select input reads based on alignment recipe
def get_fqs(w):
    try:
        if "rRNA.STAR" in w.RECIPE:  # Use trimmed & STAR-rRNA-filtered .fq's
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/final_filtered_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/STARsolo/final_filtered_R2.fq.gz"
        elif "rRNA.bwa" in w.RECIPE:  # TODO Use trimmed & bwa-rRNA-filtered .fq's
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/final_filtered_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/rRNA/bwa/final_filtered_R2.fq.gz"
        elif "rRNA" not in w.RECIPE:  # just trimmed .fq's
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R1.fq.gz"
            R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R2.fq.gz"
        else:
            print("I just don't know what to do with myself...")

        if "internalTrim" in w.RECIPE or "hardTrim" in w.RECIPE:
            R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_trimmed_R1.fq.gz"
    except:
        R1 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R1.fq.gz"
        R2 = f"{w.OUTDIR}/{w.SAMPLE}/tmp/twiceCut_R2.fq.gz"

    # return fq path(s) as a list
    return [R1, R2]


# TODO move these values to recipe_sheet
def get_ont_barcode_pattern(w):
    ## SlideSeq/Seeker: R1="C"*22 | BC1="C"*8 | Linker="C"*18 | BC2="C"*6 | UMI="N"* 7
    try:
        if "stomics" in w.RECIPE:
            BC_PATTERN = "C" * 25 + "N" * 10
        elif "visium" in w.RECIPE:
            BC_PATTERN = "C" * 16 + "N" * 12
        elif "seeker" in w.RECIPE and "noTrim" in w.RECIPE or "matchLinker" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 6 + "N" * 7
        elif "seeker" in w.RECIPE and "internalTrim" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 6 + "N" * 7
        elif "seeker" in w.RECIPE and "adapterInsert" in w.RECIPE:
            BC_PATTERN = "C" * 8 + "C" * 18 + "C" * 6 + "N" * 7
        else:
            BC_PATTERN = "C" * 16 + "N" * 12
    except:
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


def get_STAR_ref(w, mode="genome"):
    try:
        if mode == "genome":
            star_ref = REF_DICT[w.SAMPLE]
        elif mode == "rRNA":
            star_ref = rRNA_STAR_DICT[w.SAMPLE]
    except:
        star_ref = "No reference given! Check your sample sheet!"

    return star_ref

#TODO- add STAR_rRNA functinoality (no w.RECIPE accessibility)
def get_STAR_extra_params(w):
    star_info = {        
        "STAR.soloType":"",
        "STAR.soloUMI":"",
        "STAR.soloCB":"",
        "STAR.soloCBmatchWLtype":"",
        "STAR.soloAdapter":"",
        "STAR.extra":"--outFilterMultimapNmax 50 --outFilterMismatchNoverLmax 0.05  --outFilterMatchNmin 12  --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0",
    }

    # Iterate over each key
    for key in star_info.keys():
        try:
            star_info[key] = RECIPE_SHEET[key][w.RECIPE]
        except:
            # values for rRNA/default?
            recipe = get_recipe(w,mode="concatenate")

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


def get_recipe(w, mode="ONT"):
    if mode == "ONT":
        try:
            return RECIPE_ONT_DICT[w.SAMPLE]
        except:
            print("No ONT recipe given! Check your sample sheet!")
            return ""
    elif mode == "illumina":
        try:
            return RECIPE_DICT[w.SAMPLE]
        except:
            print("No ILMN recipe given! Check your sample sheet!")
            return ""
    else:
        out = []
        if "RECIPE_ONT_DICT" in globals():
            if w.SAMPLE in RECIPE_ONT_DICT:
                out.extend(RECIPE_ONT_DICT[w.SAMPLE])
        elif "RECIPE_DICT" in globals():
            if w.SAMPLE in RECIPE_DICT:
                out.extend(RECIPE_DICT[w.SAMPLE])
        if len(out) > 0:
            if mode == "list":
                return out
            else:
                return " ".join(out)
        else:
            print(f"No recipe found for {w.SAMPLE}! Check your sample sheet!")
            return ""


# Pull info from recipe sheet
def get_recipe_info(recipe, info_col):
    return RECIPE_SHEET[info_col][recipe]


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
