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
    return(whitelist)


# Build full path for one or more files; return delimited string for multiple files
def build_abs_path(files, abs_path, sep=" "):
    if isinstance(files, str):
        file_abs_path = f"{abs_path}/{files}"
    elif isinstance(whitelist, list):
        file_abs_path = sep.join([f"{abs_path}/{f}" for f in files])


# Pull info from recipe sheet
def get_recipe_info(recipe, info_col):
    return(RECIPE_SHEET[recipe, info_col])


### helper functions ########################################################################
# Borrowed from sockeye - https://github.com/nanoporetech/sockeye
def get_barcode_length(w):
    """
    Get barcode length based on the recipe(s) passed.
    """
    bc_lengths = [RECIPE_SHEET["BC.length"][R] for R in RECIPE_DICT[w.SAMPLE]]
    if len(bc_lengths) == 1:
        out = bc_lengths[0]
    else:
        #TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(bc_lengths), key=bc_lengths.count) # return mode
    return out


def get_umi_length(w):
    """
    Get UMI length based on the recipe(s) passed.
    """
    umi_lengths = [RECIPE_SHEET["UMI.length"][R] for R in RECIPE_DICT[w.SAMPLE]]
    if len(umi_lengths) == 1:
        out = umi_lengths[0]
    else:
        #TODO: there is probably a better way to handle multi-recipe than this
        out = max(set(umi_lengths), key=umi_lengths.count) # return mode
    return out


def get_split_ont_align_mem_gb(wildcards, threads):
    return config["RESOURCES_MM2_MEM_GB"] / threads
