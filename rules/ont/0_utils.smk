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
