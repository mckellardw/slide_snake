# Recipes for `slide_snake`

The default receipe sheet is stored in `resources/recipe_sheet.csv`. Change the path in `config/config.yaml` if you would like to point to a custom version. If you want to add a custom recipe, read below!

## **Recipe descriptions**:
"Recipes" are descriptions for the alignment workflow- how to trim the barcode read, whether or not to filter out ribosomal RNA reads, additional alignment parameters, etc. One goal of `slide_snake` is to make the alignment preprocessing modular so that all of these parameters can be compared directly and rigorously. Please note the following:
  - Multiple recipes can be passed for each sample. Include as many as you would like, each separated by a space in the sample sheet.
  - You can also add a new recipe! Just add a new line to `resources/recipe_sheet.csv` and give it a unique name in the 1st (0th for you pythoners) column. 
  - Recipe naming convention:
    - Include the preprocessing steps that you want to use- i.e., including "rRNA.bwa" in the name means that rRNA filtering with bwa alignment will be done prior to genomic alignment. 

## `kallisto` technology speifications
`kb_python` uses `ngs-tools` for tech specifications - [link](https://github.com/Lioscro/ngs-tools/tree/aa3e864e59ae78467a331f671967c93d62a6e2ad)
```
name            description                            on-list    barcode                    umi        cDNA
------------    -----------------------------------    -------    -----------------------    -------    -----------------------
10XV1           10x version 1                          yes        0,0,14                     1,0,10     2,None,None
10XV2           10x version 2                          yes        0,0,16                     0,16,26    1,None,None
10XV3           10x version 3                          yes        0,0,16                     0,16,28    1,None,None
10XV3_ULTIMA    10x version 3 sequenced with Ultima    yes        0,22,38                    0,38,50    0,62,None
BDWTA           BD Rhapsody                            yes        0,0,9 0,21,30 0,43,52      0,52,60    1,None,None
BULK            Bulk (single or paired)                                                                 0,None,None 1,None,None
CELSEQ          CEL-Seq                                           0,0,8                      0,8,12     1,None,None
CELSEQ2         CEL-SEQ version 2                                 0,6,12                     0,0,6      1,None,None
DROPSEQ         DropSeq                                           0,0,12                     0,12,20    1,None,None
INDROPSV1       inDrops version 1                                 0,0,11 0,30,38             0,42,48    1,None,None
INDROPSV2       inDrops version 2                                 1,0,11 1,30,38             1,42,48    0,None,None
INDROPSV3       inDrops version 3                      yes        0,0,8 1,0,8                1,8,14     2,None,None
SCRUBSEQ        SCRB-Seq                                          0,0,6                      0,6,16     1,None,None
SMARTSEQ2       Smart-seq2  (single or paired)                                                          0,None,None 1,None,None
SMARTSEQ3       Smart-seq3                                                                   0,11,19    0,11,None 1,None,None
SPLIT-SEQ       SPLiT-seq                                         1,10,18 1,48,56 1,78,86    1,0,10     0,None,None
STORMSEQ        STORM-seq                                                                    1,0,8      0,None,None 1,14,None
SURECELL        SureCell for ddSEQ                                0,0,6 0,21,27 0,42,48      0,51,59    1,None,None
Visium          10x Visium                             yes        0,0,16                     0,16,28    1,None,None
```

Info on custom technology string from [kb-python preprint](https://www.biorxiv.org/content/10.1101/2023.11.21.568164v2.full.pdf):
```
The custom technology string (supplied to -x) contains the format barcode:UMI:DNA,
representing the locational information of the barcode, UMI, and the DNA (where DNA is the
biological read to be mapped):

-x a,b,c:d,e,f:g,h,i
● a: barcode file number, b: barcode start position, c: barcode end position
● d: UMI file number, e: UMI start position, f: UMI end position
● g: DNA file number, h: DNA start position, i: DNA end position

Important notes: File numbers and positions are zero-indexed. If no specific end position
exists (i.e. the end position is the very end of the read), the end position should be set to 0. If
cell barcodes and/or UMIs are not supported by the technology, the barcode and/or UMI field
can be set to -1,0,0.

Thus, for 10xv3:
-x 0,0,16:0,16,28:1,0,0

Sequences can be stitched together by specifying multiple locations; for example, a
SPLiT-seq45 assay, which contains three separate unlinked barcodes, each of length 8, and a
UMI of length 10 in the second file and the DNA in the first file would look as follows:

-x 1,10,18,1,48,56,1,78,86:1,0,10:0,0,0

Final note about multiple locations: If the paired-end read mapping option is enabled, exactly
two DNA locations should be specified (for the first and second read in the pair).
If a technology does not fit into this format (e.g. due to barcodes or UMIs of variable lengths
and positions), preprocessing of the FASTQ file should be performed beforehand to reformat
the reads into a structure that can be handled by this format.
```

<details close>
<summary> SlideSeq/Seeker (Curio) </summary>
[manuscript link]()
[Curio website](https://curiobioscience.com/)

Because of read quality issues (indels, low Q scores, etc.) in the SlideSeq barcode read, I have added a few custom strategies for handling these data:
- `seeker` - No hard trimming, and use the base positions for barcode/UMI (*Note*, this recipe doesn't work well w/ Curio Seeker b/c of in/del issues w/ the barcode synthesis)
- `seeker_hardTrim` - Hard trim the adapter read positions in R1, and use the best barcode correction algorithms in STARsolo
- `seeker_MatchLinker` - Match the adapter sequence on R1 (w/ 2 mismatches allowed) and infer barcodes/UMIs from that position (*Note* best performer w/ Curio data)
- `seeker_MatchLinker_total` - Same as `seeker_noTrimMatchLinker`, but with additional STAR parameters for total RNAseq alignment (more multimappers, looser alignment)
</details>

<details close>
<summary> StereoSeq/STOmics (BGI) </summary>
  
- `stomics` - Standard alignment for StereoSeq/STOmics (BGI) data
- `stomics_total` - Total RNA alignment for StereoSeq/STOmics (BGI) data
- `stomics_rRNA.STAR` - Standard alignment performed after filtering rRNA with STAR alignment
- `stomics_total_rRNA.STAR` - Total RNA alignment performed after filtering rRNA with STAR alignment

</details>


<details close>
<summary> Visium (10x Genomics) </summary>

- `visium` - #description
- `visium_total` - #description
- `visium_total_rRNA.STAR` - #description

</details>


<details close>
<summary> DecoderSeq (Cao et al, Nat Biotechnol, 2024) </summary>
[manuscript link](https://www.nature.com/articles/s41587-023-02086-y)  

- `decoder` - #description
- `decoder_total` - #description
- `decoder_total_rRNA.STAR` - #description

</details>
