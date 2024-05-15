# Recipes for `slide_snake`

The default receipe sheet is stored in `resources/recipe_sheet.csv`. Change the path in `config/config.yaml` if you would like to point to a custom version. If you want to add a custom recipe, read below!

## **Recipe descriptions**:
"Recipes" are descriptions for the alignment workflow- how to trim the barcode read, whether or not to filter out ribosomal RNA reads, additional alignment parameters, etc. One goal of `slide_snake` is to make the alignment preprocessing modular so that all of these parameters can be compared directly and rigorously. Please note the following:
  - Multiple recipes can be passed for each sample. Include as many as you would like, each separated by a space in the sample sheet.
  - You can also add a new recipe! Just add a new line to `resources/recipe_sheet.csv` and give it a unique name in the 1st (0th for you pythoners) column. 
  - Recipe naming convention:
    - Include the preprocessing steps that you want to use- i.e., including "rRNA.bwa" in the name means that rRNA filtering with bwa alignment will be done prior to genomic alignment. 


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
