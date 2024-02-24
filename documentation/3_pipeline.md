# Pipeline information  

## **Trimming:**
<details close>
<summary> Five Prime adapter(s): </summary>
- Reverse complement of the SlideSeq TSO to reduce strand invasion artifacts [**CCCATTCACTCTGCGTTGATACCAGCTT**]  
</details>

<details close>
<summary> Three Prime adapter(s): </summary>

- "A" homopolymers [**100-A**]  
- "G" homopolymers, important for Illumina 2-color sequencers [**100-G**]  
- "T" homopolymers [**100-T**]  
- Nextera adapter sequence [**CTGTCTCTTATA**]  
- Reverse complement of Nextera sequence [**TATAAGAGACAG**]  
- Curio template switch olgo (TSO) - remove any polyadenylated TSOs [**AAGCTGGTATCAACGCAGAGTGAATGGG**]  
- Curio R1 internal adapter - shows up in R2, and used for R1 trimming [**TCTTCAGCGTTCCCGAGA**]  
- Reverse of Curio R1 adapter [**AGAGCCCTTGCGACTTCT**]  
- Illumina unversal sequence [**AGATCGGAAGAG**]  
</details>

## **Alignment:**
<details close>
<summary> STAR </summary>
  
- After adapter trimming, reads are aligned with `STARsolo` and `kallisto`/`bustools` to generate count matrices  
- Outputs are in `{SAMPLE_ID}/STARsolo/Solo.out` & `{SAMPLE_ID}/kb/counts_unfiltered`  
- Different recipes are written out in `resources/recipe_sheet.csv`, and must be specified for each sample within the sample sheet  

### Generating references:
#### rRNA STAR reference for in silico rRNA depletion/quantification
Ribosomal RNA (rRNA) molecules can make alignment/quantification very difficult because of the number of genomic copies of these genes. We added a first-pass-alignment just to rRNA sequences to enable stratified parameterization for these sequences, but maintain the ability to count and analyze them.  

Check out `scripts/GRCm39_GENCODEM31_STAR_rRNA.sh` for an example script showing how to generate a rRNA-only STAR reference using GENCODE annotations.  

Other option- download rRNA sequences from  with the following query:
`(expert_db:"SILVA" AND TAXONOMY:"10090") AND entry_type:"Sequence"`

#### Genomic STAR reference
This is a typical STAR reference that you would use for any other alignment job. Here is an example code snippet:
```
FASTA_DIR="/path/to/GENCODE_M31/GRCm39.genome.fa"
GENES_DIR="/path/to//GENCODE_M31/gencode.vM31.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_DIR} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
```
*You can find the reference files on [GENCODE's website](https://www.gencodegenes.org/mouse/)*
</details>


<details close>
<summary> kallisto/BUStools </summary>
  
TODO
</details>  


## **Recipe descriptions**:
"Recipes" are descriptions for the alignment workflow- how to trim the barcode read, whether or not to filter out ribosomal RNA reads, additional alignment parameters, etc. One goal of `slide_snake` is to make the alignment preprocessing modular so that all of these parameters can be compared directly and rigorously. Please note the following:
  - Multiple recipes can be passed for each sample. Include as many as you would like, each separated by a space in the sample sheet.
  - You can also add a new recipe! Just add a new line to `resources/recipe_sheet.csv` and give it a unique name in the 1st (0th for you pythoners) column. 
  - Recipe naming convention:
    - Include the preprocessing steps that you want to use- i.e., including "rRNA.bwa" in the name means that rRNA filtering with bwa alignment will be done prior to genomic alignment. 
<details close>
<summary> SlideSeq/Seeker (Curio) </summary>

Because of read quality issues (indels, low Q scores, etc.) in the SlideSeq barcode read, I have added a few custom strategies for handling these data:
- `seeker_v3.1_hardTrim` - Hard trim the adapter read positions in R1, and use the best barcode correction algorithms in STARsolo
- `seeker_v3.1` - No hard trimming, and use the base positions for barcode/UMI (*Note*, this recipe doesn't work well w/ Curio Seeker b/c of in/del issues w/ the barcode synthesis)
- `seeker_v3.1_MatchLinker` - Match the adapter sequence on R1 (w/ 2 mismatches allowed) and infer barcodes/UMIs from that position (*Note* best performer w/ Curio data)
- `seeker_v3.1_MatchLinker_total` - Same as `seeker_v3.1_noTrimMatchLinker`, but with additional STAR parameters for total RNAseq alignment (more multimappers, looser alignment)
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

## **Barcode handling:**
<details close>
<summary> STAR </summary>
  
- Removed the linker sequence in R1 so that the `1MM_multi` barcode correction in `STARsolo` can be used
- Barcode & UMI paramters for [`STAR`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md):
```
--soloType CB_UMI_Simple	\
--soloUMIstart 14 \
--soloUMIlen 7 \
--soloCBstart 1 \
--soloCBlen 14
```
</details>

<details close>
<summary> kallisto/bustools </summary>
  
- Barcode & UMI paramters for `kallisto bus`
```
-x 0,0,14:0,14,21:1,0,0
```
- Outputs are in

</details>

## **QC:**
- `fastqc` is run on R2 files before (`preTrim`) and after (`postTrim`) adapter trimming
- `qualimap` is used to assess `STAR` alignment

## **Unmapped read analysis:**
- `fastqc` is run on all unmapped reads from `STAR`
- The most abundant unmapped reads are also aligned with `blast`


## **small/micro RNA analysis:**
Work-in-progress...
- Currently have rules set up to align to databases listed below (with `bowtie2`) and generate a count matrix (note, this uses `STAR`'s cell/bead/spot barcoding)
- Need to better optimize alignment params
### miRge3.0
- [Link to documentation](https://mirge3.readthedocs.io/en/latest/quick_start.html)
- [Link to library download](https://sourceforge.net/projects/mirge3/files/miRge3_Lib/)
- Bulk miRg3.0 is implemented, but it is quite slow. I recommend commenting out the target rule(s) unless you have small/total RNA data for which you need this analysis.

### small RNA reference data bases:
- **piRNAs** - [pirbase](http://bigdata.ibp.ac.cn/piRBase/), v3.0
  - *Note* - use "gold standard" piRs b/c there are wayyyyy too many predicted sequences...
- **miRNAs** - [mirbase](https://mirbase.org/), v22
  - Use both hairpin & mature sequences, in that order 
