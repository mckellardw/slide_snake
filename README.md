# `slide_snake`
Preprocessing, alignment, QC, and quantification workflow for SlideSeq data  

The goal of this project is to build a workflow for assessing different alignment parameterizations, and digging into artifacts. ***Contributions are welcome!*** See the companion workflow for 10x Genomics' chemistries [here](https://github.com/mckellardw/txg_snake).  

## Dependencies:
- `cutadapt` [v3.4](https://cutadapt.readthedocs.io/en/stable/)
- `fastqc` [v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v2.2.2a](http://qualimap.conesalab.org/)
- `vsearch` [v2.17.0_linux_x86_64](https://github.com/torognes/vsearch)
- `BLAST`


## Trimming
Five Prime adapter(s):
- Reverse complement of the SlideSeq TSO to reduce strand invasion artifacts [**CCCATTCACTCTGCGTTGATACCAGCTT**]

Three Prime adapter(s):
- "A" homopolymers [**100-A**]
- "G" homopolymers, important for Illumina 2-color sequencers [**100-G**]
- "T" homopolymers [**100-T**]
- Nextera adapter sequence [**CTGTCTCTTATA**]
- Reverse complement of Nextera sequence [**TATAAGAGACAG**] 
- SlideSeq template switch olgo (TSO) - remove any polyadenylated TSOs [**AAGCTGGTATCAACGCAGAGTGAATGGG**]
- Illumina unversal sequence [**AGATCGGAAGAG**]


## Alignment:
- After adapter trimming, reads are aligned with `STARsolo` and `kallisto`/`bustools` to generate count matrices
- Outputs are in `{SAMPLE_ID}/STARsolo/Solo.out` & `{SAMPLE_ID}/kb/counts_unfiltered`
- Different recipes are written out in `resources/chemistry_sheet.csv`, and must be specified for each sample within the sample sheet

### Chemistry/recipe descriptions:
- `seeker_v3.1` - Hard trim the adapter read positions in R1, and use the best barcode correction algorithms in STARsolo
- `seeker_v3.1_noTrim` - No hard trimming, and use the base positions for barcode/UMI (*Note*, this recipe doesn't work well w/ Curio Seeker b/c of in/del issues w/ the barcode synthesis)
- `seeker_v3.1_noTrimMatchLinker` - Match the adapter sequence on R1 (w/ 2 mismatches allowed) and infer barcodes/UMIs from that position (*Note* best performer w/ Curio data)
- `seeker_v3.1_noTrim_total` - Same as `seeker_v3.1_noTrimMatchLinker`, but with additional STAR parameters for total RNAseq alignment (more multimappers, looser alignment)


### Barcode handling:
#### STAR
- Removed the linker sequence in R1 so that the `1MM_multi` barcode correction in `STARsolo` can be used
- Barcode & UMI paramters for [`STAR`](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md):
```
--soloType CB_UMI_Simple	\
--soloUMIstart 14 \
--soloUMIlen 7 \
--soloCBstart 1 \
--soloCBlen 14
```

#### kallisto/bustools
- Barcode & UMI paramters for `kallisto bus`
```
-x 0,0,14:0,14,21:1,0,0
```
- Outputs are in

## QC:
- `fastqc` is run on R2 files before (`preTrim`) and after (`postTrim`) adapter trimming
- `qualimap` is used to assess `STAR` alignment

## Unmapped read analysis
- `fastqc` is run on all unmapped reads from `STAR`
- The most abundant unmapped reads are also aligned with `blast`

## Output file tree:
```
{SAMPLE_ID}/
├── log.cutadapt.json
├── kb
│   ├── counts_unfiltered
│   │   ├── output.barcodes.txt
│   │   ├── output.genes.txt
│   │   └── output.mtx
│   ├── inspect.corrected.bus.json
│   ├── kallisto_align.log
│   ├── matrix.ec
│   ├── output.bus
│   ├── output.corrected.bus
│   ├── output.sorted.bus
│   ├── run_info.json
│   └── transcripts.txt
├── postTrim_fastqc_R1_out
│   ├── {SAMPLE_ID}_R1_adapterTrim_fastqc.html
│   └── {SAMPLE_ID}_R1_adapterTrim_fastqc.zip
├── postTrim_fastqc_R2_out
│   ├── {SAMPLE_ID}_R2_final_fastqc.html
│   └── {SAMPLE_ID}_R2_final_fastqc.zip
├── preTrim_fastqc_R1_out
│   ├── {SAMPLE_ID}_R1_fastqc.html
│   └── {SAMPLE_ID}_R1_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── {SAMPLE_ID}_R2_fastqc.html
│   └── {SAMPLE_ID}_R2_fastqc.zip
├── qualimap
│   ├── css
│   │   ├── agogo.css
│   │   ├── ajax-loader.gif
│   │   ├── basic.css
│   │   ├── bgfooter.png
│   │   ├── bgtop.png
│   │   ├── comment-bright.png
│   │   ├── comment-close.png
│   │   ├── comment.png
│   │   ├── doctools.js
│   │   ├── down.png
│   │   ├── down-pressed.png
│   │   ├── file.png
│   │   ├── jquery.js
│   │   ├── minus.png
│   │   ├── plus.png
│   │   ├── pygments.css
│   │   ├── qualimap_logo_small.png
│   │   ├── report.css
│   │   ├── searchtools.js
│   │   ├── underscore.js
│   │   ├── up.png
│   │   ├── up-pressed.png
│   │   └── websupport.js
│   ├── images_qualimapReport
│   │   ├── Coverage Profile Along Genes (High).png
│   │   ├── Coverage Profile Along Genes (Low).png
│   │   ├── Coverage Profile Along Genes (Total).png
│   │   ├── Junction Analysis.png
│   │   ├── Reads Genomic Origin.png
│   │   └── Transcript coverage histogram.png
│   ├── qualimapReport.html
│   ├── raw_data_qualimapReport
│   │   ├── coverage_profile_along_genes_(high).txt
│   │   ├── coverage_profile_along_genes_(low).txt
│   │   └── coverage_profile_along_genes_(total).txt
│   └── rnaseq_qc_results.txt

├── STARsolo
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   ├── Solo.out
│   │   ├── Barcodes.stats
│   │   ├── Gene
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   ├── GeneFull
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   └── Velocyto
│   │       ├── Features.stats
│   │       ├── filtered
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       ├── raw
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       └── Summary.csv
│   ├── Unmapped.out.mate1.fastq.gz
│   └── Unmapped.out.mate2.fastq.gz
├── tmp
│   ├── whitelist_1.txt
│   ├── whitelist_2.txt
│   └── whitelist.txt
├── Unmapped_fastqc_out
│   ├── Unmapped.out.mate1_fastqc.html
│   ├── Unmapped.out.mate1_fastqc.zip
│   ├── Unmapped.out.mate2_fastqc.html
│   └── Unmapped.out.mate2_fastqc.zip
└── Unmapped.out.mate2_blastResults.txt
```

## Helpful links
- [Barcode download from Curio](http://3.222.73.200/)
