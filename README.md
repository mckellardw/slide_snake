# `curio_STARsolo`
Preprocessing, alignment, QC, and quantification workflow for Curio Bioscience data (SlideSeq/Seeker)
**David W. McKellar**


## Dependencies:
- `cutadapt` [v3.4](https://cutadapt.readthedocs.io/en/stable/)
- `fastqc` [v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v2.2.2a](http://qualimap.conesalab.org/)
- `vsearch` [v2.17.0_linux_x86_64](https://github.com/torognes/vsearch)
- `BLAST`

## Alignment:
- After adapter trimming, reads are aligned with `STARsolo` and `kallisto`/`bustools` to generate count matrices
- Outputs are in `SAMPLE_ID/STARsolo/Solo.out` & `SAMPLE_ID/kb/counts_unfiltered`

## Barcode handling:
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
SAMPLE_ID/
├── cutadapt_polyA_report.txt
├── cutadapt_polyG_report.txt
├── postTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_final_fastqc.html
│   └── SAMPLE_ID_R2_final_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_fastqc.html
│   └── SAMPLE_ID_R2_fastqc.zip
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
