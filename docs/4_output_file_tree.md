# **Output file tree:**
*Note*- whereever you see `{***RECIPE_NAME***}` written, these outputs are provided for every recipe requiested for that sample.
```
{SAMPLE_ID}/
├── bb
│   ├── whitelist_1.txt
│   ├── whitelist_2.txt
│   ├── whitelist_adapter.txt
│   └── whitelist.txt
├── cutadapt2.json
├── cutadapt.json
├── cutadapt.log
├── cutadapt_round2.log
├── fastqc
│   ├── postCutadapt_R1
│   │   ├── cut_R1_fastqc.html
│   │   └── cut_R1_fastqc.zip
│   ├── postCutadapt_R2
│   │   ├── cut_R2_fastqc.html
│   │   └── cut_R2_fastqc.zip
│   ├── preCutadapt_R1
│   │   ├── merged_R1_fastqc.html
│   │   └── merged_R1_fastqc.zip
│   ├── preCutadapt_R2
│   │   ├── merged_R2_fastqc.html
│   │   └── merged_R2_fastqc.zip
│   ├── rRNA_bwa_R1
│   │   ├── final_filtered_R1_fastqc.html
│   │   └── final_filtered_R1_fastqc.zip
│   ├── rRNA_bwa_R2
│   │   ├── final_filtered_R2_fastqc.html
│   │   └── final_filtered_R2_fastqc.zip
│   ├── rRNA_STAR_R1
│   │   ├── final_filtered_R1_fastqc.html
│   │   └── final_filtered_R1_fastqc.zip
│   ├── rRNA_STAR_R2
│   │   ├── final_filtered_R2_fastqc.html
│   │   └── final_filtered_R2_fastqc.zip
│   ├── twiceCutadapt_R1
│   │   ├── twiceCut_R1_fastqc.html
│   │   └── twiceCut_R1_fastqc.zip
│   ├── twiceCutadapt_R2
│   │   ├── twiceCut_R2_fastqc.html
│   │   └── twiceCut_R2_fastqc.zip
│   └── unmapped
│       └── {***RECIPE_NAME***}
│           ├── Unmapped.out.mate1_fastqc.html
│           ├── Unmapped.out.mate1_fastqc.zip
│           ├── Unmapped.out.mate2_fastqc.html
│           └── Unmapped.out.mate2_fastqc.zip
├── kb
│   └── {***RECIPE_NAME***}
│       ├── inspect.corrected.bus.json
│       ├── kallisto_align.log
│       ├── output.sorted.bus
│       ├── raw
│       │   ├── bustools_count.log
│       │   ├── output.barcodes.txt.gz
│       │   ├── output.genes.txt.gz
│       │   ├── output.h5ad
│       │   └── output.mtx.gz
│       └── run_info.json
├── kbpython
│   └── {***RECIPE_NAME***}
│       ├── counts_unfiltered
│       │   ├── cells_x_genes.barcodes.txt.gz
│       │   ├── cells_x_genes.genes.names.txt
│       │   ├── cells_x_genes.genes.txt.gz
│       │   ├── cells_x_genes.mtx.gz
│       │   └── output.h5ad
│       ├── inspect.json
│       ├── kb_info.json
│       ├── kbpython_standard.log
│       ├── run_info.json
│       └── transcripts.txt
├── miRge_bulk
│   └── {***RECIPE_NAME***}
│   │   ├── a2IEditing.detail.txt
│   │   ├── a2IEditing.report.csv
│   │   ├── a2IEditing.report.newform.csv
│   │   ├── annotation.report.csv
│   │   ├── annotation.report.html
│   │   ├── index_data.js
│   │   ├── mapped.csv
│   │   ├── miR.Counts.csv
│   │   ├── miRge3_visualization.html
│   │   ├── miR.RPM.csv
│   │   ├── run.log
│   │   ├── sample_miRge3.gff
│   │   ├── twiceCut_R2_novel_miRNAs
│   │   │   └── twiceCut_R2_novel_miRNAs_report.csv
│   │   ├── unmapped.csv
│   │   └── unmapped.log
├── ont
│   ├── adapter_scan.tsv
│   ├── merged.fq.gz
│   ├── merged.log
│   ├── merged_stranded.fq.gz
│   ├── minimap2.log
│   ├── sorted.bam
│   ├── sorted.bam.bai
│   ├── tmp
│   │   ├── 20231213_1452_3B_PAM09861_415471fe.fq
│   │   └── 20231215_0943_3B_PAM09861_de42ca77.fq
│   └── tmp.sam
├── qualimap
│   └── rRNA
│       └── bwa
│           ├── css
│           │   ├── agogo.css
│           │   ├── ajax-loader.gif
│           │   ├── basic.css
│           │   ├── bgfooter.png
│           │   ├── bgtop.png
│           │   ├── comment-bright.png
│           │   ├── comment-close.png
│           │   ├── comment.png
│           │   ├── doctools.js
│           │   ├── down.png
│           │   ├── down-pressed.png
│           │   ├── file.png
│           │   ├── jquery.js
│           │   ├── minus.png
│           │   ├── plus.png
│           │   ├── pygments.css
│           │   ├── qualimap_logo_small.png
│           │   ├── report.css
│           │   ├── searchtools.js
│           │   ├── underscore.js
│           │   ├── up.png
│           │   ├── up-pressed.png
│           │   └── websupport.js
│           ├── images_qualimapReport
│           │   ├── Coverage\ Profile\ Along\ Genes\ (High).png
│           │   ├── Coverage\ Profile\ Along\ Genes\ (Low).png
│           │   ├── Coverage\ Profile\ Along\ Genes\ (Total).png
│           │   ├── Reads\ Genomic\ Origin.png
│           │   └── Transcript\ coverage\ histogram.png
│           ├── qualimapReport.html
│           ├── raw_data_qualimapReport
│           │   ├── coverage_profile_along_genes_(high).txt
│           │   ├── coverage_profile_along_genes_(low).txt
│           │   └── coverage_profile_along_genes_(total).txt
│           ├── rnaseq_qc.log
│           ├── rnaseq_qc_results.csv
│           └── rnaseq_qc_results.txt
├── R1_trimming.log
├── rRNA
│   ├── bwa
│   │   ├── aligned.bam
│   │   ├── aligned_sorted.bam
│   │   ├── aligned_sorted.bam.bai
│   │   ├── bwa_mem.log
│   │   ├── final_filtered_R1.fq.gz
│   │   ├── final_filtered_R2.fq.gz
│   │   └── rRNA_readID.list
│   └── STARsolo
│       ├── Aligned.sortedByCoord.out.bam
│       ├── final_filtered_R1.fq.gz
│       ├── final_filtered_R2.fq.gz
│       ├── Log.final.out
│       ├── Log.out
│       ├── Log.progress.out
│       ├── SJ.out.tab
│       └── Solo.out
│           ├── Barcodes.stats
│           └── GeneFull
│               ├── Features.stats
│               ├── filtered
│               │   ├── barcodes.tsv
│               │   ├── features.tsv
│               │   └── matrix.mtx
│               ├── raw
│               │   ├── barcodes.tsv
│               │   ├── features.tsv
│               │   ├── matrix.mtx
│               │   └── UniqueAndMult-EM.mtx
│               ├── Summary.csv
│               └── UMIperCellSorted.txt
├── STARsolo
│   └── {***RECIPE_NAME***}
│       ├── Aligned.sortedByCoord.out.bam
│       ├── Log.final.out
│       ├── Log.out
│       ├── Log.progress.out
│       ├── SJ.out.tab
│       ├── Solo.out
│       │   ├── Barcodes.stats
│       │   ├── Gene
│       │   │   ├── Features.stats
│       │   │   ├── filtered
│       │   │   │   ├── barcodes.tsv.gz
│       │   │   │   ├── features.tsv.gz
│       │   │   │   └── matrix.mtx.gz
│       │   │   ├── raw
│       │   │   │   ├── barcodes.tsv.gz
│       │   │   │   ├── features.tsv.gz
│       │   │   │   ├── matrix.h5ad
│       │   │   │   ├── matrix.mtx.gz
│       │   │   │   ├── UniqueAndMult-EM.h5ad
│       │   │   │   ├── UniqueAndMultEM.h5ad
│       │   │   │   └── UniqueAndMult-EM.mtx.gz
│       │   │   ├── Summary.csv
│       │   │   └── UMIperCellSorted.txt
│       │   ├── GeneFull
│       │   │   ├── Features.stats
│       │   │   ├── filtered
│       │   │   │   ├── barcodes.tsv.gz
│       │   │   │   ├── features.tsv.gz
│       │   │   │   └── matrix.mtx.gz
│       │   │   ├── raw
│       │   │   │   ├── barcodes.tsv.gz
│       │   │   │   ├── features.tsv.gz
│       │   │   │   ├── matrix.h5ad
│       │   │   │   ├── matrix.mtx.gz
│       │   │   │   ├── UniqueAndMult-EM.h5ad
│       │   │   │   ├── UniqueAndMultEM.h5ad
│       │   │   │   └── UniqueAndMult-EM.mtx.gz
│       │   │   ├── Summary.csv
│       │   │   └── UMIperCellSorted.txt
│       │   └── Velocyto
│       │       ├── Features.stats
│       │       ├── filtered
│       │       │   ├── ambiguous.mtx.gz
│       │       │   ├── barcodes.tsv.gz
│       │       │   ├── features.tsv.gz
│       │       │   ├── spliced.mtx.gz
│       │       │   └── unspliced.mtx.gz
│       │       ├── raw
│       │       │   ├── ambiguous.mtx.gz
│       │       │   ├── barcodes.tsv.gz
│       │       │   ├── features.tsv.gz
│       │       │   ├── spliced.mtx.gz
│       │       │   └── unspliced.mtx.gz
│       │       └── Summary.csv
│       ├── Unmapped.out.mate1.fastq.gz
│       └── Unmapped.out.mate2.fastq.gz
└── tmp
    ├── cut_R1.fq
    ├── cut_R1.fq.gz
    ├── cut_R2.fastq.gz
    ├── cut_R2.fq.gz
    ├── merged_R1.fq.gz
    ├── merged_R2.fq.gz
    ├── merged_trimmed_R1.fq.gz
    ├── twiceCut_R1.fq.gz
    ├── twiceCut_R2.fastq.gz
    └── twiceCut_R2.fq.gz
```