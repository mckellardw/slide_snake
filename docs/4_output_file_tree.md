# **Output file tree:**
*Note*- wherever you see `{***RECIPE_NAME***}` written, these outputs are provided for every recipe requested for that sample.
```
{SAMPLE_ID}/
├── bc
│   ├── map.txt
│   ├── whitelist_1.txt
│   ├── whitelist_2.txt
│   ├── whitelist_adapter.txt
│   └── whitelist.txt
├── short_read
│   ├── cutadapt2.json
│   ├── cutadapt.json
│   ├── cutadapt.log
│   ├── cutadapt_round2.log
│   ├── fastqc
│   │   ├── postCutadapt_R1
│   │   │   ├── cut_R1_fastqc.html
│   │   │   └── cut_R1_fastqc.zip
│   │   ├── postCutadapt_R2
│   │   │   ├── cut_R2_fastqc.html
│   │   │   └── cut_R2_fastqc.zip
│   │   ├── preCutadapt_R1
│   │   │   ├── merged_R1_fastqc.html
│   │   │   └── merged_R1_fastqc.zip
│   │   ├── preCutadapt_R2
│   │   │   ├── merged_R2_fastqc.html
│   │   │   └── merged_R2_fastqc.zip
│   │   ├── rRNA_bwa_R1
│   │   │   ├── final_filtered_R1_fastqc.html
│   │   │   └── final_filtered_R1_fastqc.zip
│   │   ├── rRNA_bwa_R2
│   │   │   ├── final_filtered_R2_fastqc.html
│   │   │   └── final_filtered_R2_fastqc.zip
│   │   ├── twiceCutadapt_R1
│   │   │   ├── twiceCut_R1_fastqc.html
│   │   │   └── twiceCut_R1_fastqc.zip
│   │   ├── twiceCutadapt_R2
│   │   │   ├── twiceCut_R2_fastqc.html
│   │   │   └── twiceCut_R2_fastqc.zip
│   │   └── unmapped
│   │       └── {***RECIPE_NAME***}
│   │           ├── Unmapped.out.mate1_fastqc.html
│   │           ├── Unmapped.out.mate1_fastqc.zip
│   │           ├── Unmapped.out.mate2_fastqc.html
│   │           └── Unmapped.out.mate2_fastqc.zip
│   ├── kb
│   │   └── {***RECIPE_NAME***}
│   │       ├── inspect.corrected.bus.json
│   │       ├── kallisto_align.log
│   │       ├── output.sorted.bus
│   │       ├── raw
│   │       │   ├── bustools_count.log
│   │       │   ├── output.barcodes.txt.gz
│   │       │   ├── output.genes.txt.gz
│   │       │   ├── output.h5ad
│   │       │   └── output.mtx.gz
│   │       └── run_info.json
│   ├── kbpython
│   │   └── {***RECIPE_NAME***}
│   │       ├── counts_unfiltered
│   │       │   ├── cells_x_genes.barcodes.txt.gz
│   │       │   ├── cells_x_genes.genes.names.txt
│   │       │   ├── cells_x_genes.genes.txt.gz
│   │       │   ├── cells_x_genes.mtx.gz
│   │       │   └── output.h5ad
│   │       ├── inspect.json
│   │       ├── kb_info.json
│   │       ├── kbpython_standard.log
│   │       ├── run_info.json
│   │       └── transcripts.txt
│   ├── miRge_bulk
│   │   └── {***RECIPE_NAME***}
│   │       ├── a2IEditing.detail.txt
│   │       ├── a2IEditing.report.csv
│   │       ├── annotation.report.csv
│   │       ├── annotation.report.html
│   │       ├── mapped.csv
│   │       ├── miR.Counts.csv
│   │       ├── miRge3_visualization.html
│   │       ├── miR.RPM.csv
│   │       ├── run.log
│   │       ├── unmapped.csv
│   │       └── unmapped.log
│   ├── qualimap
│   │   └── rRNA
│   │       └── bwa
│   │           ├── qualimapReport.html
│   │           ├── rnaseq_qc.log
│   │           ├── rnaseq_qc_results.csv
│   │           └── rnaseq_qc_results.txt
│   ├── R1_trimming.log
│   ├── rRNA
│   │   ├── bwa
│   │   │   ├── aligned.bam
│   │   │   ├── aligned_sorted.bam
│   │   │   ├── aligned_sorted.bam.bai
│   │   │   ├── bwa_mem.log
│   │   │   ├── no_rRNA_R1.fq.gz
│   │   │   ├── no_rRNA_R2.fq.gz
│   │   │   └── rRNA_readID.list
│   │   └── ribodetector
│   │       ├── no_rRNA_R1.fq.gz
│   │       └── no_rRNA_R2.fq.gz
│   ├── STARsolo
│   │   └── {***RECIPE_NAME***}
│   │       ├── Aligned.sortedByCoord.out.bam
│   │       ├── Aligned.sortedByCoord.out.bam.bai
│   │       ├── Aligned.sortedByCoord.out.dedup.bam
│   │       ├── Log.final.out
│   │       ├── Log.out
│   │       ├── Log.progress.out
│   │       ├── SJ.out.tab
│   │       ├── Solo.out
│   │       │   ├── Barcodes.stats
│   │       │   ├── Gene
│   │       │   │   ├── Features.stats
│   │       │   │   ├── filtered
│   │       │   │   │   ├── barcodes.tsv.gz
│   │       │   │   │   ├── features.tsv.gz
│   │       │   │   │   └── matrix.mtx.gz
│   │       │   │   ├── raw
│   │       │   │   │   ├── barcodes.tsv.gz
│   │       │   │   │   ├── features.tsv.gz
│   │       │   │   │   ├── matrix.h5ad
│   │       │   │   │   ├── matrix.mtx.gz
│   │       │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │       │   │   ├── Summary.csv
│   │       │   │   └── UMIperCellSorted.txt
│   │       │   ├── GeneFull
│   │       │   │   ├── Features.stats
│   │       │   │   ├── filtered
│   │       │   │   │   ├── barcodes.tsv.gz
│   │       │   │   │   ├── features.tsv.gz
│   │       │   │   │   └── matrix.mtx.gz
│   │       │   │   ├── raw
│   │       │   │   │   ├── barcodes.tsv.gz
│   │       │   │   │   ├── features.tsv.gz
│   │       │   │   │   ├── matrix.h5ad
│   │       │   │   │   ├── matrix.mtx.gz
│   │       │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │       │   │   ├── Summary.csv
│   │       │   │   └── UMIperCellSorted.txt
│   │       │   └── Velocyto
│   │       │       ├── Features.stats
│   │       │       ├── filtered
│   │       │       │   ├── ambiguous.mtx.gz
│   │       │       │   ├── barcodes.tsv.gz
│   │       │       │   ├── features.tsv.gz
│   │       │       │   ├── spliced.mtx.gz
│   │       │       │   └── unspliced.mtx.gz
│   │       │       ├── raw
│   │       │       │   ├── ambiguous.mtx.gz
│   │       │       │   ├── barcodes.tsv.gz
│   │       │       │   ├── features.tsv.gz
│   │       │       │   ├── spliced.mtx.gz
│   │       │       │   └── unspliced.mtx.gz
│   │       │       └── Summary.csv
│   │       ├── Unmapped.out.mate1.fastq.gz
│   │       └── Unmapped.out.mate2.fastq.gz
│   └── tmp
│       ├── cut_R1.fq.gz
│       ├── cut_R2.fq.gz
│       ├── merged_R1.fq.gz
│       ├── merged_R2.fq.gz
│       ├── twiceCut_R1.fq.gz
│       └── twiceCut_R2.fq.gz
├── ont
│   ├── adapter_scan.tsv
│   ├── merged.fq.gz
│   ├── merged.log
│   ├── merged_stranded.fq.gz
│   ├── minimap2
│   │   └── {***RECIPE_NAME***}
│   │       ├── sorted.bam
│   │       └── sorted.bam.bai
│   └── tmp
│       └── ...
└── logs
    └── ...
```
