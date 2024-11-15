# Sample sheets in `slide_snake`

In this pipeline, sample sheets are used to keep things organized, and to add flexibility to the alignment and preprocessing options. Here are a few important notes about how to prepare your sample sheet:

## Sample Sheet Columns

1. **sampleID**: Unique identifier for each sample. Used for labeling and tracking samples throughout the workflow.
2. **fastq_R1**: Path to the Read 1 FASTQ file(s). Multiple files can be separated by spaces, or represented with regular expressions.
3. **fastq_R2**: Path to the Read 2 FASTQ file(s). Multiple files can be separated by spaces, or represented with regular expressions.
4. **recipe**: Recipe(s) to use for each sample. Multiple recipes can be separated by spaces. Specifies the preprocessing and alignment strategies to be applied.
5. **recipe_ONT**: Recipe(s) for ONT data. Multiple recipes can be separated by spaces. Specifies the preprocessing and alignment strategies for Oxford Nanopore Technologies data.
6. **whitelist**: Path to the barcode whitelist file. Used for demultiplexing and identifying valid barcodes.
7. **species**: Species name (e.g., human, mouse). Used to select the appropriate reference genome and annotation files.
8. **STAR_ref**: Path to the STAR reference genome directory. Used for alignment with the STAR aligner.
10. **kb_idx**: Path to the kallisto index file. Used for pseudo-alignment with kallisto.
11. **kb_t2g**: Path to the kallisto transcript-to-gene mapping file. Used for converting transcript-level counts to gene-level counts.
12. **rRNA_ref**: Path to the rRNA reference file for BWA. Used for aligning and quantifying rRNA reads.
13. **rRNA_gtf**: Path to the GTF file for rRNA annotations. Used for annotating rRNA reads.
14. **genome_fa**: Path to the genome FASTA file. Used for various steps requiring the reference genome sequence.
9. **genes_gtf**: Path to the GTF file for gene annotations. Used for annotating aligned reads and generating gene counts.
15. **cDNA_fa**: Path to the transcripts/cDNA FASTA file. Used for long-read transcriptome alignment w/ minimap2.


## Example Sample Sheet

Here is an example of a sample sheet in CSV format:

```csv
sampleID,fastq_R1,fastq_R2,recipe,recipe_ONT,whitelist,species,STAR_ref,genes_gtf,kb_idx,kb_t2g,rRNA_ref,rRNA_gtf,genome_fa
test_Human_Lymphoma_DLBCL_Region2_10um,data/test_pathodbit/SAMN43151516_10k_R2.fq.gz,data/test_pathodbit/SAMN43151516_10k_R1.fq.gz,dbit-pretrim,,/gpfs/commons/groups/innovation/dwm/slide_snake/resources/dbit_whitelist/Spatial_barcode_100x100.txt,human,/gpfs/commons/groups/innovation/dwm/genomes/homo_sapiens/STAR/GRCh38_GENCODEv47,/gpfs/commons/groups/innovation/dwm/genomes/homo_sapiens/gencode/release_v47/gencode.v47.annotation.gtf,/gpfs/commons/groups/innovation/dwm/genomes/homo_sapiens/kallisto/GRCh38_GENCODEv47/transcriptome.idx,/gpfs/commons/groups/innovation/dwm/genomes/homo_sapiens/kallisto/GRCh38_GENCODEv47/t2g.txt,/gpfs/commons/groups/innovation/dwm/ref_snake/out/homo_sapiens/rRNA/bwa_mem2/ref.fa.gz,/gpfs/commons/groups/innovation/dwm/ref_snake/out/homo_sapiens/rRNA/raw/annotations.gtf,/gpfs/commons/groups/innovation/dwm/genomes/homo_sapiens/gencode/release_v47/GRCh38.primary_assembly.genome.fa
```

## Notes
- Ensure that all file paths are correct and accessible.
- Multiple recipes can be passed for each sample, separated by spaces. This is useful for benchmarking different preprocessing and alignment strategies.
- The sample sheet should be saved as a CSV file and the path to this file should be specified in the config.yaml file under SAMPLE_SHEET_PATH.