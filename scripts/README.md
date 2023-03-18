# scripts

- `cache_h5ad.py` - load in .mtx file and save a .h5ad for easy scanpy analysis later
- `hardTrimfq.awk` - hard trimming on the adapter sequence in R1
- `internal_adapter_trim_R1.py` - find and remove the adapter sequence (between the two bead barcodes, BB_1 & BB_2) in read 1
- `kb.sh` - run kallisto/bustools
- `split_dedup.sh` - parallelized umi_tools deduplication (split reads across chromosomes, then dedup!). Also filters the .bam file based on mapq scores
- `splitNfqs.py` - split a fastq file into N chunks (parallelized, uses an easy `sed` command)