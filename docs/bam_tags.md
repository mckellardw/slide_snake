# .bam tags

Doing our best to be consistent with the .bam tags used by cellranger & STAR...
From [Dave Tang's excellent blog](https://davetang.org/muse/2018/06/06/10x-single-cell-bam-files/):
| Line | Tag | Description |
|------|-----|-------------|
| 1    | NA  | Query name |
| 2    | NA  | Flag |
| 3    | NA  | Reference name |
| 4    | NA  | Position |
| 5    | NA  | Mapping quality |
| 6    | NA  | Cigar string |
| 7    | NA  | Reference name of mate |
| 8    | NA  | Position of mate |
| 9    | NA  | Template length |
| 10   | NA  | Sequence |
| 11   | NA  | Sequence quality |
| 12   | NH  | Number of reported alignments for query |
| 13   | HI  | Query hit index |
| 14   | AS  | Alignment score |
| 15   | nM  | Number of mismatches per pair |
| 16   | RE  | Region type (E = exonic, N = intronic, I = intergenic) |
| 17   | BC  | Sample index read |
| 18   | QT  | Sample index read quality |
| 19   | CR  | Cell barcode |
| 20   | CY  | Cell barcode read quality |
| 21   | CB  | Cell barcode that is error-corrected and confirmed against a list of known-good barcode sequences |
| 22   | UR  | Unique Molecular Identifier (UMI) |
| 23   | UY  | UMI read quality |
| 24   | UB  | UMI that is error-corrected among other molecular barcodes with the same cellular barcode and gene alignment |
| 25   | RG  | Read group |

Other tags:

| ##   | XB  | CB + UMI (combined tag for simple deduplication) |