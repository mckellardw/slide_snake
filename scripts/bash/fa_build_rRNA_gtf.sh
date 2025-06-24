#!/bin/bash

# Check if the correct number of arguments was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: scripts/bash/extract_rRNA_gtf.sh FASTA_cDNA GTF_rRNA"
    exit 1
fi

# Assign the first argument to FASTA_cDNA and the second to GTF_rRNA
FASTA_cDNA=$1
GTF_rRNA=$2

# Manually build gtf with this hideous awk code...
zcat ${FASTA_cDNA} \
| awk \
    -v RS=">" \
    -v FS=" " \
    -vOFS='' \
    '$5=="gene_biotype:rRNA" || $5=="gene_biotype:Mt_rRNA" {
        seqname = $1;
        source = "custom";
        feature = "exon";
        start = "1";

        n=split($0, lines, "\n");
        len = 0;
        for (i = 2; i <= n; ++i) {
            len += length(lines[i]);
        }
        end = len;

        score = ".";
        strand = "+";
        frame = ".";
        
        split($4, b, ":");
        gene_id = b[2];
        split($6, c, ":");
        transcript_type = c[2];
        split($7, d, ":");
        transcript_name = d[2];
        split($7, e, ":");
        gene_name = e[2];
        split($5, f, ":");
        gene_type = f[2];

        attributes = "gene_id \"" gene_id "\"; transcript_id \"" seqname "\"; gene_type \"" gene_type "\"; gene_name \"" gene_name "\"; transcript_type \"" transcript_type "\"; transcript_name \"" transcript_name "\";";

        print seqname,"\t", source,"\t", feature,"\t", start,"\t", end,"\t", score,"\t", strand,"\t", frame,"\t", attributes;
    }' \
| gzip \
> ${GTF_rRNA}