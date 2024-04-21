
# Notes for rRNA:


# VASAseq alignment code:
```
# mapping short reads
${p2bwa}/bwa aln ${ref} ${folder}/${fq} > ${folder}/aln_${fq%.f*q}.sai 
${p2bwa}/bwa samse ${ref} ${folder}/aln_${fq%.f*q}.sai ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.aln-ribo.bam &

# mapping normal reads
${p2bwa}/bwa mem -t 8 -h 15 ${ref} ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.mem-ribo.bam & 
```
- Split into short/less short reads?
- Use bowtie instead?