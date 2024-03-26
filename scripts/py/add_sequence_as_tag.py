import argparse
import pysam
import gzip

# Usage:
## python scripts/py/add_sequence_as_tag.py /gpfs/commons/groups/innovation/dwm/ontxg/data/slide_out/SHVN1/ont/minimap2/sorted.bam /gpfs/commons/groups/innovation/dwm/ontxg/data/slide_out/SHVN1/tmp/ont/cut_R1.fq.gz shvn1.bam 

def add_sequence_as_tag(bam_file, fastq_file, output_bam_file, sequence_tag='BC', quality_tag='BY'):
    bam = pysam.AlignmentFile(bam_file, "rb")
    new_bam = pysam.AlignmentFile(output_bam_file, "wb", template=bam)

    # Check if the FASTQ file is gzipped
    if fastq_file.endswith('.gz'):
        # Open the gzipped FASTQ file
        with gzip.open(fastq_file, "rt") as fastq:
            fastq_lines = fastq.readlines()
    else:
        # Open the FASTQ file normally
        with open(fastq_file, "r") as fastq:
            fastq_lines = fastq.readlines()

    for i in range(0, len(fastq_lines), 4):
        read_id = fastq_lines[i].strip().split(' ')[0][1:]
        sequence = fastq_lines[i+1].strip()
        quality_scores = fastq_lines[i+3].strip()
        
        for read in bam.fetch():
            if read.query_name == read_id:
                read.set_tag(sequence_tag, sequence, 'Z')
                read.set_tag(quality_tag, quality_scores, 'Z')
                new_bam.write(read)

    bam.close()
    new_bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add sequence from FASTQ as tag in BAM file.')
    parser.add_argument('bam_file', help='Path to the input BAM file')
    parser.add_argument('fastq_file', help='Path to the FASTQ file containing barcode information')
    parser.add_argument('output_bam_file', help='Path to the output BAM file')
    parser.add_argument('--sequence_tag', default='BC', help='Tag name for the sequence (default: BC)')
    parser.add_argument('--quality_tag', default='BY', help='Tag name for the quality scores (default: BY)')

    args = parser.parse_args()

    add_sequence_as_tag(args.bam_file, args.fastq_file, args.output_bam_file, args.sequence_tag, args.quality_tag)
