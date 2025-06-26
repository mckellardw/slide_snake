# fq_hardTrim.awk
## Hard trim FASTQ reads to specified positions

# Usage:
## awk -v s=start -v S=startTrim -v E=endTrim -f fq_hardTrim.awk input.fastq > trimmed.fastq

# Process FASTQ records
(FNR-1) % 2 == 0 { 
    name=$1
    tag=$2
    next 
}

(FNR-2) % 4 == 0 { 
    substr($0,1,s) substr($0,S) 
}

{ 
    print name tag
    print substr($0,1,s) substr($0,S,E) 
}
