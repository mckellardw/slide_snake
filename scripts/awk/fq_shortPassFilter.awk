# fq_shortPassFilter.awk
##  Filter out reads that are *shorter* than the specified length (`maxLength`)

# Usage:
## awk -v maxLength=650 -f fq_shortPassFilter.awk input.fastq > shortReads.fastq

#TODO - fix this so that the extra "+" and read IDs are not also written...
{
    pos = NR % 4;
    if (pos == 2) {
        if (length($0) >= maxLength) {
            print $0;
        }
    } else if (pos == 4) {
        if (length($0) >= maxLength) {
            print $0;
        }
    } else{
        print $0;
    }
}

# BEGIN {
#     RS="@"; # Set the record separator to "@"
#     # maxLength = 100; # Define your maxLength here
# }

# {
#     if (length($2) <= maxLength) { 
#         print "@"$1;    # Print the read identifier
#         print $2;       # Print the sequence
#         print $3;       # Print the separator line
#         print $4;       # Print the quality scores
#     }
# }


