# trim N bases from either side of a read sequence in a .bam file

# BEGIN {
#     if (ARGC != 3) {
#         print "Usage: cat <input.sam> | awk -f bam_trimNBases.awk <N> <trim_option>"
#         print "option: start, end, both"
#         exit 1
#     }

#     N = ARGV[1]
#     option = ARGV[2]
#     ARGV[1] = ARGV[2] = ""
# }

BEGIN { 
    FS = OFS = "\t" 
}

{
    if ($0 !~ /^@/) {
        if (option == "start") {
            $10 = substr($10, N + 1)
            $11 = substr($11, N + 1)
        } else if (option == "end") {
            $10 = substr($10, 1, length($10) - N)
            $11 = substr($11, 1, length($11) - N)
        } else if (option == "both") {
            $10 = substr($10, N + 1, length($10) - 2*N)
            $11 = substr($11, N + 1, length($11) - 2*N)
        }
    }
    print $0
}
