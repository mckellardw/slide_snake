BEGIN {
    quality = ARGV[1]
    ARGV[1] = ""
}

{
    if ($5 < quality) {
        print
    }
}
