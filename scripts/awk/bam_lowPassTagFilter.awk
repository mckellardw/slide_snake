BEGIN {

}

# Skip header lines
/^@/ { print; next }

{
    for (i = 12; i <= NF; i++) {
        if ($i ~ tag":i:") {
            split($i, arr, ":")
            if (arr[3] <= threshold) {
                print
            }
            break
        }
    }
}
