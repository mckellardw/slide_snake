# remove_ansi.awk
## Remove ANSI escape sequences from input text

# Usage:
## cat input.txt | awk -f remove_ansi.awk > clean_output.txt

{
    gsub(/\x1B\[[0-9;]*[mK]/, "")  # Remove ANSI escape sequences
    print                          # Print the cleaned line
}