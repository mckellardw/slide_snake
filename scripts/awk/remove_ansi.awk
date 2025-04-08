# remove_ansi.awk
{
    gsub(/\x1B\[[0-9;]*[mK]/, ""); # Remove ANSI escape sequences
    print; # Print the cleaned line
}