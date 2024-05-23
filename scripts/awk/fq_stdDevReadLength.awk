BEGIN {
    sum = 0;
    sum_sq = 0;
    count = 0;
    mean = 0;
}

{
    if (NR % 4 == 2) {
        len = length($0);
        sum += len;
        sum_sq += len * len;
        count++;
    }
}

END {
    if (count > 0) {
        mean = sum / count;
        variance = sum_sq / count - mean * mean;
        std_dev = sqrt(variance);
        print std_dev;
    } else {
        print "No reads found.";
    }
}
