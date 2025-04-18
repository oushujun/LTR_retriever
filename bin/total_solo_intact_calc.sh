#!/bin/bash

# Input file
input_file=$1

# Use awk to sum columns 2 and 3, then calculate the ratio
awk -v input_file="$input_file" '
NR > 1 {  # Skip the header row
    sum2 += $2
    sum3 += $3
}
END {
    if (sum3 == 0) {
        ratio = "inf"
    } else {
        ratio = sum2 / sum3
    }
    print "Genome\tTotal_solo\tTotal_intact\tOverall_SI_ratio"
    printf "%s\t%d\t%d\t%s\n", input_file, sum2, sum3, ratio
}
' "$input_file"
