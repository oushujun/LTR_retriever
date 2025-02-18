#!/usr/bin/env bash

# Shujun Ou (shujun.ou.1@gmail.com)
# 02/18/2025

# usage: A wrapper script to calculate solo intact ratio
# sh solo_intact_ratio_wrap.sh TElib.fa genome.RM.out

# Get the directory of the current script
script_dir=$(dirname "$0")

# inputs
lib=$1 # LTR_retriever or EDTA library
RM_out=$2 # RepeatMasker out file

echo "Library: $lib"
echo "RepeatMasker out: $RM_out\n"

echo "1. Analyze library..."
perl "$script_dir/find_LTR.pl" -lib "$lib" > "$lib.info"

echo "2. Find solo LTRs..."
perl "$script_dir/solo_finder.pl" -i "$RM_out" -info "$lib.info" > "$RM_out.solo"

echo "3. Find intact LTRs..."
perl "$script_dir/intact_finder_coarse.pl" "$RM_out" > "$RM_out.intact"

echo "4. Calculate solo:intact ratios..."
perl "$script_dir/solo_intact_ratio.pl" "$RM_out.solo" "$RM_out.intact" > "$RM_out.SI.ratio"

echo "5. Overall solo:intact ratio:"
sh "$script_dir/total_solo_intact_calc.sh" "$RM_out.SI.ratio"

echo "\nDone!"
