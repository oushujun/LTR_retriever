#!/usr/bin/bash perl
# Created with the help by DeepSeek v3.0
# Shujun Ou (shujun.ou.1@gmail.com)
# 01/05/2024
# Usage: generate $index.retriever.scn.adj from the $index.defalse file.
# perl $script_path/bin/defalse2scn_adj.pl $index.defalse > $index.retriever.scn.adj

use strict;
use warnings;

open my $defalse_fh, "<$ARGV[0]" or die "Cannot open input";

while (<$defalse_fh>) {
    chomp;
    next unless /motif/;  # Skip lines that do not start with "Chr"

    # Split the line into fields
    my @fields = split /\t/;

    # Extract relevant information
    my ($chr_range, $status, $motif, $tsd, $tsd_start_end, $ltr_start_end, $in_range, $score, $strand, $type, $ltr, $id) = @fields;

    # Parse the chromosome range
    my ($chr, $start, $end) = $chr_range =~ /(\S+):(\d+)\.\.(\d+)/;

    # Ensure start is less than end
    if ($start > $end) {
        ($start, $end) = ($end, $start);  # Swap start and end if necessary
    }

    # Parse the TSD start and end positions (if available)
    my ($tsd_start, $tsd_end) = ($tsd_start_end =~ /(\d+)\.\.(\d+)/) ? ($1, $2) : ("NA", "NA");

    # Parse the LTR start and end positions (if available)
    my ($ltr_start, $ltr_end) = ($ltr_start_end =~ /(\d+)\.\.(\d+)/) ? ($1, $2) : ("NA", "NA");

    # Parse the IN range
    my ($in_start, $in_end) = $in_range =~ /IN:(\d+)\.\.(\d+)/;

    # Ensure in_start is less than in_end
    if ($in_start > $in_end) {
        ($in_start, $in_end) = ($in_end, $in_start);  # Swap if necessary
    }

    # Calculate lengths
    my $total_length = $end - $start + 1;
    my $left_flank_length = $in_start - $start;
    my $right_flank_length = $end - $in_end;

    # Ensure lengths are non-negative
    $left_flank_length = 0 if $left_flank_length < 0;
    $right_flank_length = 0 if $right_flank_length < 0;

    # Prepare the output line
    my $output_line = join("\t",
        $start, $end, $total_length,
        $start, $in_start - 1, $left_flank_length,
        $in_end + 1, $end, $right_flank_length,
        $score, 0, $chr, $strand,
        $tsd =~ s/^TSD://r, "$tsd_start..$tsd_end", "$ltr_start..$ltr_end",
        $motif =~ s/motif://r, $ltr, $type, $id
    );

    # Write to the output file
    print "$output_line\n";
}

close $defalse_fh;
