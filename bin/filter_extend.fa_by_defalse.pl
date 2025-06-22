#!/usr/bin/env perl
use strict;
use warnings;

# Shujun Ou (shujun.ou.1@gmail.com) and ChatGPT 4o-mini
# 06/22/2025
# Description: Reads an overlap file (bedtools intersect -a fa.bed -b defalse.bed -wao) and a fasta file, output sequences whose regions
# 		differ from overlap by more than 5bp (including zero overlaps).
# Usage: filter_fasta_by_overlap.pl overlap.bed output.fa > filtered.fa

my ($overlap_file, $fasta_file) = @ARGV;
unless ($overlap_file && $fasta_file) {
    die "Usage: $0 <overlap_file> <fasta_file> > remained.fasta_file\n";
}

# Step 1: parse overlap file, build remove hash
my %remove;
open my $olh, '<', $overlap_file or die "Cannot open overlap file '$overlap_file': $!";
while (<$olh>) {
    chomp;
    next if /^\s*$/;
    my @fields = split /\t|\s+/;
    # fields: chr all_start all_end ... overlap_length is last
    my ($chr, $start, $end, @rest) = @fields;
    my $overlap_len = $rest[-1];
    # normalize coordinates
    ($start, $end) = ($end, $start) if $start > $end;
    my $all_len = $end - $start + 1;
    # record if difference > 5
    if (abs($all_len - $overlap_len) <= 5) {
        my $key = join ":", $chr, $start, $end;
        $remove{$key} = 1;
    }
}
close $olh;

# Step 2: traverse fasta, output sequences in remove
open my $fh, '<', $fasta_file or die "Cannot open fasta file '$fasta_file': $!";

local $/ = "\n>";
while (my $entry = <$fh>) {
    chomp $entry;
    # restore '>' if missing
    $entry = ">" . $entry unless $entry =~ /^>/;
    my ($header, @seq_lines) = split /\n/, $entry;
    my $seq = join('', @seq_lines);
    # extract the "all" region coords from header
    # header format: >chr:start..end|chr:start..end
    my ($chrA, $sA, $eA) = $header =~ /\|(\S+):(\d+)\.\.(\d+)$/;
    next unless defined $chrA;
    # normalize
    ($sA, $eA) = ($eA, $sA) if $sA > $eA;
    my $keyA = join ":", $chrA, $sA, $eA;
    # print if in remove
    print "$header\n$seq\n" unless $remove{$keyA};
}
close $fh;

