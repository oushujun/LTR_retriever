#!/usr/bin/env perl
use strict;
use warnings;

# Shujun Ou (shujun.ou.1@gmail.com) and ChatGPT o4-mini
# 06/22/2025
# Description: filter out LTRharvest scn file with LTR_retriever defalse file. output scn lines not exist in defalse

# Usage: filter_scn.pl defalse_file scn_file > remained_scn_file

my ($def_f, $scn_f) = @ARGV;
unless ($def_f && $scn_f) {
    die "\nUsage: $0 defalse_file scn_file > remained_scn_file\n\n";
}

# Read defalse file into a hash of keys "seqid:start:end"
open my $DEF, '<', $def_f or die "Cannot open '$def_f': $!\n";
my %seen;
while (<$DEF>) {
    chomp;
    next unless /^\S+:[0-9]+..[0-9]+/;
    # First column is like NC_053489:39826..43100
    my $col1 = (split)[0];
    my ($seq, $range) = split /:/, $col1, 2;
    my ($a, $b) = split /\.\./, $range, 2;
    # Normalize so smaller first
    ($a, $b) = ($b, $a) if $a > $b;
    $seen{"$seq:$a:$b"} = 1;
}
close $DEF;

# Process scn file
open my $SCN, '<', $scn_f or die "Cannot open '$scn_f': $!\n";
while (<$SCN>) {
    if (/^\#/) {
        # Print any header/comment lines
        print;
        next;
    }
    chomp;
    next if /^\s*$/;
    my @F = split;
    # first two cols are start and end, last col is seqid
    my ($s, $e) = @F[0,1];
    my $seq = $F[-1];
    # Normalize
    ($s, $e) = ($e, $s) if $s > $e;
    my $key = "$seq:$s:$e";
    # If not seen in defalse, print it
    print "$_\n" unless $seen{$key};
}
close $SCN;

