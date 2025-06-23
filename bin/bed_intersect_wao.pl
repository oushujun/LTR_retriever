#!/usr/bin/env perl
use strict;
use warnings;

# check args
die "Usage: $0 <A.bed> <B.bed> > result\n" unless @ARGV == 2;
my ($fileA, $fileB) = @ARGV;

# peek first non-track line in B to count columns
open my $peekB, '<', $fileB or die "open $fileB: $!";
my $firstB;
while (defined($firstB = <$peekB>)) {
    next if $firstB =~ /^track/;
    chomp $firstB;
    last;
}
close $peekB;
my $nBcols = scalar split /\t/, $firstB;

# open sorted streams
open my $A, "-|", "sort -k1,1 -k2,2n $fileA" or die "sort A: $!";
open my $B, "-|", "sort -k1,1 -k2,2n $fileB" or die "sort B: $!";

# in‐memory buffer of “active” B intervals: [chr, start, end, \@fields]
my @buffer;
my $lineB = <$B>;
chomp $lineB if defined $lineB;

while (defined(my $lineA = <$A>)) {
    chomp $lineA;
    next if $lineA =~ /^track/;
    my @a = split /\t/, $lineA;
    my ($chrA, $stA, $enA) = @a[0,1,2];

    # 1) purge any buffered B intervals that are on the wrong chr or end ≤ A.start
    @buffer = grep { $_->[0] eq $chrA && $_->[2] > $stA } @buffer;

    # 2) read ahead in B until we pass A.end or hit a later chr
    while (defined $lineB) {
        my @b = split /\t/, $lineB;
        my ($chrB, $stB, $enB) = @b[0,1,2];

        last if $chrB gt $chrA
             || ($chrB eq $chrA && $stB >= $enA);

        # only buffer intervals on same chr that could overlap
        if ($chrB eq $chrA && $enB > $stA) {
            push @buffer, [ $chrB, $stB, $enB, \@b ];
        }

        $lineB = <$B>;
        chomp $lineB if defined $lineB;
    }

    # 3) scan buffer for real overlaps
    my $found = 0;
    for my $ent (@buffer) {
        # each $ent = [ chr, start, end, \@b_fields ]
        next unless $ent->[0] eq $chrA;
        my ($stB, $enB, $bf) = @{$ent}[1,2,3];

        # half-open overlap
        my $ov_s = $stA > $stB ? $stA : $stB;
        my $ov_e = $enA < $enB ? $enA : $enB;
        next if $ov_e <= $ov_s;

        $found = 1;
        my $ov_len = $ov_e - $ov_s;
        print join("\t", @a, @$bf, $ov_len), "\n";
    }

    # 4) if no overlaps, emit the “wao” dummy line
    unless ($found) {
        print join("\t",
            @a,
            ('.') x $nBcols,
            0,
        ), "\n";
    }
}

close $A;
close $B;

