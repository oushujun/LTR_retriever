#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl bed_intersect_wao_fixed.pl A.bed B.bed

my $fileA = shift or die "Usage: $0 <A.bed> <B.bed>\n";
my $fileB = shift or die "Usage: $0 <A.bed> <B.bed>\n";

# Constants
my $LOOKBACK = 20;  # number of previous B lines to keep

# Open sorted streams
open my $fhA, "-|", "sort -k1,1 -k2,2n $fileA" or die "Cannot open sorted stream for $fileA: $!";
open my $fhB, "-|", "sort -k1,1 -k2,2n $fileB" or die "Cannot open sorted stream for $fileB: $!";

# Read first lines
my $lineA = <$fhA>;
my $lineB = <$fhB>;

# Circular buffer for last $LOOKBACK B intervals
my @b_buffer;

while (defined $lineA) {
    chomp($lineA);
    next if $lineA =~ /^track/;
    my @a = split /\t/, $lineA;
    my ($chromA, $startA, $endA) = @a[0,1,2];

    my $found_overlap = 0;

    # Save current position in B file
    my $pos_before_scan = tell($fhB);

    # Try matching A with buffered B lines
    foreach my $i (0 .. $#b_buffer) {
        my ($startB, $endB, $lineB) = @{$b_buffer[$i]};
        my $overlap_start = $startA > $startB ? $startA : $startB;
        my $overlap_end   = $endA < $endB   ? $endA   : $endB;
        my $overlap_len   = $overlap_end > $overlap_start ? $overlap_end - $overlap_start : 0;

        if ($overlap_len > 0) {
            my @b_fields = split /\t/, $lineB;
            print join("\t", @a, @b_fields, $overlap_len), "\n";
            $found_overlap = 1;
        }
    }

    # Now scan forward through B
    while (defined $lineB) {
        chomp($lineB);
        next if $lineB =~ /^track/;
        my @b = split /\t/, $lineB;
        my ($chromB, $startB, $endB) = @b[0,1,2];

        # Chromosome mismatch
        if ($chromA lt $chromB) {
            last;
        } elsif ($chromA gt $chromB) {
            $lineB = <$fhB>;
            next;
        }

        # No overlap on start/end
        if ($endA <= $startB) {
            last;
        } elsif ($endB <= $startA) {
            save_to_buffer();
            $lineB = <$fhB>;
            next;
        }

        # Overlap found!
        my $overlap_start = $startA > $startB ? $startA : $startB;
        my $overlap_end   = $endA < $endB   ? $endA   : $endB;
        my $overlap_len   = $overlap_end - $overlap_start;

        print join("\t", @a, @b, $overlap_len), "\n";
        $found_overlap = 1;

        save_to_buffer();
        if ($endB < $endA) {
            $lineB = <$fhB>;
        } else {
            last;
        }
    }

    # If no match, rewind B and check again
    unless ($found_overlap) {
        seek($fhB, $pos_before_scan, 0);
        $lineB = <$fhB>;

        my @backup_lines;
        for (1..$LOOKBACK) {
            last unless defined $lineB;
            push @backup_lines, $lineB;
            $lineB = <$fhB>;
        }

        seek($fhB, $pos_before_scan, 0);
        $lineB = <$fhB>;

        while (defined $lineB) {
            chomp($lineB);
            next if $lineB =~ /^track/;
            my @b = split /\t/, $lineB;
            my ($chromB, $startB, $endB) = @b[0,1,2];

            if ($chromA eq $chromB) {
                my $overlap_start = $startA > $startB ? $startA : $startB;
                my $overlap_end   = $endA < $endB   ? $endA   : $endB;
                my $overlap_len   = $overlap_end > $overlap_start ? $overlap_end - $overlap_start : 0;

                if ($overlap_len > 0) {
                    my @b_fields = split /\t/, $lineB;
                    print join("\t", @a, @b_fields, $overlap_len), "\n";
                    $found_overlap = 1;
                }
            }

            $lineB = <$fhB>;
        }

        # Restore file pointer after rewind
        foreach my $l (@backup_lines) {
            $lineB = $l;
            last;
        }
    }

    # Output dummy line if no overlap
    unless ($found_overlap) {
        print join("\t", @a, ('.') x 3, (-1) x 3, '.', 0), "\n";
    }

    $lineA = <$fhA>;
}

close $fhA;
close $fhB;

# Subroutines

sub save_to_buffer {
    my @b = split /\t/, $lineB;
    my ($startB, $endB) = @b[1,2];
    push @b_buffer, [ $startB, $endB, $lineB ];
    shift @b_buffer if @b_buffer > $LOOKBACK;
}
