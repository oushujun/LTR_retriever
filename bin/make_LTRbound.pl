#!/usr/bin/env perl

## Generate LTR boundary map file from a whole-element LTR library and scn.adj file
## Usage: perl make_LTRbound.pl $scn_adj $LTRlib.fa > $LTRlib.fa.LTRbound
## Output format: element_name\ttotal_len\tlLTR_len\trLTR_len
## Shujun Ou (oushujun@msu.edu) 03-06-2026

use warnings;
use strict;

die "Usage: perl make_LTRbound.pl <scn.adj> <LTRlib.fa>\n" unless @ARGV == 2;

open SCN, "<$ARGV[0]" or die "ERROR: Cannot open $ARGV[0]: $!";
open LIB, "<$ARGV[1]" or die "ERROR: Cannot open $ARGV[1]: $!";

# Read scn.adj to get LTR boundary info keyed by element_start..element_end
my %bound;
while (<SCN>){
	next if /^#/;
	s/^\s+//;
	my @f = split /\s+/;
	next unless @f >= 9;
	my ($estart, $eend, $lstart, $lend, $llen, $rstart, $rend, $rlen) = @f[0,1,3,4,5,6,7,8];
	$bound{"$estart..$eend"} = [$llen, $rlen];
}
close SCN;

# Read library FASTA headers and look up boundary info
while (<LIB>){
	next unless /^>/;
	chomp;
	s/^>//;
	my $name = (split /\s+/)[0];
	# Extract element coordinates: chr:start..end#class or chr:start..end_REGION#class
	my ($coords) = $name =~ /:([0-9]+\.\.[0-9]+)/;
	next unless defined $coords;
	# Normalize coordinates to ascending order for lookup
	my $lookup = $coords;
	if ($lookup =~ /^(\d+)\.\.(\d+)$/ && $1 > $2) {
		$lookup = "$2..$1";
	}
	if (defined $bound{$lookup}){
		my ($llen, $rlen) = @{$bound{$lookup}};
		# Get total sequence length from the coordinates
		my ($s, $e) = split /\.\./, $coords;
		my $total_len = abs($e - $s) + 1;
		print "$name\t$total_len\t$llen\t$rlen\n";
	}
}
close LIB;
