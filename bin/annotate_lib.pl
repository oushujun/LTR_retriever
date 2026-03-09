#!/usr/bin/env perl

##To annotate library with family and direction information
##Usage: perl annotate_lib.pl $index.scn.adj $index.LTRlib.fa
##Shujun Ou (oushujun@msu.edu) 06/23/2016


use warnings;
use strict;

open List, "<$ARGV[0]" or die "ERROR: $!";
open Lib, "<$ARGV[1]" or die "ERROR: $!";

my %info;
while (<List>){
	next if /^#/;
	s/^\s+//;
	my ($start, $end, $lstart, $lend, $rstart, $rend, $direction, $fam)=(split)[0,1,3,4,6,7,12,18];
	my ($istart, $iend)=($lend+1, $rstart-1);
	next unless defined $fam;
	$info{"$start..$end"}=[$direction, $fam];
	$info{"$lstart..$lend"}=[$direction, $fam];
	$info{"$rstart..$rend"}=[$direction, $fam];
	$info{"$istart..$iend"}=[$direction, $fam];
	}

$/="\n>";
while (<Lib>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$id=~s/\s+//g;
	$seq=~s/\s+//g;
	my ($chr, $region)=("NA", "NA");
	($chr=$1, $id=$2, $region=$3) if $id=~/^(.*):([0-9]+\.\.[0-9]+)\|(.*)$/;
	$region="LTR" if $region=~/LTR_[1-2]/;
	$region="INT" if $region=~/LTR_IN/;
	my $is_whole = ($region=~/^[0-9]+\.\.[0-9]+$/) ? 1 : 0; # whole element: region is coordinates
	my ($direction, $fam)=("?", "unknown");
	# Normalize coordinate order for lookup (headers may have descending coords for minus-strand)
	my $lookup_id = $id;
	if ($lookup_id =~ /^(\d+)\.\.(\d+)$/ && $1 > $2) {
		$lookup_id = "$2..$1";
	}
	($direction, $fam)=@{$info{$lookup_id}}[0,1] if defined $info{$lookup_id};
	if ($direction eq "-"){
		$seq=~tr/ATGCatgc/TACGtacg/;
		$seq=reverse $seq;
		}

	if ($is_whole){
		print ">$chr:${id}#LTR/$fam\n$seq\n";
	} else {
		print ">$chr:${id}_$region#LTR/$fam\n$seq\n";
	}
	}
close Lib;
close List;
