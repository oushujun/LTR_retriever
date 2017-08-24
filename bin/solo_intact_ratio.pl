#!/usr/bin/perl -w
use strict;

#compute solo-intact LTR ratio for given lists
#Shujun Ou (oushujun@msu.edu)
#07-25-2017

my $usage="perl this_script.pl solo_count intact_count > solo_intact_ratio";

open Solo, "<$ARGV[0]" or die $!;
open Intact, "<$ARGV[1]" or die $!;

my %all;

while (<Intact>){
	chomp;
	my ($id, $intact)=(split);
	$all{$id}=["0", $intact];
	}
close Intact;

while (<Solo>){
	chomp;
	my ($id, $solo)=(split);
	if (exists $all{$id}){
		$all{$id}[0]=$solo;
		} else {
		$all{$id}=[$solo, "0"];
		}
	}
foreach my $id (sort {$a cmp $b} keys %all){
	my $ratio="NA";
	my ($solo, $intact)=($all{$id}[0], $all{$id}[1]);
	if ($intact==0){
		$ratio="inf";
		} 
	elsif ($solo==0){
		$ratio=0;
		}
	else {
		$ratio=sprintf ("%.1f", $solo/$intact);
		}
	print "$id\t$solo\t$intact\t$ratio\n";
	}
close Solo;

