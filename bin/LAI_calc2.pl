#!/usr/bin/perl -w
use strict;

#Description: This is the core script to calculate the LTR Assembly Index
#Author: Shujun Ou (oushujun@msu.edu)
#Last updated: 02/05/2018
#Note: Memory consumption of this scrip is approx. 2X the size of the input genome

#usage: perl LAI_calc.pl -i Intact_LTR.bed -a All_LTR.bed -g genome.fa [options]
my $window="3000000"; #3Mb/window
my $step="300000"; #300Kb/step, not fully implemented function
my $iden="94"; #mean identity (%) of LTR sequences in the haploid genome; default=94 (means no adjustment for age)
my $iden_slope="3.13310889698647"; #the correction slope for LTR sequence identity (age)
my $intact="";
my $all="";
my $genome="";

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-g|genome$/i;
	$intact=$ARGV[$k+1] if /^-i|intact$/i;
	$all=$ARGV[$k+1] if /^-a|allLTR$/i;
	$window=$ARGV[$k+1] if /^-w|window$/i;
	$step=$ARGV[$k+1] if /^-s|step$/i;
	$iden=$ARGV[$k+1] if /^-d$/i;
	$iden_slope=$ARGV[$k+1] if /^-k$/i;
	$k++;
	}

open INTACT,"sort -suV -k1,3 $intact |" or die "ERROR: $!";
open ALL,"sort -suV -k1,3 $all |" or die "ERROR: $!";
open Genome, "<$genome" or die $!;

my %length; #store sequence length
my %intact; #store intact LTR-RT info
my %total; #store all LTR sequence info
my @seqID; #store chr ID names in order
my $genome_len=0; #length of the genome
my $output=''; #stores output info

$/="\n>";
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	my ($chr, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	$chr=~s/\s+$//; #remove space at the end of the seq ID
	push @seqID, $chr;
	$seq=length($seq);
	$genome_len+=$seq;
	$length{$chr}=$seq;
	$intact{$chr}='0' x $seq; #create a '0' string that has the same length as the chr
	$total{$chr}='0' x $seq; #create another '0' string that has the same length as the chr
	}
$/="\n";
close Genome;

while (<INTACT>){
	my ($chr, $from, $to)=(split)[0,1,2];
	my $len=$to-$from+1;
	substr($intact{$chr}, $from-1, $len)="i" x $len; #substitute '0' with 'i' where intact LTR-RT is occurred
	}
close INTACT;

while (<ALL>){
	my ($chr, $from, $to)=(split);
	my $len=$to-$from+1;
	substr($total{$chr}, $from-1, $len)="a" x $len; #substitute '0' with 'a' where LTR sequence is occurred
	}
close ALL;

my ($tot_int_count, $tot_all_count, $tot_int_per, $tot_all_per, $tot_LAI, $tot_LAI_adj)=(0, 0, 0, 0, 0, 0);
foreach my $chr (@seqID){
	my $int_count = $intact{$chr} =~ tr/i/i/;
	my $all_count = $total{$chr} =~ tr/a/a/;
	$tot_int_count += $int_count;
	$tot_all_count += $all_count;
#estimate LAI based on windows and steps
	my $win_len = $window;
	for (my $start=1; $win_len == $window; $start += $step){
		my $end = $start+$window-1;
		$end = $length{$chr} if $end > $length{$chr};
		$win_len = $end-$start+1; #the actual size of the window (chromosome end may have win_len < window)
		my $win_seq_int = substr($intact{$chr}, $start, $win_len); #intact LTR-RT information in the window
		my $win_seq_all = substr($total{$chr}, $start, $win_len); #total LTR sequence information in the win
		my $win_int_count = $win_seq_int =~ tr/i/i/; #intact LTR-RT length
		my $win_all_count = $win_seq_all =~ tr/a/a/; #total LTR sequence length
		my $win_int_per = sprintf("%.4f", $win_int_count/$win_len); #propotion of intact LTR-RT in the window
		my $win_all_per = sprintf("%.4f", $win_all_count/$win_len); #propotion of total LTR sequence in the window
		my $win_LAI = 0;
		$win_LAI = sprintf("%.2f", ($win_int_count*$win_len*100)/($win_all_count*$window)) if $win_all_count != 0; #calculate LTR Assembly Index = intact LTR length / all LTR length, LAI is also weighted by window length
		$win_LAI = 100.1 if $win_LAI > 100; #LAI>100 could happen in chance if using the reduced redundenty library to find all LTRs
		$win_LAI *= 0.1 if $win_all_per < 0.01; #scale down to 10% if total LTR content less than 1%
		my $win_LAI_adj = sprintf("%.2f", $win_LAI + $iden_slope * (94 - $iden));
		$output .= "$chr\t$start\t$end\t$win_int_per\t$win_all_per\t$win_LAI\t$win_LAI_adj\n";
		}
	}
#estimate genome-wide LAI
$tot_int_per = sprintf("%.4f", $tot_int_count/$genome_len);
$tot_all_per = sprintf("%.4f", $tot_all_count/$genome_len);
$tot_LAI = sprintf("%.2f", ($tot_int_count*100)/$tot_all_count);
$tot_LAI = 100.1 if $tot_LAI > 100;
$tot_LAI *= 0.1 if $tot_all_per < 0.01;
$tot_LAI_adj = sprintf("%.2f", $tot_LAI + $iden_slope * (94 - $iden));
print "#Chromosome\tFrom\tTo\tIntact\tTotal\tLAI\tLAI_adj\n";
print "whole_genome\tbegin\tend\t$tot_int_per\t$tot_all_per\t$tot_LAI\t$tot_LAI_adj\n$output"; #print out all LAI info

