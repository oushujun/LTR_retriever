#!/usr/bin/perl -w
#use strict;
#Description: This is the core script to calculate the LTR Assembly Index
#Author: Shujun Ou (oushujun@msu.edu)
#Last updated: 09/08/2017

#usage: perl LAI_calc.pl Intact_LTR.bed All_LTR.bed
my $window="1000000"; #1Mb window
my $intact="";
my $all="";
my $genome="";

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-g|genome$/i;
	$intact=$ARGV[$k+1] if /^-i|intact$/i;
	$all=$ARGV[$k+1] if /^-a|allLTR$/i;
	$window=$ARGV[$k+1] if /^-w|window$/i;
	$k++;
	} 

open Genome, "<$genome" or die $!;
open INTACT,"sort -suV $intact |" or die "ERROR: $!";
open ALL,"sort -suV $all |" or die "ERROR: $!";

my $low=1; #lower bound of a window
my $high=$low+$window; #higher bound of a window
my $win_len=0; #object (intact LTR or LTR) length in the window
my $tail=0; #leftover object lenth ahead of a window
my $seq=''; #sequence id
my $len=0; #length of an object
my $total_intact=0; #total length of intact LTR
my $total_LTR=0; #total length of all LTR
my $average_LAI=0; #whole genome average LAI
my $output=''; #to store output info
my %length; #store sequence length
my $genome_len=0; #length of the genome
my $curr_win_len=0; #length of the current window

$/="\n>";
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	my ($chr, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	$seq=length($seq);
	$genome_len+=$seq;
	$length{$chr}=$seq;
	}
$length{"total_genome_length"}=$genome_len;
$/="\n";
close Genome;

while (<INTACT>){
	my ($chr, $from, $to)=(split);
	$seq=$chr if $seq eq '';
	if ($seq ne $chr){ #switch to another chromosome
		$win_len=0 if $win_len==1; #correct the 0bp case
		$win_len+=$tail;
		${$seq}{$low}=$win_len;
		$seq=$chr;
		$low=1;
		$high=$low+$window;
		$win_len=0;
		$tail=0;
		}
	$len=$to-$from+1;
NEXTWIN:
	if ($from>$low and $to<$high){
		$win_len+=$len;
		} 
	elsif ($from>$low and $from<$high and $to>$high) {
		$win_len+=$high-$from+1; #left bound of the right overlapping entry
		$tail=$to-$high+1; #right bound of the right overlapping entry
		}
	elsif ($from<$low and $to>$low and $to<$high){
		$win_len+=$to-$low+1; #right bound of the left overlapping entry
		}
	elsif ($from<$low and $to>$high){
		$win_len+=$window;
		$tail=$to-$high+1; #right bound of the right overlapping entry
		}
	elsif ($from>$high){
		$win_len=0 if $win_len==1; #correct the 0bp case
		${$seq}{$low}=$win_len;
		$win_len=0;
		$low+=$window;
		$high=$low+$window;
		goto NEXTWIN;
		}
	}
$win_len=0 if $win_len==1; #correct the 0bp case
${$seq}{$low}=$win_len; #recycle the second last window
$tail=$window if $tail>$window;
${$seq}{$high}=$tail if $tail>0; #recycle the tail of last window

$low=1;
$high=$low+$window;
$win_len=0;
$tail=0;
$len=0;
$seq='';
print "#Chromosome\tFrom\tTo\tIntact\tTotal\tLAI\n";
while (<ALL>){
	my ($chr, $from, $to)=(split);
	$seq=$chr if $seq eq '';
	if ($seq ne $chr){ #switch to another chromosome
		$win_len=0 if $win_len==1; #correct the 0bp case
		$win_len+=$tail;
		$curr_win_len=$high-$low+1;
		&LAI($low, $high, $seq, $win_len, $curr_win_len) if $win_len>0; #recycle the last entry of a chromosome
		$seq=$chr;
		$low=1;
		$high=$low+$window;
		$high=$length{$seq} if $high > $length{$seq};
		$win_len=0;
		$tail=0;
		}
	$len=$to-$from+1;
NEXTWIN:
	if ($from>$low and $to<$high){
		$win_len+=$len;
		}
	elsif ($from>$low and $from<$high and $to>$high) {
		$win_len+=$high-$from+1; #left bound of the overlapping entry
		$tail=$to-$high+1; #right bound of the overlapping entry
		}
	elsif ($from<$low and $to>$low and $to<$high){
		$win_len+=$to-$low+1; #right bound of the left overlapping entry
		}
	elsif ($from<$low and $to>$high){
		$win_len+=$window;
		$tail=$to-$high+1; #right bound of the right overlapping entry
		}
	elsif ($from>$high){
		$win_len=0 if $win_len==1; #correct the 0bp case
		$curr_win_len=$high-$low+1;
		&LAI($low, $high, $seq, $win_len, $curr_win_len) if $win_len>0;
		$low+=$window;
		$high=$low+$window;
		$high=$length{$seq} if $high > $length{$seq};
		$win_len=0;
		goto NEXTWIN;
		}
	}
$win_len=0 if $win_len==1; #correct the 0bp case
$curr_win_len=$high-$low+1;
&LAI($low, $high, $seq, $win_len, $curr_win_len) if $win_len>0;
$total_intact=sprintf("%.5f", $total_intact/$length{"total_genome_length"});
$total_LTR=sprintf("%.5f", $total_LTR/$length{"total_genome_length"});
if ($total_LTR != 0){
	$average_LAI=sprintf("%.3f", $total_intact/$total_LTR)*100;  #whole genome average LAI
	$average_LAI*=0.1 if $total_LTR<0.01; #scale down to 10% if total LTR content less then 1%
	}
print "whole_genome\tbegin\tend\t$total_intact\t$total_LTR\t$average_LAI\n$output"; #print out all LAI info

sub LAI {
	my ($low, $high, $seq, $win_len, $curr_win_len)=@_[0,1,2,3,4];
	${$seq}{$low}=0 unless exists ${$seq}{$low};
	my $intact=sprintf("%.5f", ${$seq}{$low}/$curr_win_len); #propotion of intact LTR in the window
	my $total=sprintf("%.5f", $win_len/$curr_win_len); #propotion of total LTR in the window
	my $LAI=sprintf("%.3f", (${$seq}{$low}*$curr_win_len)/($win_len*$window))*100; #calculate LTR Assembly Index = intact LTR length / all LTR length, LAI is also weighted by window length 
	$LAI=100.1 if $LAI>100; #LAI>100 could happen in chance if using the reduced redundenty library to find all LTRs
	$LAI*=0.1 if $total<0.01; #scale down to 10% if total LTR content less then 1%
	$total_intact+=${$seq}{$low};
	$total_LTR+=$win_len;
	$output.="$seq\t$low\t$high\t$intact\t$total\t$LAI\n";
	}

