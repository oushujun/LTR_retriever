#!/usr/bin/perl -w
#use strict;
#Description: This is the core script to calculate the LTR Assembly Index
#Author: Shujun Ou (oushujun@msu.edu)
#Last updated: 09/08/2017

#usage: perl LAI_calc.pl Intact_LTR.bed All_LTR.bed
my $window="2000000"; #1Mb window
my $average_age="0.8"; #mean age (MY) of intact LTR-RTs; default=0.8
my $intact="";
my $all="";
my $genome="";

my $k=0;
foreach (@ARGV){
	$genome=$ARGV[$k+1] if /^-g|genome$/i;
	$intact=$ARGV[$k+1] if /^-i|intact$/i;
	$all=$ARGV[$k+1] if /^-a|allLTR$/i;
	$window=$ARGV[$k+1] if /^-w|window$/i;
	$average_age=$ARGV[$k+1] if /^-t|time$/i;
	$k++;
	}

open Genome, "<$genome" or die $!;
open INTACT,"sort -suV -k1,3 $intact |" or die "ERROR: $!";
open ALL,"sort -suV -k1,3 $all |" or die "ERROR: $!";

my $low=1; #lower bound of a window
my $high=$low+$window; #higher bound of a window
my $win_len=0; #object (intact LTR or LTR) length in the window
my $tail=0; #leftover object lenth ahead of a window
my $seq=''; #sequence id
my $len=0; #length of an object
my $total_intact=0; #total length of intact LTR
my $total_LTR=0; #total length of all LTR
my $average_LAI=0; #whole genome average LAI
my $average_LAI_adj=0; #whole genome average LAI adjusted by mean age of all intact LTR-RTs
my $output=''; #to store output info
my %length; #store sequence length
my $genome_len=0; #length of the genome
my $curr_win_len=0; #length of the current window
my $age_tot=0; #total age of intact LTR-RTs in the window
#my $age_avg=0; #average age of intact LTR-RTs in the window
my $num_intact=0; #number of intact LTR-RTs in the window

$/="\n>";
while (<Genome>){
	next if /^>\s?$/;
	chomp;
	s/>//g;
	my ($chr, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	$chr=~s/\s+$//; #remove space at the end of the seq ID
	$seq=length($seq);
	$genome_len+=$seq;
	$length{$chr}=$seq;
	}
$length{"total_genome_length"}=$genome_len;
$/="\n";
close Genome;

while (<INTACT>){
	my ($chr, $from, $to, $age)=(split);
	$seq=$chr if $seq eq '';
	if ($seq ne $chr){ #switch to another chromosome
		$win_len=0 if $win_len==1; #correct the 0bp case
		$win_len+=$tail;
		${$seq}{$low}=[$win_len, $age_tot/$num_intact];
		$seq=$chr;
		$low=1;
		$high=$low+$window;
		$win_len=0;
		$tail=0;
		$age_tot=0;
		$num_intact=0;
		}
	$len=$to-$from+1;
NEXTWIN:
	if ($from>$low and $to<$high){
		$win_len+=$len;
		$age_tot+=$age;
		$num_intact++;
		} 
	elsif ($from>$low and $from<$high and $to>$high) {
		$win_len+=$high-$from+1; #left bound of the right overlapping entry
		$proportion=($high-$from+1)/($to-$from+1);
		$age_tot+=$age*$proportion; #the proportional age of the LTR-RT falls into the current window is counted
		$num_intact+=$proportion; #the proportion of LTR-RT falls into the current window is counted
		$tail=$to-$high+1; #right bound of the right overlapping entry
		}
	elsif ($from<$low and $to>$low and $to<$high){
		$win_len+=$to-$low+1; #right bound of the left overlapping entry
		$proportion=($to-$low+1)/($to-$from+1);
		$age_tot+=$age*$proportion; #the proportional age of the LTR-RT falls into the current window is counted
		$num_intact+=$proportion; #the proportion of LTR-RT falls into the current window is counted
		}
	elsif ($from<$low and $to>$high){
		$win_len+=$window;
		$proportion=$window/($to-$from+1);
		$age_tot+=$age*$proportion; #the proportional age of the LTR-RT falls into the current window is counted
		$num_intact+=$proportion; #the proportion of LTR-RT falls into the current window is counted
		$tail=$to-$high+1; #right bound of the right overlapping entry
		}
	elsif ($from>$high){
		$win_len=0 if $win_len==1; #correct the 0bp case
#		${$seq}{$low}=$win_len;
		if ($num_intact!=0){
			${$seq}{$low}=[$win_len, $age_tot/$num_intact];
			} else {
			${$seq}{$low}=[$win_len, 4];
			}
		$win_len=0;
		$age_tot=0;
		$num_intact=0;
		$low+=$window;
		$high=$low+$window;
		goto NEXTWIN;
		}
	}
$win_len=0 if $win_len==1; #correct the 0bp case
${$seq}{$low}=[$win_len, $age_tot/$num_intact]; #recycle the second last window
$tail=$window if $tail>$window;
${$seq}{$high}=$tail if $tail>0; #recycle the tail of last window

$low=1;
$high=$low+$window;
$win_len=0;
$tail=0;
$len=0;
$seq='';
print "#Chromosome\tFrom\tTo\tIntact\tTotal\tLAI\tLAI_adj\tIntact_LTR_Age(MY)\n";
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
	$average_LAI_adj=sprintf("%.1f", $average_LAI-16.48*(0.8-$average_age)); #adjust the raw LAI score by mean age of LTR-RT
	$average_LAI_adj=0 if $average_LAI_adj<0;
	$average_age=sprintf("%.1f", $average_age);
	}
print "whole_genome\tbegin\tend\t$total_intact\t$total_LTR\t$average_LAI\t$average_LAI_adj\t$average_age\n$output"; #print out all LAI info

sub LAI {
	my ($low, $high, $seq, $win_len, $curr_win_len)=@_[0,1,2,3,4];
	${$seq}{$low}=[0,0] unless exists ${$seq}{$low}[0];
	my $intact=sprintf("%.5f", ${$seq}{$low}[0]/$curr_win_len); #propotion of intact LTR in the window
	my $total=sprintf("%.5f", $win_len/$curr_win_len); #propotion of total LTR in the window
	my $age=${$seq}{$low}[1]/1000000; #mean age (MY) of LTR-RTs in the window
	my $LAI=sprintf("%.3f", (${$seq}{$low}[0]*$curr_win_len)/($win_len*$window))*100; #calculate LTR Assembly Index = intact LTR length / all LTR length, LAI is also weighted by window length 
	$LAI=100.1 if $LAI>100; #LAI>100 could happen in chance if using the reduced redundenty library to find all LTRs
	$LAI*=0.1 if $total<0.01; #scale down to 10% if total LTR content less than 1%
	my $LAI_adj=sprintf("%.1f", $LAI-16.48*(0.8-$age)); #adjust the raw LAI score by mean age of LTR-RT
	$LAI_adj=0 if $LAI_adj<0;
	$age=sprintf("%.1f", $age);
	$total_intact+=${$seq}{$low}[0];
	$total_LTR+=$win_len;
	$output.="$seq\t$low\t$high\t$intact\t$total\t$LAI\t$LAI_adj\t$age\n";
	}

