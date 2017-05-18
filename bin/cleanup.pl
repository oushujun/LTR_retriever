#!usr/bin/perl -w
use strict;

my $usage="
        Usage: perl cleanup.pl -f sample.fa -e [int] -r [0, 1] > sample.cln.fa 
        \n";
my $version="
cleanup.pl
cleanup: clean up tandem repeat sequence and sequencing gap sequence in the fasta file
Author: Shujun Ou (oushujun\@msu.edu), Department of Horticulture, Michigan State University, East Lansing, MI, 48823, USA
Version: 	1.6 Add alignment score control (-e [int]) and missing count control (-c [int])
		1.5 Add missing rate control (-r [0, 1]) and modify missing count control (-nc not control)
		1.0 2014/06/02
\n";

my $target="n";
my $func_nc=1; #1, do $n_count screening
my $n_count=0; #count the $target in each sequence, if it exceeds this number, it will be discarted.
my $n_rate=1; #count the $target in each sequence, if it exceeds this percentage, it will be discarted.
my $align_score=1000; #-e para, dft:1000 
my $max_seed=2000; #maximum period size to report
my $cleanN=0; #1 will remove $target="n" in output sequence
my $trf=1; #1 will enable tandem repeat finder (default), 0 will not
my $file;

my $k=0;
foreach (@ARGV){
	$target=$ARGV[$k+1] if /^-misschar$/i;
	$func_nc=0 if /^-Nscreen$/i;
	$n_count=$ARGV[$k+1] if /^-nc$/i;
	$n_rate=$ARGV[$k+1] if /^-nr$/i;
	$align_score=$ARGV[$k+1] if /^-minscore$/i;
	$file=$ARGV[$k+1] if /^-f$/i;
	$cleanN=1 if /^-cleanN$/i;
	$trf=0 if /^-trf$/i;
	$k++;
	}

my $workdir=`dirname "$0"`;
$workdir=~s/\n$//;

my %tandem;
my $tandem='';
$tandem=`$workdir/trf409.legacylinux64 $file 2 7 7 80 10 $align_score $max_seed -ngs -h -l 6` if $trf==1;
while ($tandem=~s/\@(.*)\n?//){
	$tandem{$1}=$1;
	}

open Info, ">$file.cleanup";
open File, "<$file" or die $!;
$/="\n>";
while (<File>){
	if (/^$/){next}
	s/>//g;
	my ($id, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $length=length $seq;
	my $mark=0;
	my $count=0;
	$count++ while $seq=~/$target/gi;
	my $count_rate=$count/$length;
	if ($count_rate>=$n_rate){
		print Info "$id\t$count_rate missing\n";
		$mark=1;
		}
	elsif ($count>=$n_count && $func_nc==1 && $n_count>0){
		print Info "$id\t$count sequence gap\n";
		$mark=1;
		}
	if (exists $tandem{$id}){
		print Info "$id\ttandem sequence\n";
		$mark=1;
		}
	unless ($mark==1){
		$seq=~s/$target//gi if $cleanN==1;
		print ">$id\n$seq\n";
		}
	}
$/="\n";
close Info;
close File;
