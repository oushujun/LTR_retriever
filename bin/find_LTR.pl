#!/usr/bin/env perl
use strict;
use warnings;

my $usage = "#######################

Report LTR regions of TE libraries.

Usage:	perl find_LTR.pl -lib [file]

Credit: Shujun Ou (shujun.ou.1\@gmail.com) 12/10/2020

#######################
";
# Known bug: if the LTR region cannot align as a whole (eg, break into two pieces), this script will fail

my $blastplus = '';
my $min_align = 100; #minimal alignment length (bp)
my $max_off = 20; #maximum terminal offshift for LTR-LTR alignments
my $lib = "NA"; # library sequences

my $k = 0;
foreach (@ARGV){
	$min_align = $ARGV[$k+1] if /^-min_aln$/i;
	$max_off = $ARGV[$k+1] if /^-max_off$/i;
	$lib = $ARGV[$k+1] if $lib eq "NA" and /^-lib$/i;
	$blastplus = $ARGV[$k+1] if /^-blastplus$/i;
	$k++;
	}
$blastplus = '' unless defined $blastplus;

die $usage if $lib eq "NA";

$/="\n>";
open Lib, "<$lib" or die $!;
while (<Lib>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	next unless $id =~ /\#LTR\//;
	next if $id =~ /_INT\#LTR\//;
	my ($ltr_len, $seq_len) = (0, 0);
	$seq_len = length $seq;
	if ($id =~ /_LTR#LTR/){
		print "$id\t1\t$seq_len\t$seq_len\n";
	} else {
		my $rand = int(rand(100000));
		open Temp1, ">temp_$rand.fa" or die $!;
		print Temp1 ">$id\n$seq";
		close Temp1;
		my ($Blast, @Blast);
		$Blast = `${blastplus}blastn -subject temp_$rand.fa -query temp_$rand.fa -dust no -outfmt 6`;
		@Blast = (split /\n/, $Blast);
		`rm temp_$rand.fa`;
		exit if @Blast < 1;
		foreach (@Blast){
			chomp;
			my ($qseqid, $sseqid, $iden, $length, $snp, $gap, $qstart, $qend, $sstart, $send) = (split);
			my $qlen = abs($qstart - $qend);
			next if $length < $min_align; #control alignment length
			next if $qlen > $seq_len*0.5; #skip self alignment
#print "$_\n";
			next unless ($sstart < $max_off and $seq_len - $qend < $max_off) or ($qstart < $max_off and $seq_len - $send < $max_off); #requires hit within 20bp of the input start/end
			print "$id\t$sstart\t$send\t$seq_len\n";
			print "$id\t$qstart\t$qend\t$seq_len\n";
			last;
			}
		}
	}
close Lib;

