#!/usr/bin/env perl
use warnings;
use strict;

##Generate gff3 file from RepeatMasker .out file of LTR_retriever
##Usage: perl make_gff3.pl genome.fa.out
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version: 	2.0 05-19-2020
#		1.0 01-27-2018


my $usage = "\n\tperl make_gff3_with_RMout.pl RM.out\n\n";
die $usage unless -s $ARGV[0];
my $date=`date -u`;
chomp ($date);
open RMout, "sort -sV -k5,5 $ARGV[0]|" or die "ERROR: $!";
open GFF, ">$ARGV[0].gff" or die "ERROR: $!";
print GFF "##gff-version 3\n##date $date
##seqid source repeat_class/superfamily start end sw_score strand phase attributes\n";

#defined sequence ontology
my %class = ("Cent/CentC" => "Cent", "Centro/tandem" => "Cent", "DNAauto/CACTA" => "TIR", "DNAauto/CACTG" => "TIR", "DNAauto/hAT" => "TIR", "DNAauto/Helitron" => "Helitron", "DNAauto/MLE" => "TIR", "DNAauto/MULE" => "TIR", "DNAauto/PILE" => "TIR", "DNAauto/POLE" => "TIR", "DNA/DTA" => "TIR", "DNA/DTC" => "TIR", "DNA/DTH" => "TIR", "DNA/DTM" => "TIR", "DNA/DTT" => "TIR", "DNA/Helitron" => "Helitron", "DNAnona/CACTA" => "TIR", "DNAnona/CACTG" => "TIR", "DNAnona/hAT" => "TIR", "DNAnona/Helitron" => "Helitron", "DNAnona/MLE" => "TIR", "DNAnona/MULE" => "TIR", "DNAnona/MULEtir" => "TIR", "DNAnona/PILE" => "TIR", "DNAnona/POLE" => "TIR", "DNAnona/Tourist" => "TIR", "DNAnona/unknown" => "DNATE", "Evirus/ERTBV-A" => "TE", "Evirus/ERTBV-B" => "TE", "Evirus/ERTBV-C" => "TE", "Evirus/ERTBV" => "TE", "knob/knob180" => "knob", "knob/TR-1" => "knob", "LINE/L1" => "LINE", "LINE/RTE" => "LINE", "LINE/unknown" => "LINE", "LTR/Copia" => "LTRRT", "LTR/CRM" => "LTRRT", "LTR/Gypsy" => "LTRRT", "LTR/Solo" => "LTRRT", "LTR/TRIM" => "LTRRT", "LTR/unknown" => "LTRRT", "MITE/DTA" => "MITE", "MITE/DTC" => "MITE", "MITE/DTH" => "MITE", "MITE/DTM" => "MITE", "MITE/DTT" => "MITE", "MITE/Stow" => "MITE", "MITE/Tourist" => "MITE", "Satellite/rice" => "satellite_DNA", "SINE/unknown" => "SINE", "subtelomere/4-12-1" => "subtelomere");

my %SO = (repeat_region => "SO:0000657", TE => "SO:0000101", DNATE => "SO:0000182", TIR => "SO:0000208", MITE => "SO:0000338", Helitron => "SO:0000544", retrotransposon => "SO:0000180", nonLTR => "SO:0000189", LINE => "SO:0000194", SINE => "SO:0000206", LTRRT => "SO:0000186", target_site_duplication => "SO:0000434", primer_binding_site => "SO:0005850", long_terminal_repeat => "SO:0000286", Low_complexity => "SO:0001005", "rDNA/spacer" => "SO:0001860", telomeric_repeat => "SO:0001496", subtelomere => "SO:0001997", Cent => "SO:0001797", satellite_DNA => "SO:0000005");

my %seq_flag;
my $annotator="LTR_retriever";
my $i = 0; #annotation count
while (<RMout>){
	s/^\s+//;
	my ($SW_score, $div, $iden, $chr, $chr_len, $element_start, $element_end, $element_length, $left_len, $strand, $TE_ID, $TE_class);
	($SW_score, $div, $chr, $element_start, $element_end, $left_len, $strand, $TE_ID, $TE_class)=(split)[0,1,4,5,6,7,8,9,10];
	next unless defined $SW_score and $SW_score =~ /[0-9]+/;
	my $so = $TE_class;
	$so = $class{$TE_class} if exists $class{$TE_class};

	$TE_ID=~s/_INT-int$/_INT/;
	$element_length=$element_end-$element_start+1;
	next if $SW_score<300 and $element_length<80;
	$strand="-" if $strand eq "C";
	unless (exists $seq_flag{$chr}){
		$seq_flag{$chr}=$chr;
		$left_len=~s/[\(\) ]+//g;
		$left_len=0 unless $left_len=~/^[0-9]+$/;
		$chr_len=$element_end+$left_len;
		print GFF "##sequence-region $chr 1 $chr_len\n";
		}
	$iden = 1 - $div/100;
	my $info = "ID=TE_annot_$i;Sequence_ontology=$SO{$so};Name=$TE_ID;Identity=$iden;Method=homology";
	print GFF "$chr\t$annotator\t$TE_class\t$element_start\t$element_end\t$SW_score\t$strand\t.\t$info\n";
	$i++; #annotation count
	}
close RMout;
close GFF;
