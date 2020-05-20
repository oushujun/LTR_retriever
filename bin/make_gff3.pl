#!/usr/bin/env perl

##Generate gff3 file from pass list of LTR_retriever
##Usage: perl make_gff3.pl genome.fa LTR.pass.list
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version: 
#	1.0 04-09-2015
#	2.0 04-11-2018 Add ltr_digest support
#	3.0 05-19-2020 Reformat following GFF3 standard



use warnings;
use strict;

open FA, "<$ARGV[0]" or die "ERROR: $!";
open List, "<$ARGV[1]" or die "ERROR: $!";
open GFF, ">$ARGV[1].gff3" or die "ERROR: $!";

my $date=`date -u`;
chomp ($date);
print GFF "##gff-version 3\n##date $date
##seqid source repeat_class/superfamily start end score strand phase attributes\n";

#defined sequence ontology
my %SO = (Cent => "SO:0001797", TIR => "SO:0000208", MITE => "SO:0000338", Helitron => "SO:0000544", LINE => "SO:0000194", SINE => "SO:0000206", long_terminal_repeat => "SO:0000286", Low_complexity => "SO:0001005", LTRRT => "SO:0000186", "rDNA/spacer" => "SO:0001860", repeat_region => "SO:0000657", telomeric_repeat => "SO:0001496", subtelomere => "SO:0001997", nonLTR => "SO:0000189", target_site_duplication => "SO:0000434", primer_binding_site => "SO:0005850");

$/="\n>";
my $chr_info='';
my %seq_flag;
my $j=0;
while (<FA>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq)=(split /\n/, $_, 2);
	$seq=~s/\s+//g;
	my $len=length ($seq);
	$seq_flag{$id}=[$j, $len];
	$j++;
	}
$/="\n";

my $annotator="LTR_retriever";
my $i=1;
while (<List>){
	next if /^#/;
	next if /^\s+/;
	my ($chr, $element_start, $element_end, $element_length, $lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end, $rLTR_length, $seq_ID, $loc, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam);
	($loc, undef, $motif, $TSD, $lTSD, $rTSD, $IN, $sim, $strand, $supfam) = (split);
#chr03:19323647..19328532        pass    motif:TGCA      TSD:GTCGC       19323642..19323646      19328533..19328537      IN:19323854..19328325   99.03
	($chr, $lLTR_start, $rLTR_end)=($1, $2, $3) if $loc=~/^(.*):([0-9]+)..([0-9]+)$/;
	($lLTR_start, $rLTR_end)=($rLTR_end, $lLTR_start) if $rLTR_end<$lLTR_start;
	($lLTR_end, $rLTR_start)=($1-1, $2+1) if $IN=~/^IN:([0-9]+)..([0-9]+)$/;
	$motif=~s/motif://gi;
	$TSD=~s/TSD://gi;
	$lTSD=~s/\.\./\t/;
	$rTSD=~s/\.\./\t/;
	$element_start=(split /\s+/, $lTSD)[0];
	$element_end=(split /\s+/, $rTSD)[1];
	my $id = "$chr:$lLTR_start..$rLTR_end#LTR/$supfam";
	if ($TSD eq "NA"){
		$element_start=$lLTR_start;
		$element_end=$rLTR_end;
		}
	my $chr_ori=$chr;
	if (exists $seq_flag{$chr_ori}){
		print GFF "##sequence-region   $chr_ori 1 $seq_flag{$chr_ori}[1]\n"; #seqence length
		delete $seq_flag{$chr_ori};
		}
	my $info = "Name=$id;motif=$motif;tsd=$TSD;ltr_identity=$sim;Method=structural";
	print GFF "$chr\t$annotator\trepeat_region\t$element_start\t$element_end\t.\t$strand\t.\tID=repeat_region$i;Sequence_ontology=$SO{'repeat_region'};$info\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tID=lTSD_$i;Parent=repeat_region$i;Sequence_ontology=$SO{'target_site_duplication'};$info\n" unless $TSD eq "NA";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$lLTR_start\t$lLTR_end\t.\t$strand\t.\tID=lLTR_$i;Parent=repeat_region$i;Sequence_ontology=$SO{'long_terminal_repeat'};$info\n";
	print GFF "$chr\t$annotator\tLTR/$supfam\t$lLTR_start\t$rLTR_end\t.\t$strand\t.\tID=LTRRT_$i;Parent=repeat_region$i;Sequence_ontology=$SO{'LTRRT'};$info\n";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$rLTR_start\t$rLTR_end\t.\t$strand\t.\tID=rLTR_$i;Parent=repeat_region$i;Sequence_ontology=$SO{'long_terminal_repeat'};$info\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tID=rTSD_$i;Parent=repeat_region$i;Sequence_ontology=$SO{'target_site_duplication'};$info\n" unless $TSD eq "NA";
	print GFF "###\n";
	$i++;
	}
