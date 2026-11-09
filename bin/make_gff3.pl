#!/usr/bin/env perl

##Generate gff3 file from pass list of LTR_retriever
##Usage: perl make_gff3.pl genome.fa LTR.pass.list [LTR.pass.list.pbs_ppt]
##Author: Shujun Ou (oushujun@msu.edu), Department of Horticulture, Michigan State University
##Version:
#	1.0 04-09-2015
#	2.0 04-11-2018 Add ltr_digest support
#	3.0 05-19-2020 Reformat following GFF3 standard
#	4.0 11-08-2026 Add PBS, PPT, and strand evidence from find_PBS_PPT.pl



use warnings;
use strict;

my $usage = "\n\tperl make_gff3.pl genome.fa LTR.pass.list [LTR.pass.list.pbs_ppt]\n\n";
open FA, "<$ARGV[0]" or die "ERROR: $usage";
open List, "<$ARGV[1]" or die "ERROR: $usage";
open GFF, ">$ARGV[1].gff3" or die "ERROR: $!";

#PBS, PPT, and strand evidence of the elements from find_PBS_PPT.pl (optional), by element location
my %pbsppt;
if (defined $ARGV[2]){
	open PP, "<$ARGV[2]" or die "ERROR: $usage";
	my @col;
	while (<PP>){
		chomp;
		if (s/^#//){
			@col = split /\t/;
			next;
			}
		my %r;
		@r{@col} = split /\t/;
		next unless defined $r{LTR_loc} and $r{LTR_loc} =~ /^(.*):([0-9]+)\.\.([0-9]+)$/;
		my ($chr, $s, $e) = ($1, $2, $3);
		($s, $e) = ($e, $s) if $s > $e;
		$pbsppt{"$chr:$s..$e"} = \%r;
		}
	close PP;
	}

my $date=`date -u`;
chomp ($date);
print GFF "##gff-version 3\n##date $date
##ltr_identity: Sequence identity (0-1) between the left and right LTR region.
##tsd: target site duplication.
##seqid source sequence_ontology start end score strand phase attributes\n";

#defined sequence ontology
my %SO = (repeat_region => "SO:0000657", long_terminal_repeat => "SO:0000286", target_site_duplication => "SO:0000434", primer_binding_site => "SO:0005850", RR_tract => "SO:0000435", "LTR_retrotransposon" => "SO:0000186", "Copia_LTR_retrotransposon" => "SO:0002264", "Gypsy_LTR_retrotransposon" => "SO:0002265", "Bel_Pao_LTR_retrotransposon" => "SO:0002266");

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
	my $id = $strand eq '-' ? "$chr:$rLTR_end..$lLTR_start" : "$chr:$lLTR_start..$rLTR_end"; # strand sensitive $id
	my $so = "LTR_retrotransposon";
	$so = "Copia_LTR_retrotransposon" if $supfam eq "Copia";
	$so = "Gypsy_LTR_retrotransposon" if $supfam eq "Gypsy";
	$so = "Bel_Pao_LTR_retrotransposon" if $supfam =~ /Bel.Pao/;
	if ($TSD eq "NA"){
		$element_start=$lLTR_start;
		$element_end=$rLTR_end;
		}
	my $chr_ori=$chr;
	if (exists $seq_flag{$chr_ori}){
		print GFF "##sequence-region   $chr_ori 1 $seq_flag{$chr_ori}[1]\n"; #seqence length
		delete $seq_flag{$chr_ori};
		}
	my $info = "Name=$id;classification=LTR/$supfam";
	my $info2 = "ltr_identity=$sim;method=structural;motif=$motif;tsd=$TSD";

	#strand evidence (domain, PBS, PPT), and the PBS and PPT that support the strand, if the table has the list's strand
	my ($info3, @pp) = ('');
	my $r = $pbsppt{"$chr:$lLTR_start..$rLTR_end"};
	if (defined $r and defined $r->{strand} and $r->{strand} eq $strand){
		($info3 = ";strand_evidence=$r->{strand_evidence}") =~ s/\+/,/g;
		($info3 .= ";strand_conflict=$r->{strand_conflict}") =~ s/\+/,/g unless $r->{strand_conflict} eq 'NA';
		my $h = $strand eq '+' ? 'plus' : 'minus';
		if ($strand ne '?' and $r->{call_pbs} eq $strand){
			(my $trna = $r->{"pbs_${h}_trna"}) =~ s/_[0-9]+$//; #tRNA isotype
			push @pp, [$r->{"pbs_${h}_start"}, "$chr\t$annotator\tprimer_binding_site\t".$r->{"pbs_${h}_start"}."\t".$r->{"pbs_${h}_end"}."\t.\t$strand\t.\tID=PBS_$i;Parent=LTRRT_$i;$info;sequence_ontology=$SO{'primer_binding_site'};$info2;trna=$trna;trna_offset=".$r->{"pbs_${h}_trna_offset"}.";ltr_offset=".$r->{"pbs_${h}_gap"}.";edist=".$r->{"pbs_${h}_edist"}."\n"];
			}
		if ($strand ne '?' and $r->{call_ppt} eq $strand){
			push @pp, [$r->{"ppt_${h}_start"}, "$chr\t$annotator\tRR_tract\t".$r->{"ppt_${h}_start"}."\t".$r->{"ppt_${h}_end"}."\t.\t$strand\t.\tID=PPT_$i;Parent=LTRRT_$i;$info;sequence_ontology=$SO{'RR_tract'};$info2;ltr_offset=".$r->{"ppt_${h}_gap"}."\n"];
			}
		}

	print GFF "$chr\t$annotator\trepeat_region\t$element_start\t$element_end\t.\t$strand\t.\tID=repeat_region_$i;$info;sequence_ontology=$SO{'repeat_region'};$info2$info3\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$lTSD\t.\t$strand\t.\tID=lTSD_$i;Parent=repeat_region_$i;$info;sequence_ontology=$SO{'target_site_duplication'};$info2\n" unless $TSD eq "NA";
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$lLTR_start\t$lLTR_end\t.\t$strand\t.\tID=lLTR_$i;Parent=repeat_region_$i;$info;sequence_ontology=$SO{'long_terminal_repeat'};$info2\n";
	print GFF "$chr\t$annotator\t$so\t$lLTR_start\t$rLTR_end\t.\t$strand\t.\tID=LTRRT_$i;Parent=repeat_region_$i;$info;sequence_ontology=$SO{$so};$info2$info3\n";
	print GFF $_->[1] foreach sort { $a->[0] <=> $b->[0] } @pp;
	print GFF "$chr\t$annotator\tlong_terminal_repeat\t$rLTR_start\t$rLTR_end\t.\t$strand\t.\tID=rLTR_$i;Parent=repeat_region_$i;$info;sequence_ontology=$SO{'long_terminal_repeat'};$info2\n";
	print GFF "$chr\t$annotator\ttarget_site_duplication\t$rTSD\t.\t$strand\t.\tID=rTSD_$i;Parent=repeat_region_$i;$info;sequence_ontology=$SO{'target_site_duplication'};$info2\n" unless $TSD eq "NA";
	print GFF "###\n";
	$i++;
	}
