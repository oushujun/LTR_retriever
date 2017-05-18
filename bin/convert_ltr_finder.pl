#convert LTR_finder screen output to LTR-modeler input format
#Usage: perl convert.pl LTR_finder_scn > out.scn
#Author:Shujun OU (oushujun@msu.edu) 03/07/2015


#!/usr/bin/perl -w
use strict;

print "#this summary file is created from LTR_finder screen output by convert_ltr_finder.pl (Author: Shujun Ou, email: oushujun\@msu.edu 2/14/2015)
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif superfamily family age(ya)\n";
my ($from, $to, $chr, $LTR_len, $lLTR_len, $rLTR_len, $lLTR_end, $rLTR_str, $similarity, $TSD, $motif, $direction, $seq_id);
$from=$to=$LTR_len=$from=$lLTR_end=$lLTR_len=$rLTR_str=$to=$rLTR_len=$similarity=$chr=$direction=$TSD=$motif=' ';
$seq_id=-1;
while (<>){
	($chr=$1 and $seq_id++) if /^>Sequence:\s?(.+)\s+Len:[0-9]+/i;
	($from, $to, $LTR_len, $direction)=($1, $2, $3, $4) if /Location :\s+([0-9]+)\s+-\s+([0-9]+)\s+Len:\s+([0-9]+)\s+Strand:(\+|-)/i;
	$similarity=$1 if /Score\s+:\s+[0-9]+\s+\[LTR region similarity:([0-9.]+)\]/i;
	$lLTR_len=$1 and $lLTR_end=$from+$lLTR_len-1 if /5'-LTR.*Len:\s+([0-9]+)/i;
	$rLTR_len=$1 and $rLTR_str=$to-$rLTR_len+1 if /3'-LTR.*Len:\s+([0-9]+)/i;
	$motif=$1 if /5'-TG\s+:\s+(.*)/;
	$motif="$motif".",$1" if /3'-CA\s+:\s+(.*)/i;
	$motif=~s/ //g if defined $motif;
	$TSD="$5\t$1..$2, $3..$4" if /TSR\s+:\s+([0-9]+)\s+-\s+([0-9]+)\s+,\s+([0-9]+)\s+-\s+([0-9]+)\s+\[([ATCGN]+)\]/i;
	if (/^PPT\s+:\s+\[/i){
		print "$from $to $LTR_len $from $lLTR_end $lLTR_len $rLTR_str $to $rLTR_len $similarity $chr $chr $direction $TSD $motif\n";
		$from=$to=$LTR_len=$lLTR_len=$rLTR_len=$lLTR_end=$rLTR_str=$similarity=$TSD=$motif=$direction="NA";
		}
	}

