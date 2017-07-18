##To annotate the gff file with information from the library
##Usage: perl annotate_gff.pl lib.fa gff > anno.gff
##Shujun Ou (oushujun@msu.edu)

#!/usr/bin/perl -w
use strict;

my $usage="perl annotate_gff.pl lib.fa gff > anno.gff\n";
my $lib=$ARGV[0];
my $gff=$ARGV[1];

open GFF, "<$gff" or die "ERROR: $!";
while (<GFF>){
	my $id='';
	$id=$1 if /\"(.*)\"/;
	$id=~s/Motif://;
	$id=~s/-int//;
	my $name='';
	$name=`grep $id $lib` if $id ne '';
	$name=(split /\s+/, $name)[0];
	my $anno=$1 if $name=~/#(.*)/;
	$_=~s/Motif/$anno/;
	print "$_";
	}
