#!/usr/bin/perl -w
use strict;

my ($index, $LINE, $DNA, $PlantP)=('', '', '', '');
($index, $LINE, $DNA, $PlantP)=@ARGV[0,1,2,3];
die "ERROR: Usage: perl cleanOut.pl index\n" if $index eq '';

`rm dummy060817.fa* $LINE.* $DNA.* $PlantP.* alluniRefprexp082813 Tpases020812DNA Tpases020812LINE`;

`rm $index.cat.gz $index.LTRlib.clust.clstr $index.LTRlib $index.LTRlib.fa.n* $index.ltrTE* $index.nmtf $index.prelib* $index.retriever.scn*`;

