#!/usr/bin/perl -w
use strict;

my ($index, $LINE, $DNA, $PlantP)=('', '', '', '');
($index, $LINE, $DNA, $PlantP)=@ARGV[0,1,2,3];
die "ERROR: Usage: perl cleanOut.pl index\n" if $index eq '';

`rm $LINE.* $DNA.* $PlantP.* alluniRefprexp082813* Tpases020812DNA* Tpases020812LINE* 2>/dev/null`;

`rm $index.cat.gz $index.LTRlib.clust.clstr $index.LTRlib $index.LTRlib.raw $index.LTRID.list $index.LTRlib.fa.n* $index.ltrTE* $index.nmtf $index.prelib* $index.nmtf.prelib $index.retriever.scn* $index.retriever.all.scn.list $index.*exclude* 2>/dev/null`;

