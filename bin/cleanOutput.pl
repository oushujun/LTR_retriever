#!/usr/bin/perl -w
use strict;

my $index='';
$index=$ARGV[0];
die "ERROR: Usage: perl cleanOut.pl index\n" if $index eq '';

`rm $index.cat.gz $index.LTRlib.clust.clstr $index.LTRlib $index.LTRlib.fa.n* $index.ltrTE* $index.nmtf $index.prelib* $index.retriever.scn*`;

