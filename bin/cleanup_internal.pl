#!/usr/bin/env perl
use warnings;
use strict;
#Shujun Ou (shujun.ou.1@gmail.com) 03/07/2026

my $usage = "\n
Clean up solo LTR insertions nested inside internal regions of whole LTR elements.

Further info:
	Each LTR sequence from the LTR library is used as query to search the internal
	regions of whole elements. Matches covering >=95% of the LTR query within an
	internal region are removed. Terminal LTR regions are never searched or modified.
	After removal, elements with internal regions shorter than the threshold are discarded.

Usage: cleanup_internal.pl -in prelib.fa -scn scn.adj -ltrlib LTRonly.clust > prelib.cln
Options:	-in		Input whole-element library in FASTA format
		-scn		SCN file with LTR boundary info (retriever.all.scn format)
		-ltrlib		Clustered LTR-only sequences for masking
		-cov		Minimum coverage of the LTR query to consider as nested. Default: 0.95
		-minlen		Minimum length of the cleaned internal region to retain. Default: 100 (bp)
		-blastplus	Path to the blastn and makeblastdb programs.
		-threads	Threads to run blastn
\n";

my $IN = "";
my $SCN = "";
my $LTRLIB = "";
my $coverage = 0.95;
my $minlen = 100;
my $blastplus = "";
my $threads = 4;

my $k = 0;
foreach (@ARGV){
	$IN = $ARGV[$k+1] if /^-in$/i;
	$SCN = $ARGV[$k+1] if /^-scn$/i;
	$LTRLIB = $ARGV[$k+1] if /^-ltrlib$/i;
	$coverage = $ARGV[$k+1] if /^-cov$/i;
	$minlen = $ARGV[$k+1] if /^-minlen$/i;
	$blastplus = $ARGV[$k+1] if /^-blastplus$/i;
	$threads = $ARGV[$k+1] if /^-threads$/i;
	$k++;
}

die $usage unless -s $IN and -s $SCN and -s $LTRLIB;

# Step 1: Read SCN file to get LTR boundary info per element
# Keys use ascending coordinate order: "$chr:$start..$end"
my %bound; # "$chr:$start..$end" => [lLTR_len, rLTR_len]
open SCN, "<$SCN" or die "ERROR: Cannot open $SCN: $!\n";
while (<SCN>){
	next if /^#/;
	s/^\s+//;
	my @f = split /\s+/;
	next unless @f >= 9;
	my ($start, $end, $lLTR_len, $rLTR_len, $chr) = @f[0, 1, 5, 8, 11];
	next unless defined $chr and defined $lLTR_len and defined $rLTR_len;
	# Ensure ascending order
	($start, $end) = ($end, $start) if $start > $end;
	my $key = "$chr:$start..$end";
	$bound{$key} = [$lLTR_len, $rLTR_len] unless exists $bound{$key}; # first entry wins
}
close SCN;

# Step 2: Read whole-element FASTA (prelib)
my %seq;     # id => sequence
my @order;   # preserve input order
my %lLTR;    # id => left LTR subsequence
my %rLTR;    # id => right LTR subsequence
my %internal;# id => internal region subsequence

open IN, "<$IN" or die "ERROR: Cannot open $IN: $!\n";
$/ = "\n>";
while (<IN>){
	s/>//g;
	s/^\s+//;
	my ($id, $seq) = (split /\n/, $_, 2);
	next unless defined $id and defined $seq;
	$id =~ s/\s+//g;
	$seq =~ s/\s+//g;
	$seq{$id} = $seq;
	push @order, $id;

	# Parse header to get chr and coordinates
	my ($chr, $coords) = ("", "");
	if ($id =~ /^([^:]+):([0-9]+\.\.[0-9]+)/) {
		($chr, $coords) = ($1, $2);
	} else {
		# No coordinates in header — keep element as-is (no boundary info)
		$internal{$id} = $seq;
		$lLTR{$id} = "";
		$rLTR{$id} = "";
		next;
	}

	# Normalize coordinates to ascending for lookup
	my $lookup_key = $coords;
	if ($lookup_key =~ /^(\d+)\.\.(\d+)$/ && $1 > $2) {
		$lookup_key = "$2..$1";
	}
	$lookup_key = "$chr:$lookup_key";

	if (exists $bound{$lookup_key}) {
		my ($llen, $rlen) = @{$bound{$lookup_key}};
		my $total = length($seq);
		my $int_len = $total - $llen - $rlen;

		if ($int_len < 1) {
			# Element is all LTR, no internal region — keep as-is
			$lLTR{$id} = $seq;
			$rLTR{$id} = "";
			$internal{$id} = "";
			next;
		}

		$lLTR{$id} = substr($seq, 0, $llen);
		$internal{$id} = substr($seq, $llen, $int_len);
		$rLTR{$id} = substr($seq, $total - $rlen);
	} else {
		# No boundary info found — keep element as-is
		$internal{$id} = $seq;
		$lLTR{$id} = "";
		$rLTR{$id} = "";
	}
}
$/ = "\n";
close IN;

# Step 3: Write internal regions to temp file for blast
my $tmp_int = "$IN.internal.tmp.$$";
open TMP, ">$tmp_int" or die "ERROR: Cannot write $tmp_int: $!\n";
my $has_internal = 0;
for my $id (@order) {
	next unless exists $internal{$id} and length($internal{$id}) > 0;
	print TMP ">$id\n$internal{$id}\n";
	$has_internal = 1;
}
close TMP;

# If no internal regions to search, output all elements unchanged
unless ($has_internal) {
	for my $id (@order) {
		next unless exists $seq{$id};
		print ">$id\n$seq{$id}\n";
	}
	unlink $tmp_int;
	exit 0;
}

# Step 4: Build blast database from internal regions and search with LTR library
`${blastplus}makeblastdb -in $tmp_int -dbtype nucl 2>/dev/null`;

# Get LTR query lengths
my %qlen;
open LTR, "<$LTRLIB" or die "ERROR: Cannot open $LTRLIB: $!\n";
$/ = "\n>";
while (<LTR>) {
	s/>//g;
	s/^\s+//;
	my ($qid, $qseq) = (split /\n/, $_, 2);
	next unless defined $qid and defined $qseq;
	$qid =~ s/\s+//g;
	$qseq =~ s/\s+//g;
	$qlen{$qid} = length($qseq);
}
$/ = "\n";
close LTR;

# Step 5: Blast LTR library against internal regions
my $exec = "timeout -s KILL 300s ${blastplus}blastn -query $LTRLIB -db $tmp_int -outfmt 6 -num_threads $threads 2>/dev/null";
my @blast = qx($exec);

# Collect regions to remove from each internal sequence
# For each subject (internal region), track ranges to excise
my %remove; # subject_id => [[start, end], ...]
for my $hit (@blast) {
	chomp $hit;
	my @f = split /\t/, $hit;
	next unless @f >= 12;
	my ($query, $subject, $iden, $len, $sbj_start, $sbj_end) = @f[0, 1, 2, 3, 8, 9];

	# Check coverage: alignment length / query length >= threshold
	next unless exists $qlen{$query};
	next unless $len / $qlen{$query} >= $coverage;

	# Normalize subject coordinates
	($sbj_start, $sbj_end) = ($sbj_end, $sbj_start) if $sbj_start > $sbj_end;

	push @{$remove{$subject}}, [$sbj_start, $sbj_end];
}

# Step 6: Remove matching portions from internal regions
for my $id (keys %remove) {
	next unless exists $internal{$id};
	my $int_seq = $internal{$id};
	my $int_len = length($int_seq);

	# Sort ranges by start position (descending) to remove from end first
	my @ranges = sort { $b->[0] <=> $a->[0] } @{$remove{$id}};

	for my $range (@ranges) {
		my ($rs, $re) = @$range;
		# Clamp to sequence bounds (1-based coords from blast)
		$rs = 1 if $rs < 1;
		$re = $int_len if $re > $int_len;
		next if $rs > $int_len;

		# Remove the region (convert 1-based to 0-based for substr)
		my $before = ($rs > 1) ? substr($int_seq, 0, $rs - 1) : "";
		my $after = ($re < length($int_seq)) ? substr($int_seq, $re) : "";
		$int_seq = $before . $after;
	}

	$internal{$id} = $int_seq;
}

# Step 7: Reassemble and output
for my $id (@order) {
	next unless exists $seq{$id};

	my $new_seq;
	if (exists $lLTR{$id} and (length($lLTR{$id}) > 0 or length($rLTR{$id}) > 0)) {
		# Has boundary info — reassemble
		my $cleaned_int = $internal{$id} // "";

		# Check if cleaned internal is long enough
		if (length($cleaned_int) < $minlen and length($lLTR{$id}) > 0) {
			# Internal too short after cleanup — discard element
			next;
		}

		$new_seq = $lLTR{$id} . $cleaned_int . $rLTR{$id};
	} else {
		# No boundary info — output original
		$new_seq = $seq{$id};
	}

	# Final length check
	next if length($new_seq) < $minlen;
	print ">$id\n$new_seq\n";
}

# Cleanup temp files
unlink $tmp_int;
unlink "$tmp_int.ndb", "$tmp_int.nhr", "$tmp_int.nin", "$tmp_int.njs",
       "$tmp_int.not", "$tmp_int.nsq", "$tmp_int.ntf", "$tmp_int.nto";