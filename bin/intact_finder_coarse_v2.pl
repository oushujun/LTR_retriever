#!/usr/bin/env perl
use warnings;
use strict;

#Indentify intact LTRs in RM out files
#Shujun Ou (shujun.ou.1@gmail.com)
#07-25-2017
#update: 02-22-2025, facilitated by ChatGPT

my $usage = "Usage: perl $0 lib.info RepeatMasker.out > intact_list\n";
die "\n$usage\n\n" unless @ARGV==2;

my $lib_file = shift @ARGV;
my $rm_file  = shift @ARGV;

# parameters
my $min_internal = 1; # min # of INT segments
my $max_internal = 5; # max # of INT segments
my $max_size = 100000; # max size of an intact element (nested included)


#############################################
# 1. Build the lib.info Database
#############################################
# Each record is keyed by the element name.
# We expect intact elements to have exactly two entries and their name should NOT contain "_LTR"
my %libdb;
open(my $lib_fh, "<", $lib_file) or die "Cannot open $lib_file: $!";
while (<$lib_fh>) {
    chomp;
    next if /^\s*$/;
    my ($name, $start, $end, $full_length) = split;
    $name =~ s/#.*//;
    push @{ $libdb{$name} }, { start => $start, end => $end, full_length => $full_length };
}
close $lib_fh;
# For elements with two entries, sort them so that the first is the left LTR and the second is the right LTR.
for my $name (keys %libdb) {
    if (@{ $libdb{$name} } > 1) {
        @{ $libdb{$name} } = sort { $a->{start} <=> $b->{start} } @{ $libdb{$name} };
    }
}

#############################################
# 2. Read the RepeatMasker .out File and Group Lines
#############################################
# The RM file columns are:
# SW_score, perc_div, perc_del, perc_ins, query_sequence, query_begin, query_end, query_remain,
# strand, matching_repeat, repeat_class/family, repeat_begin, repeat_end, repeat_remain, ID
#
# We separate RM lines for intact candidates (those that have two lib.info entries and whose name does not include "_LTR")
# from the rest (LTR-only or split cases).
my %rm_intact;  # key: consensus name, value: array ref of RM lines
my @rm_split;   # for split LTR-INT candidates (or LTR-only)
open(my $rm_fh, "<", $rm_file) or die "Cannot open $rm_file: $!";
while (<$rm_fh>) {
    next unless /LTR/;
    s/[()]//g;
    s/^\s+//g;
    next unless /^[0-9]+/;
    s/-int\s+/\t/;
    my @fields = split;
    my $score  = $fields[0];
    my $qstart = $fields[5];
    my $qend   = $fields[6];

    # Filter out low-score or very short alignments
    next if $score <= 300;
    #next if ($qend - $qstart) <= 100;

    # Entries annotated by intact elements must be in libdb with two entries and must not include "_LTR" in the ID.
    my $match = $fields[9];  # matching_repeat (the consensus ID)
    if ( exists $libdb{$match} and @{ $libdb{$match} } == 2 and $match !~ /_LTR/ ) {
         push @{ $rm_intact{$match} }, $_;
    }
    else {
         push @rm_split, $_;
    }
}
close $rm_fh;

#####################################################
# 3. Process Intact Candidates Using RM Repeat Coordinates
#####################################################
# For each intact candidate we require that (collectively) the RM pieces cover some part of both LTR regions.
# Since the RM coordinates (repeat_begin, repeat_end) are in the same space as the lib.info file,
# we check for an overlap of at least 1 bp.
print "TE_ID\tChr\tStart\tEnd\tLTR_start\tLTR_end\tINT_count\tLTR_left_start\tLTR_left_end\tLTR_right_start\tLTR_right_end\n";
my $INT_count = 1;
foreach my $consensus (keys %rm_intact) {
    # Get the two lib.info entries: left and right LTR regions.
    my ($ltr_left, $ltr_right) = @{ $libdb{$consensus} };
    my ($found_left, $found_right) = (0, 0);
    foreach my $line (@{ $rm_intact{$consensus} }) {
         my @fields = split(/\s+/, $line);
         my $strand = $fields[8];
	 my $chr = $fields[4];
	 my $start = $fields[5];
	 my $end = $fields[6];
	 my ($rm_start, $rm_end);
	 if ($strand =~ /[C\-]/){
	         ($rm_start, $rm_end) = ($fields[13], $fields[12]);
	 } else {
		 ($rm_start, $rm_end) = ($fields[11], $fields[12]);
	 }
	 ($rm_start, $rm_end) = ($rm_end, $rm_start) if ($rm_start > $rm_end);
	 #	 print "${line}$consensus\t$rm_start\t$rm_end\t$ltr_left->{start}\t$ltr_left->{end}\t$ltr_right->{start}\t$ltr_right->{end}\n"; #test

	 if ($rm_start <= $ltr_left->{end} and $rm_end >= $ltr_right->{start}) {
		 print "$consensus\t$chr\t$start\t$end\t$rm_start\t$rm_end\t$INT_count\t$ltr_left->{start}\t$ltr_left->{end}\t$ltr_right->{start}\t$ltr_right->{end}\n";
	 }
    }
}

#####################################################
# 4. Process Split Candidates (LTR-INT[-INT...]-LTR)
#####################################################
# Revised approach:
#   1. Group RM lines (from @rm_split) by chromosome and gap (<300 bp)
#      regardless of strand.
#   2. Within each group, search every contiguous sub‐array for an intact element.
#      The candidate must:
#         - Have at least three segments (LTR, INT(s), LTR).
#         - Begin and end with an LTR segment (non-INT).
#         - Have gaps between consecutive segments <300 bp.
#         - Have first and last segments on the same strand.
#         - Have matching LTR IDs (after stripping any “-int” suffix).
#         - Have divergence (field[1]) difference between the LTRs <4.
#         - Contain an allowed number of internal segments (default: 1 to 3).
#         - Total length less than 100kb
#         - The candidate's left LTR RM line's query_end must be within 50 bp of the left library LTR's end,
#		and the right LTR RM line's query_begin must be within 50 bp of the right library LTR's start.
#

# First, group the RM lines by chromosome and contiguous distance (<300bp).
# Note: We assume @rm_split is already sorted by chromosome and coordinates.
my @groups;
my @current_group;
my $prev_fields;
foreach my $line (@rm_split) {
    chomp $line;
    my @fields = split(/\s+/, $line);
    # Fields: [4]=chr, [5]=query_begin, [6]=query_end, [8]=strand, [9]=matching_repeat, [1]=perc_div
    if (!@current_group) {
        @current_group = ($line);
        $prev_fields = \@fields;
        next;
    }
    # Check if same chromosome and gap (<300 bp between previous query_end and current query_begin)
    if ($fields[4] eq $prev_fields->[4] and (($fields[5] - $prev_fields->[6]) < 300)) {
        push @current_group, $line;
    } else {
        push @groups, [ @current_group ];
        @current_group = ($line);
    }
    $prev_fields = \@fields;
}
push @groups, [ @current_group ] if @current_group;

# Process each group to search for candidate intact elements.
foreach my $group_ref (@groups) {
    process_candidate_group($group_ref);
}

sub process_candidate_group {
    my ($group_ref) = @_;
    # Convert each line in the group into an array of fields.
    my @lines;
    foreach my $line (@$group_ref) {
        chomp $line;
        my @f = split(/\s+/, $line);
        push @lines, \@f;
    }
    my $n = scalar(@lines);
    # Iterate over all possible contiguous sub-arrays.
    for (my $i = 0; $i < $n; $i++) {
        # Candidate must start with an LTR (non-INT) segment.
        my $ltr_left = $lines[$i]->[9];
        next unless ($ltr_left =~ /LTR/i and $ltr_left !~ /int/i);
	my $num_internal = 0;
        for (my $j = $i+1; $j < $n; $j++) {
            # Candidate details.
            my $chr   = $lines[$i]->[4];
            my $start = $lines[$i]->[5];
            my $end   = $lines[$j]->[6];
	    # Search distance less than 100 kb from the left LTR
	    next if $end - $start + 1 > $max_size;
	    # count INT number
	    $num_internal++ if $lines[$j]->[9] =~ /int/i;
            # Candidate must end with an LTR segment.
            my $ltr_right = $lines[$j]->[9];
            next unless ($ltr_right =~ /LTR/i and $ltr_right !~ /int/i);
            # Check that LTR IDs match
            next unless ($ltr_left eq $ltr_right);
            # Enforce strand consistency for the candidate.
            next unless ($lines[$i]->[8] eq $lines[$j]->[8]);
	    # The candidate's left LTR RM line's query_end must be within 50 bp of the left library LTR's end,
            # and the right LTR RM line's query_begin must be within 50 bp of the right library LTR's start.
	    # works for both + or - strand
	    my ($left_end_off, $right_start_off) = ($lines[$i]->[13], $lines[$j]->[11]);
	    next if $left_end_off > 50;
	    next if $right_start_off > 50;
	    my ($ltr_left_start, $ltr_left_end, $ltr_right_start, $ltr_right_end);
	    if ($lines[$i]->[8] =~ /[C\-]/){
		    ($ltr_left_start, $ltr_right_end) = ($lines[$i]->[12], $lines[$j]->[13]);
	    } else {
		    ($ltr_left_start, $ltr_right_end) = ($lines[$i]->[11], $lines[$j]->[12]);
	    }
            # Check divergence difference between the left and right LTR (<4).
            my $div_left  = $lines[$i]->[1];
            my $div_right = $lines[$j]->[1];
            next if (abs($div_left - $div_right) >= 4);
	    #print "\n@{$lines[$i]}\n@{$lines[$j]}\n"; #test
	    #print "($num_internal < $min_internal or $num_internal > $max_internal\n"; #test
            # Number of internal segments must be within the allowed range.
	    next if ($num_internal < $min_internal or $num_internal > $max_internal);
	    my ($lib_start, $lib_end) = (@{ $libdb{$ltr_left} }[0]->{start}, @{ $libdb{$ltr_left} }[0]->{end});
            print "$ltr_left\t$chr\t$start\t$end\t$ltr_left_start\t$ltr_right_end\t$num_internal\t$lib_start\t$lib_end\t$lib_start\t$lib_end\n";
        }
    }
}

