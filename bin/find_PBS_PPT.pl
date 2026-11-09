#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Temp qw(tempdir);
use POSIX qw(ceil _exit);
use IO::Handle;

##Find the primer binding site (PBS) and polypurine tract (PPT) of intact LTR retrotransposons and call their strand.
##
##Credits: the PBS and PPT searches are LTRdigest's algorithms, reimplemented in Perl from GenomeTools 1.6.5
##(src/ltr/ltrdigest_pbs_visitor.c, src/ltr/ltrdigest_ppt_visitor.c, src/extended/swalign.c, src/extended/hmm.c,
##src/extended/alignment.c) with LTRdigest's default parameters and tie-breaking rules, so GenomeTools is not needed.
##	Steinbiss S, Willhoeft U, Gremme G, Kurtz S (2009). Fine-grained annotation and classification of de novo predicted
##	LTR retrotransposons. Nucleic Acids Res 37(21):7002-7013. https://doi.org/10.1093/nar/gkp759
##	Gremme G, Steinbiss S, Kurtz S (2013). GenomeTools: a comprehensive software library for efficient processing of
##	structured genome annotations. IEEE/ACM Trans Comput Biol Bioinform 10(3):645-656.
##The GenomeTools code this follows carries these notices:
##	Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
##	Copyright (c) 2005-2009 Gordon Gremme <gordon@gremme.org>
##	Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
##	Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
##	Copyright (c) 2005-2015 Center for Bioinformatics, University of Hamburg
##
##	Permission to use, copy, modify, and distribute this software for any
##	purpose with or without fee is hereby granted, provided that the above
##	copyright notice and this permission notice appear in all copies.
##
##	THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
##	WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
##	MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
##	ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
##	WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
##	ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
##	OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
##The tRNA library (database/tRNA.fa) is built from the tRNA gene sets distributed with LTR_FINDER (Xu and Wang 2007,
##Nucleic Acids Res 35:W265), tRNAscan-SE predictions (Lowe and Eddy 1997, Nucleic Acids Res 25:955), plus the 3' CCA.
##
##The searches:
##	PBS: Smith-Waterman alignment (match 5, mismatch -10, indel -20) of the 61 nt around the 5' LTR end against the
##	reverse complement of each tRNA; a hit has <= 1 edit, 11-30 nt, starts <= 5 nt from the LTR end and <= 5 nt from
##	the tRNA 3' end. Only tRNAs that can give such a hit (from an 11-mer index) are aligned.
##	PPT: Viterbi decoding of LTRdigest's 4-state purine-tract HMM on the 61 nt around the 3' LTR start; a hit is a
##	tract of 8-30 nt, the one ending closest to the 3' LTR.
##Both strand hypotheses are searched. The strand call was calibrated on 52,291 intact LTR-RTs of maize B73: 99.6%
##precision at 95.9% coverage on the elements whose strand was known from their protein domains and the reference annotation.
##	A PBS counts if its tRNA 3' offset is <= 2; a PPT counts if it ends between 2 nt inside the 3' LTR and 10 nt before it.
##	If both hypotheses have one, the better (PBS: length - 2 x edit distance; PPT: length - 2 x |distance to the LTR|)
##	must lead by >= 2. PBS and PPT agreeing, or only one of them -> that strand; opposite -> no call.
##The final strand keeps the protein-domain strand of the list unless PBS and PPT agree on the other strand.
##Usage: perl find_PBS_PPT.pl -genome genome.fa -list genome.pass.list [-trna tRNA.fa] [-threads 1] > genome.pass.list.pbs_ppt

my $usage = "
Find the primer binding site (PBS) and polypurine tract (PPT) of intact LTR-RTs and call their strand.
	perl find_PBS_PPT.pl -genome genome.fa -list genome.pass.list [options] > genome.pass.list.pbs_ppt
Options:
	-genome	FILE	The genome in FASTA, with the sequence IDs of the list.
	-list	FILE	An LTR_retriever pass list. Its strand column is the protein-domain strand.
	-update		Write the final strands into the list (and -scn), as LTR.identifier.pl writes them.
	-scn	FILE	The scn.adj file of the list, to update with -update.
	-trna	FILE	Mature tRNA sequences, with the 3' CCA, in FASTA. Default: database/tRNA.fa
	-threads	INT	Number of processes. Default: 1
	-h		Print this help and exit.
\n";

# LTRdigest defaults
my $PBS_RADIUS = 30;	# the PBS window: the last 30 nt of the 5' LTR and the first 31 internal nt
my $PBS_MAXEDIST = 1;
my ($PBS_ALILEN_MIN, $PBS_ALILEN_MAX) = (11, 30);
my $PBS_OFFSET_MAX = 5;	# start of the hit from the first internal nt
my $PBS_TOFF_MAX = 5;	# unpaired nt at the tRNA 3' end
my ($SW_MATCH, $SW_MISMATCH, $SW_INDEL) = (5, -10, -20);
my $PPT_RADIUS = 30;	# the PPT window: the last 31 internal nt and the first 30 nt of the 3' LTR
my ($PPT_LEN_MIN, $PPT_LEN_MAX) = (8, 30);
my $SEED = 11;	# = $PBS_ALILEN_MIN

# the calibrated call rule
my $PBS_CALL_TOFF_MAX = 2;
my ($PPT_CALL_GAP_MIN, $PPT_CALL_GAP_MAX) = (-2, 10);
my $CALL_MARGIN = 2;

my %CODE = (A => 0, C => 1, G => 2, T => 3, U => 3);

sub revcomp {
	my $s = reverse $_[0];
	$s =~ tr/ACGTUacgtuRYKMBVDHrykmbvdh/TGCAAtgcaaYRMKVBHDyrmkvbhd/;
	return $s;
}

## the tRNA library: [{name, len, rc (reverse complement), code (codes of rc)}], in file order
sub read_trna {
	my ($file) = @_;
	open my $fh, '<', $file or die "ERROR: Cannot read the tRNA file $file: $!\n";
	my (@raw, @trna);
	while (<$fh>) {
		s/\r?\n$//;
		if (/^>\s*(\S+)/) { push @raw, [$1, ''] }
		elsif (@raw) { s/\s+//g; $raw[-1][1] .= $_ }
	}
	close $fh;
	my $skip = 0;
	foreach (@raw) {
		my ($name, $seq) = @$_;
		$seq = uc $seq;
		$seq =~ tr/U/T/;
		if ($seq =~ /[^ACGT]/ or length $seq < $SEED + $PBS_TOFF_MAX + 1) { $skip++; next }
		my $rc = revcomp($seq);
		push @trna, {name => $name, len => length $seq, rc => $rc, code => [map { $CODE{$_} } split //, $rc]};
	}
	warn "Warning: $skip tRNA sequences with non-ACGT characters or shorter than 17 nt were skipped.\n" if $skip;
	die "ERROR: No tRNA sequence was read from $file\n" unless @trna;
	return \@trna;
}

## PBS seeds: every 11-mer that can open a PBS hit (>= 11 nt, <= 1 edit, at tRNA offset 0-5) -> tRNA index * 8 + offset.
## A local alignment starts with a match, so the edit is at position 1-10 of the first 11 nt of the hit.
sub pbs_index {
	my ($trna) = @_;
	my %idx;
	foreach my $t (0..$#$trna) {
		my $v = $trna->[$t]{rc};
		foreach my $vs (0..$PBS_TOFF_MAX) {
			my $exact = substr $v, $vs, $SEED;
			my $short = substr $v, $vs, $SEED - 1;
			my $long = substr $v, $vs, $SEED + 1;
			my %keys = ($exact => 1);
			foreach my $p (1..$SEED - 1) {
				$keys{substr($long, 0, $p) . substr($long, $p + 1)} = 1;	# a tRNA nt unpaired
				foreach my $b (qw(A C G T)) {
					my $k = $exact;
					substr($k, $p, 1) = $b;
					$keys{$k} = 1;	# a mismatch
					$keys{substr($short, 0, $p) . $b . substr($short, $p)} = 1;	# an extra nt in the element
					}
				}
			push @{$idx{$_}}, $t * 8 + $vs foreach keys %keys;
			}
		}
	return \%idx;
}

## the seeds of PBS hits in the window: tRNA index -> [[window position, tRNA offset], ...]
sub pbs_seeds {
	my ($w, $idx) = @_;
	my %seeds;
	foreach my $us ($PBS_RADIUS - $PBS_OFFSET_MAX .. $PBS_RADIUS + $PBS_OFFSET_MAX) {
		my $k = substr $w, $us, $SEED;
		next unless length $k == $SEED;
		my @keys = ($k);
		if ($k =~ /[^ACGT]/) {	# one ambiguous nt can only be the edit
			my $p = $-[0];
			next if substr($k, $p + 1) =~ /[^ACGT]/;
			@keys = map { my $x = $k; substr($x, $p, 1) = $_; $x } qw(A C G T);
			}
		foreach my $key (@keys) { push @{$seeds{$_ >> 3}}, [$us, $_ & 7] foreach @{$idx->{$key} || []} }
		}
	return \%seeds;
}

## the tRNAs that can give a PBS hit in the window, in library order
sub pbs_candidate_trnas { return sort { $a <=> $b } keys %{pbs_seeds(@_)} }

## the number of identical nt of $x from $i on and $y from $j on
sub match_len {
	my ($x, $i, $y, $j) = @_;
	return 0 if $i >= length $x or $j >= length $y;
	my $d = substr($x, $i) ^ substr($y, $j);
	return $d =~ /[^\0]/ ? $-[0] : length $d;
}

## an upper bound of the LTRdigest score of tRNA $T in window $w: a hit starts at one of its seeds and is the
## identical run from there (no edit), or that run plus one edit (at the first difference) and the next identical run
sub pbs_bound {
	my ($w, $T, $seeds) = @_;
	my ($v, $ub) = ($T->{rc}, 0);
	foreach (@$seeds) {
		my ($us, $vs) = @$_;
		my $l0 = match_len($w, $us, $v, $vs) or next;
		my @ext = (1 + match_len($w, $us + $l0 + 1, $v, $vs + $l0 + 1), 1 + match_len($w, $us + $l0 + 1, $v, $vs + $l0),
			match_len($w, $us + $l0, $v, $vs + $l0 + 1));	# mismatch, extra element nt, extra tRNA nt
		my $l1 = $l0 + (sort { $b <=> $a } @ext)[0];
		my $cap = length($w) - $us < $PBS_ALILEN_MAX ? length($w) - $us : $PBS_ALILEN_MAX;
		my $f = ($T->{len} - $vs) / $T->{len};
		my $off = abs($PBS_RADIUS - $us) || 1;
		my $b0 = $l0 >= $PBS_ALILEN_MIN ? ($l0 < $cap ? $l0 : $cap) * $f : 0;
		my $b1 = $l1 >= $PBS_ALILEN_MIN ? (($l1 < $cap ? $l1 : $cap) * $f) / $off : 0;
		$ub = $b0 if $b0 > $ub;
		$ub = $b1 if $b1 > $ub;
		}
	return $ub;
}

## Smith-Waterman as gt_swalign (GenomeTools extended/swalign.c), u = window, v = reverse-complemented tRNA, in two
## phases. A hit ends in tRNA column <= 36 (offset <= 5 + 30 nt + 1 extra tRNA nt), so phase 1 fills columns 1-36
## and traces back from their first maximal cell in v-major order (preferring replacement, then deletion = an extra
## u nt, then insertion = an extra v nt). Returns [] if that alignment is not a PBS hit (nor is the full alignment then),
## else [score, u start, u end, v start, edits, column 36] (0-based). $prof->[c][i] = score of u[i] against base c.
my $PBS_PHASE1 = $PBS_TOFF_MAX + $PBS_ALILEN_MAX + $PBS_MAXEDIST;
sub pbs_phase1 {
	my ($prof, $n, $T, $J, $w) = @_;
	my $v = $T->{code};
	my @H = ([(0) x ($n + 1)]);
	my ($best, $bi, $bj) = (-1, 0, 0);
	foreach my $j (1..$J) {
		my $pc = $prof->[$v->[$j - 1]];
		my $p = $H[$j - 1];
		my @c = (0);
		my $up = 0;
		foreach my $i (1..$n) {
			my $s = $p->[$i - 1] + $pc->[$i - 1];
			my $t = $up + $SW_INDEL;
			$s = $t if $t > $s;
			$t = $p->[$i] + $SW_INDEL;
			$s = $t if $t > $s;
			$s = 0 if $s < 0;
			$c[$i] = $up = $s;
			($best, $bi, $bj) = ($s, $i, $j) if $s > $best;
			}
		$H[$j] = \@c;
		}
	return [] unless $best > 0;
	my ($i, $j, $ed, $si, $sj) = ($bi, $bj, 0);
	while ($H[$j][$i] > 0) {
		($si, $sj) = ($i, $j);
		my $h = $H[$j][$i];
		if ($h == $H[$j - 1][$i - 1] + $prof->[$v->[$j - 1]][$i - 1]) {
			$ed++ if substr($w, $i - 1, 1) ne substr($T->{rc}, $j - 1, 1);
			$i--; $j--;
		} elsif ($h == $H[$j][$i - 1] + $SW_INDEL) {
			$ed++; $i--;
		} else {
			$ed++; $j--;
			}
		}
	my ($us, $ue, $vs) = ($si - 1, $bi - 1, $sj - 1);
	my $len = $ue - $us + 1;
	return [] unless $ed <= $PBS_MAXEDIST and abs($PBS_RADIUS - $us) <= $PBS_OFFSET_MAX and $len >= $PBS_ALILEN_MIN
		and $len <= $PBS_ALILEN_MAX and $vs <= $PBS_TOFF_MAX;
	return [$best, $us, $ue, $vs, $ed, $H[$J]];
}

## phase 2: whether a cell of tRNA columns > $J scores above $M (then the full alignment ends there: no hit)
sub sw_exceeds {
	my ($prof, $n, $v, $p, $J, $M) = @_;
	foreach my $j ($J + 1 .. scalar @$v) {
		my $pc = $prof->[$v->[$j - 1]];
		my @c = (0);
		my $up = 0;
		foreach my $i (1..$n) {
			my $s = $p->[$i - 1] + $pc->[$i - 1];
			my $t = $up + $SW_INDEL;
			$s = $t if $t > $s;
			$t = $p->[$i] + $SW_INDEL;
			$s = $t if $t > $s;
			$s = 0 if $s < 0;
			return 1 if $s > $M;
			$c[$i] = $up = $s;
			}
		$p = \@c;
		}
	return 0;
}

## the best PBS hit in a 61-nt window (index 30 = the first internal nt), as gt_pbs_find + the ranking of
## pbs_attach_results_to_gff3: highest score, ties to the tRNA first in the library. Candidates are aligned in the
## order of their score bound, until no bound can beat (or tie with an earlier tRNA) the best hit.
my %PBS_CACHE;
sub pbs_window_hit {
	my ($w, $trna, $idx) = @_;
	return $PBS_CACHE{$w} if exists $PBS_CACHE{$w};
	my $seeds = pbs_seeds($w, $idx);
	my @order = sort { $b->[1] <=> $a->[1] or $a->[0] <=> $b->[0] }
		grep { $_->[1] > 0 } map { [$_, pbs_bound($w, $trna->[$_], $seeds->{$_})] } keys %$seeds;
	my ($best, @prof, %phase1);
	if (@order) {
		my @u = map { exists $CODE{$_} ? $CODE{$_} : 4 } split //, $w;
		foreach my $c (0..3) { $prof[$c] = [map { $_ == $c ? $SW_MATCH : $SW_MISMATCH } @u] }
		}
	my $n = length $w;
	foreach (@order) {
		my ($t, $ub) = @$_;
		last if $best and ($ub < $best->{score} or ($ub == $best->{score} and $t > $best->{t}));
		my $T = $trna->[$t];
		my $J = $T->{len} < $PBS_PHASE1 ? $T->{len} : $PBS_PHASE1;
		my $r = $phase1{substr($T->{rc}, 0, $J)} ||= pbs_phase1(\@prof, $n, $T, $J, $w);
		next unless @$r;
		my ($M, $us, $ue, $vs, $ed, $col) = @$r;
		my ($len, $off) = ($ue - $us + 1, abs($PBS_RADIUS - $us));
		my $pen = ($ed == 0 or $off == 0) ? 1 : $ed * $off;
		my $score = ($len * (($T->{len} - $vs) / $T->{len})) / $pen;
		next if $best and ($score < $best->{score} or ($score == $best->{score} and $t > $best->{t}));
		next if sw_exceeds(\@prof, $n, $T->{code}, $col, $J, $M);
		$best = {t => $t, us => $us, ue => $ue, toff => $vs, edist => $ed, len => $len, score => $score, trna => $T->{name}};
		}
	return $PBS_CACHE{$w} = $best;
}

## LTRdigest's PPT HMM (ltr/ltrdigest_ppt_visitor.c) in GenomeTools' log space (extended/hmm.c)
my ($IN, $OUT, $UBOX, $NST) = (0, 1, 2, 3);
my $LOG0 = -99999.0;
my $DBL_EPSILON = 2 ** -52;
my (@INIT, @TRANS, @EMIT);
{
	my $lg = sub { $_[0] == 0 ? $LOG0 : log $_[0] };
	my $ex = sub { $_[0] == $LOG0 ? 0 : exp $_[0] };
	@INIT = map { $lg->(1 / 4) } 0..3;
	foreach my $f (0..3) { foreach my $t (0..3) { $TRANS[$f][$t] = $f == $t ? $LOG0 : $lg->(0.05) } }
	foreach my $f (0..3) {	# gt_hmm_set_missing_transition_probabilities
		my ($sum, $miss) = (0, 0);
		foreach my $t (0..3) { my $p = $ex->($TRANS[$f][$t]); if ($p == 0) { $miss++ } else { $sum += $p } }
		foreach my $t (0..3) { $TRANS[$f][$t] = $lg->((1 - $sum) / $miss) if $ex->($TRANS[$f][$t]) == 0 }
		}
	my $nonu = (1 - 0.91) / 3;
	my %e = ($IN => [0.97 / 2, 0.03 / 2, 0.97 / 2, 0.03 / 2, 0], $OUT => [0.25, 0.25, 0.25, 0.25, 0],
		$UBOX => [$nonu, $nonu, $nonu, 0.91, 0], $NST => [0, 0, 0, 0, 1]);	# A C G T other
	foreach my $s (0..3) { $EMIT[$s] = [map { $lg->($_) } @{$e{$s}}] }
}

## gt_hmm_decode: Viterbi with ties to the lower state; the end state is the last one above state 0
sub viterbi {
	my ($w) = @_;
	my @x = map { exists $CODE{$_} ? $CODE{$_} : 4 } split //, $w;
	my @M = map { $INIT[$_] + $EMIT[$_][$x[0]] } 0..3;
	my @bt;
	foreach my $col (1..$#x) {
		my (@N, @b);
		foreach my $row (0..3) {
			my $em = $EMIT[$row][$x[$col]];
			my $best = $M[0] + $TRANS[0][$row] + $em;
			my $arg = 0;
			foreach my $p (1..3) {
				my $t = $M[$p] + $TRANS[$p][$row] + $em;
				if ($t - $best > $DBL_EPSILON) { $best = $t; $arg = $p }
				}
			($N[$row], $b[$row]) = ($best, $arg);
			}
		@M = @N;
		$bt[$col] = \@b;
		}
	my $s = 0;
	foreach my $row (1..3) { $s = $row if $M[$row] - $M[0] > $DBL_EPSILON }
	my @path = ($s);
	for (my $col = $#x; $col > 0; $col--) { unshift @path, $bt[$col][$path[0]] }
	return @path;
}

## the best PPT in a window of 2r+1 nt (index r = the last internal nt), as gt_group_hits + gt_ppt_score:
## tracts over indices 0..2r-2, the one ending closest to index r (ties to the first)
my %PPT_CACHE;
sub ppt_window_hit {
	my ($w, $r) = @_;
	return $PPT_CACHE{$w} if exists $PPT_CACHE{$w};
	my @d = viterbi($w);
	my ($best, $bd);
	my $start = 0;
	foreach my $i (0..2 * $r - 2) {
		next unless $i + 2 == 2 * $r or $d[$i + 1] != $d[$i];
		if ($d[$i] == $IN) {
			my $len = $i - $start + 1;
			if ($len >= $PPT_LEN_MIN and $len <= $PPT_LEN_MAX) {
				my $dist = abs($r - $i);
				($best, $bd) = ([$start, $i], $dist) if !defined $bd or $dist < $bd;
				}
			}
		$start = $i + 1;
		}
	return $PPT_CACHE{$w} = $best;
}

## the 61-nt regions that hold all search windows of an element (5' LTR first, + strand):
## A = last 30 nt of the 5' LTR + first 31 internal nt; B = last 31 internal nt + first 30 nt of the 3' LTR
sub regions {
	my ($S, $ll, $rl) = @_;
	my $L = length $S;
	my $A = ($ll >= 30 and $ll + 31 <= $L) ? substr($S, $ll - 30, 61) : undef;
	my $B = ($L - $rl - 31 >= 0 and $rl >= 30) ? substr($S, $L - $rl - 31, 61) : undef;
	return ($A, $B);
}

## PBS and PPT of both strand hypotheses, in 1-based coordinates on the + strand of the element (1..L).
## gap: nt between the LTR and the PBS start (PPT end); negative = inside the LTR.
sub evidence {
	my ($A, $B, $ll, $rl, $L, $trna, $idx) = @_;
	$A = uc $A if defined $A;
	$B = uc $B if defined $B;
	my %ev = (plus => {}, minus => {});
	if (defined $A) {	# + PBS
		my $h = pbs_window_hit($A, $trna, $idx);
		$ev{plus}{pbs} = {%$h, start => $ll - 29 + $h->{us}, end => $ll - 29 + $h->{ue}, gap => $h->{us} - $PBS_RADIUS} if $h;
		}
	if (defined $B) {	# - PBS
		my $h = pbs_window_hit(revcomp($B), $trna, $idx);
		$ev{minus}{pbs} = {%$h, start => $L - $rl + 30 - $h->{ue}, end => $L - $rl + 30 - $h->{us}, gap => $h->{us} - $PBS_RADIUS} if $h;
		}
	my $r = $rl - 1 < $PPT_RADIUS ? $rl - 1 : $PPT_RADIUS;
	if (defined $B and $r > 0) {	# + PPT
		my $h = ppt_window_hit(substr($B, $PPT_RADIUS - $r, 2 * $r + 1), $r);
		$ev{plus}{ppt} = {start => $L - $rl - $r + $h->[0], end => $L - $rl - $r + $h->[1], len => $h->[1] - $h->[0] + 1, gap => $r - $h->[1]} if $h;
		}
	$r = $ll - 1 < $PPT_RADIUS ? $ll - 1 : $PPT_RADIUS;
	if (defined $A and $r > 0) {	# - PPT
		my $h = ppt_window_hit(revcomp(substr($A, $PPT_RADIUS - $r, 2 * $r + 1)), $r);
		$ev{minus}{ppt} = {start => $ll + $r + 1 - $h->[1], end => $ll + $r + 1 - $h->[0], len => $h->[1] - $h->[0] + 1, gap => $r - $h->[1]} if $h;
		}
	return \%ev;
}

## one hypothesis: the element $S (5' LTR first) on its + strand
sub scan {
	my ($S, $ll, $rl, $trna, $idx) = @_;
	return evidence(regions($S, $ll, $rl), $ll, $rl, length $S, $trna, $idx)->{plus};
}

## both hypotheses of the element $S
sub element_evidence {
	my ($S, $ll, $rl, $trna, $idx) = @_;
	return evidence(regions($S, $ll, $rl), $ll, $rl, length $S, $trna, $idx);
}

## the calibrated call
sub pbs_ok { my ($p) = @_; return ($p and $p->{edist} <= 1 and $p->{len} >= 11 and abs($p->{gap}) <= 5 and $p->{toff} <= $PBS_CALL_TOFF_MAX) ? 1 : 0 }
sub ppt_ok { my ($q) = @_; return ($q and $q->{len} >= 8 and $q->{gap} >= $PPT_CALL_GAP_MIN and $q->{gap} <= $PPT_CALL_GAP_MAX) ? 1 : 0 }
sub pbs_q { my ($p) = @_; return $p ? $p->{len} - 2 * $p->{edist} : 0 }
sub ppt_q { my ($q) = @_; return $q ? $q->{len} - 2 * abs($q->{gap}) : 0 }
sub decide {
	my ($okp, $okm, $qp, $qm) = @_;
	return '+' if $okp and !$okm;
	return '-' if $okm and !$okp;
	return '.' unless $okp and $okm;
	return '+' if $qp - $qm >= $CALL_MARGIN;
	return '-' if $qm - $qp >= $CALL_MARGIN;
	return '.';
}
sub call_strand {
	my ($ev) = @_;
	my ($p, $m) = @$ev{qw(plus minus)};
	my $pbs = decide(pbs_ok($p->{pbs}), pbs_ok($m->{pbs}), pbs_q($p->{pbs}), pbs_q($m->{pbs}));
	my $ppt = decide(ppt_ok($p->{ppt}), ppt_ok($m->{ppt}), ppt_q($p->{ppt}), ppt_q($m->{ppt}));
	my ($call, $basis) = ('.', 'none');
	if ($pbs ne '.' and $ppt ne '.') { ($call, $basis) = $pbs eq $ppt ? ($pbs, 'PBS+PPT') : ('.', 'PBS/PPT conflict') }
	elsif ($pbs ne '.') { ($call, $basis) = ($pbs, 'PBS only') }
	elsif ($ppt ne '.') { ($call, $basis) = ($ppt, 'PPT only') }
	return {pbs => $pbs, ppt => $ppt, call => $call, basis => $basis};
}

## final strand from the domain strand and the call: PBS/PPT replace the domain strand only when they agree.
## Returns (strand, evidence = signals supporting it or 'none', conflict = signals against it; for '?', PBS and PPT
## if they disagree).
sub final_strand {
	my ($dom, $c) = @_;
	my ($pbs, $ppt) = @$c{qw(pbs ppt)};
	if ($dom ne '+' and $dom ne '-') {
		return ('?', 'none', ($pbs ne '.' and $ppt ne '.') ? 'PBS+PPT' : '') if $c->{call} eq '.';
		return ($c->{call}, join('+', ($pbs eq $c->{call} ? 'PBS' : ()), ($ppt eq $c->{call} ? 'PPT' : ())), '');
		}
	return ($pbs, 'PBS+PPT', 'domain') if $pbs ne '.' and $pbs eq $ppt and $pbs ne $dom;
	my $opp = $dom eq '+' ? '-' : '+';
	my @ev = ('domain', ($pbs eq $dom ? 'PBS' : ()), ($ppt eq $dom ? 'PPT' : ()));
	my @cf = (($pbs eq $opp ? 'PBS' : ()), ($ppt eq $opp ? 'PPT' : ()));
	return ($dom, join('+', @ev), join('+', @cf));
}

my @COLS = qw(LTR_loc domain_strand strand strand_evidence strand_conflict call_pbs call_ppt call_basis);
foreach my $h (qw(plus minus)) {
	push @COLS, map { "pbs_${h}_$_" } qw(start end len edist gap trna_offset trna score);
	push @COLS, map { "ppt_${h}_$_" } qw(start end len gap);
	}

## the final strand written back as LTR.identifier.pl writes elements: a '-' element runs from its 3' to its 5' end
## (chr:end..start) and has its motif and TSD on the minus strand; '?' and '+' read on the plus strand
sub flips { my ($old, $new) = @_; return (($old eq '-') xor ($new eq '-')) ? 1 : 0 }
sub rc_field { my ($s) = @_; return $s =~ /^NA$/i ? $s : revcomp($s) }
sub update_list_line {	# a pass list line (tab-separated): loc, category, motif:, TSD:, ..., strand (column 9)
	my ($line, $new) = @_;
	my @f = split /\t/, $line, -1;
	return $line if $f[8] eq $new;
	if (flips($f[8], $new)) {
		$f[0] =~ s/^(.*):([0-9]+)\.\.([0-9]+)$/$1:$3..$2/;
		$f[2] =~ s/^motif:(.*)$/'motif:' . rc_field($1)/e;
		$f[3] =~ s/^TSD:(.*)$/'TSD:' . rc_field($1)/e;
		}
	$f[8] = $new;
	return join "\t", @f;
}
sub scn_key { my @f = split ' ', $_[0]; return defined $f[11] ? "$f[11]:$f[0]..$f[1]" : '' }
sub update_scn_line {	# an scn.adj line: start end ... chr(12) direction(13) TSD(14) lTSD rTSD motif(17) ...; spacing kept
	my ($line, $new) = @_;
	my @t = split /(\s+)/, $line, -1;
	my $o = $t[0] eq '' ? 2 : 0;	# fields at even positions, after a leading blank if any
	return $line if $t[$o + 24] eq $new;
	if (flips($t[$o + 24], $new)) {
		$t[$o + 26] = rc_field($t[$o + 26]);
		$t[$o + 32] = rc_field($t[$o + 32]);
		}
	$t[$o + 24] = $new;
	return join '', @t;
}

## write the final strands into the list and the scn file; through temporary files, so a failure leaves them unchanged
sub update_files {
	my ($list, $scn, $el, $strand) = @_;
	my %scn_new;
	my $k = 0;
	open my $in, '<', $list or die "ERROR: Cannot read the list $list: $!\n";
	open my $out, '>', "$list.pbs_ppt.tmp" or die "ERROR: Cannot write $list.pbs_ppt.tmp: $!\n";
	while (my $line = <$in>) {
		if ($line =~ /^#/ or $line =~ /^\s*$/) { print $out $line; next }
		chomp $line;
		my ($dom, $chr, $s, $e) = @{$el->[$k]}[1, 2, 3, 6];
		my $new = $strand->[$k++];
		$scn_new{"$chr:$s..$e"} = $new if $new ne $dom;
		print $out update_list_line($line, $new), "\n";
		}
	close $in;
	close $out or die "ERROR: Cannot write $list.pbs_ppt.tmp: $!\n";
	die "ERROR: The list $list changed while it was read\n" unless $k == @$strand;
	if ($scn) {
		open $in, '<', $scn or die "ERROR: Cannot read $scn: $!\n";
		open $out, '>', "$scn.pbs_ppt.tmp" or die "ERROR: Cannot write $scn.pbs_ppt.tmp: $!\n";
		while (my $line = <$in>) {
			if ($line !~ /^#/ and $line =~ /\S/) {
				chomp $line;
				my $key = scn_key($line);
				$line = (exists $scn_new{$key} ? update_scn_line($line, $scn_new{$key}) : $line) . "\n";
				}
			print $out $line;
			}
		close $in;
		close $out or die "ERROR: Cannot write $scn.pbs_ppt.tmp: $!\n";
		rename "$scn.pbs_ppt.tmp", $scn or die "ERROR: Cannot update $scn: $!\n";
		}
	rename "$list.pbs_ppt.tmp", $list or die "ERROR: Cannot update $list: $!\n";
}

## one output line; $off = genome position of the element start - 1; with $update, the location of the updated list
sub table_line {
	my ($loc, $dom, $ev, $off, $update) = @_;
	my $c = call_strand($ev);
	my @f = ($loc, $dom, final_strand($dom, $c), @$c{qw(pbs ppt basis)});
	$f[0] =~ s/^(.*):([0-9]+)\.\.([0-9]+)$/$1:$3..$2/ if $update and flips($dom, $f[2]);
	$f[4] = 'NA' if $f[4] eq '';
	foreach my $h (qw(plus minus)) {
		my ($p, $q) = @{$ev->{$h}}{qw(pbs ppt)};
		push @f, $p ? ($p->{start} + $off, $p->{end} + $off, @$p{qw(len edist gap toff trna)}, sprintf('%.3f', $p->{score})) : ('NA') x 8;
		push @f, $q ? ($q->{start} + $off, $q->{end} + $off, @$q{qw(len gap)}) : ('NA') x 4;
		}
	return join("\t", @f) . "\n";
}

sub main {
	my ($genome, $list, $scn, $update, $trna_file, $threads) = ('', '', '', 0, dirname(__FILE__) . '/../database/tRNA.fa', 1);
	while (@ARGV) {
		my $o = shift @ARGV;
		if ($o =~ /^-genome$/i) { $genome = shift @ARGV }
		elsif ($o =~ /^-list$/i) { $list = shift @ARGV }
		elsif ($o =~ /^-scn$/i) { $scn = shift @ARGV }
		elsif ($o =~ /^-update$/i) { $update = 1 }
		elsif ($o =~ /^-trna$/i) { $trna_file = shift @ARGV }
		elsif ($o =~ /^-threads$/i or $o =~ /^-t$/i) { $threads = shift @ARGV }
		elsif ($o =~ /^-h(elp)?$/i) { print $usage; exit 0 }
		else { die "ERROR: Unknown option $o\n$usage" }
		}
	die "ERROR: -genome and -list are required\n$usage" unless $genome and $list;
	die "ERROR: -threads must be a positive integer\n" unless defined $threads and $threads =~ /^[1-9][0-9]*$/;

	# elements of the list, in list order: [loc, domain strand, seqid, start, 5' LTR end, 3' LTR start, end]
	open my $lfh, '<', $list or die "ERROR: Cannot read the list $list: $!\n";
	my (@el, %by_chr);
	while (<$lfh>) {
		next if /^#/ or /^\s*$/;
		my ($loc, $in, $strand) = (split)[0, 6, 8];
		my ($chr, $s, $e) = $loc =~ /^(.*):([0-9]+)\.\.([0-9]+)$/ or die "ERROR: Cannot read the element location of $list line $.: $_";
		($s, $e) = ($e, $s) if $s > $e;
		my ($is, $ie) = $in =~ /^IN:([0-9]+)\.\.([0-9]+)$/ or die "ERROR: Cannot read the internal region of $list line $.: $_";
		push @el, [$loc, $strand, $chr, $s, $is - 1, $ie + 1, $e];
		push @{$by_chr{$chr}}, $#el;
		}
	close $lfh;

	# cut the two regions of each element out of the genome
	my @work;	# [loc, domain strand, A, B, ll, rl, L, start - 1]
	open my $gfh, '<', $genome or die "ERROR: Cannot read the genome $genome: $!\n";
	{
		local $/ = "\n>";
		while (my $rec = <$gfh>) {
			$rec =~ s/^>//;
			$rec =~ s/>$//;
			my ($id, $seq) = split /\n/, $rec, 2;
			$id = (split /\s+/, $id)[0];
			next unless exists $by_chr{$id};
			$seq =~ s/\s+//g;
			foreach my $k (@{$by_chr{$id}}) {
				my ($loc, $strand, $chr, $s, $le, $rs, $e) = @{$el[$k]};
				die "ERROR: $loc is beyond the end of $id\n" if $e > length $seq;
				my ($ll, $rl, $L) = ($le - $s + 1, $e - $rs + 1, $e - $s + 1);
				my $A = ($ll >= 30 and $ll + 31 <= $L) ? substr($seq, $le - 30, 61) : undef;
				my $B = ($L - $rl - 31 >= 0 and $rl >= 30) ? substr($seq, $rs - 32, 61) : undef;
				$work[$k] = [$loc, $strand, $A, $B, $ll, $rl, $L, $s - 1];
				}
			delete $by_chr{$id};
			}
	}
	close $gfh;
	die "ERROR: Sequence $_ of the list is not in the genome $genome\n" foreach sort keys %by_chr;

	my $trna = read_trna($trna_file);
	my $idx = pbs_index($trna);
	my $run = sub { join '', map { my $w = $_; table_line($w->[0], $w->[1], evidence(@$w[2..6], $trna, $idx), $w->[7], $update) } @_ };

	my $table = '';
	if ($threads == 1 or @work < 2 * $threads) {
		$table = $run->(@work);
		} else {
		$table = run_parallel($run, \@work, $threads);
		}
	print '#' . join("\t", @COLS) . "\n", $table;
	update_files($list, $scn, \@el, [map { (split /\t/)[2] } split /\n/, $table]) if $update;
}

## the table lines of @$work from $threads processes, in order
sub run_parallel {
	my ($run, $work, $threads) = @_;
	my @work = @$work;
	my $tmp = tempdir("find_PBS_PPT.XXXXXX", DIR => '.', CLEANUP => 1);
	STDOUT->flush;
	my $size = ceil(@work / $threads);
	my @kids;
	for (my $i = 0; $i < @work; $i += $size) {
		my $j = $i + $size - 1;
		$j = $#work if $j > $#work;
		my $part = "$tmp/" . scalar @kids;
		my $pid = fork;
		die "ERROR: fork failed: $!\n" unless defined $pid;
		if ($pid == 0) {	# never return into the parent's loop, never run its END blocks
			my $ok = eval {
				open my $out, '>', $part or die "ERROR: Cannot write $part: $!\n";
				print $out $run->(@work[$i..$j]);
				close $out or die "ERROR: Cannot write $part: $!\n";
				1;
				};
			print STDERR $@ unless $ok;
			_exit($ok ? 0 : 1);
			}
		push @kids, [$pid, $part];
		}
	my $table = '';
	foreach (@kids) {
		my ($pid, $part) = @$_;
		waitpid($pid, 0);
		die "ERROR: a find_PBS_PPT.pl process failed\n" if $?;
		open my $in, '<', $part or die "ERROR: Cannot read $part: $!\n";
		$table .= do { local $/; <$in> };
		close $in;
		}
	return $table;
}

main() unless caller;
1;
