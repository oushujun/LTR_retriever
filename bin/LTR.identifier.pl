#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use File::Basename;
use File::Spec;
use File::Temp qw(tempdir);

my $usage="
	Usage: perl LTR.identifier.pl LTR_index
	Example: perl LTR.identifier.pl Chr1-TIGR7-sd50sm90v8mx20
	\n";
my $version="
LTR.identifier.pl
LTR.identifier: Alignment assisted examination of LTR candidates
Author: Shujun Ou (shujun.ou.1\@gmail.com), Department of Horticulture, Michigan State University, East Lansing, MI, 48823, USA
Version:
	4.8 Add the pdist and K2P model 2025/03/23
	4.7 Incorporate TEsorter results 2023/05/01
	4.6 Improvement: only consider SNPs for age estimation 2019/01/25
	4.5 Improve TSD-motif identification 2018/12/08
	4.0 Use Thread::Queue instead of Semaphore for multi-threading 2018/04/01
	3.6 Enable multi-threading 2016/6/30
	3.5 Add TE family annotation and age estimation to scn output 2016/6/16
	3.0 Combine Boundary_correction.pl and LTR.identifier.pl 2016/4/29
	2.0 Improve TSD reporting, take multiple motifs into account, default no TSD control. 2016/4/27
	1.6 Improve TSD reporting. 2015/4/9
	1.5 Improve TSD reporting. 2015/2/24
	1.0 2014/02/14
\n";

my $name=$ARGV[0];
die "ERROR: $usage" unless defined $name;
my $a_cutoff=0.6; #0.5 #minimum ac, bd alignment portion
my $s_cutoff=60; #alignment similariy (%) cutoff
my $w_size=7; #word size for boundary alignment
my $TSD_ctrl=0; #1 requires a TSD to be an authentic LTR, 0 does not require
my $boundary_ctrl=1; #1 for boundaries alignment, 0 for no alignment
my $boundary_N=25; #if any boundary has more than 25 bp missing (nN-), report as false positive
my $length_diff=15;     #boundary adjustments for length difference between adjusted LTSs higher than this value will be discarted.
my $minlen=100; #dft=100, minmum LTR region length for nmtf candidates
my $model="K2P"; #distance model, default K2P, available: K2P, JC69, pdist
my $miu="1.3e-8"; #neutral mutation rate, default: 1.3e-8 (rice) per bp per ya
my @motif=qw/TGCT TACA TACT TGGA TATA TGTA TGCA/;
my $threads="4"; #threads to run this program
my $blastplus=''; #path to the blast+ directory
my $timeout = 120; # timeout blastn after 120s.

#obtain the exact path for the program location
my $script_path = dirname(__FILE__);

my $List;
my $FA;
my $ANNO;

my $k=0;
my $argv='';
foreach (@ARGV){
	$argv.="$_ ";
        $List=$ARGV[$k+1] if /^-list$/i;
        $FA=$ARGV[$k+1] if /^-seq$/i;
	$ANNO=$ARGV[$k+1] if /^-anno$/i;
        $length_diff=$ARGV[$k+1] if /^-lendiff$/i;
	$minlen=$ARGV[$k+1] if /^-minlen$/i;
	$s_cutoff=$ARGV[$k+1] if /^-flanksim$/i;
	$a_cutoff=$ARGV[$k+1] if /^-flankaln$/i;
	$boundary_N=$ARGV[$k+1] if /^-flankmiss$/i;
	$boundary_ctrl=0 if /^-b$/i;
	$TSD_ctrl=1 if /^-tsdaln$/i;
	$model=uc $ARGV[$k+1] if /^-m$|^-model$/i;
	$miu=$ARGV[$k+1] if /^-u$/i;
	@motif=(split /\s+/, $1) if $argv=~/-motif\s+\[([atcgnx ]+)\]/i;
	$threads=$ARGV[$k+1] if /^-t$|^-threads$/i;
	$blastplus=$ARGV[$k+1] if /^-blastplus$/i;
	die $version if /^-v$/i;
	$k++;
	}
$a_cutoff-=0.10 if $boundary_ctrl==0;

open List, "<$List" or die "ERROR: No candidate list file!\n$usage";
open FA, "<$FA" or die "ERROR: No candidate sequence file!\n$usage";
die "Error: You specified -model $model, but only K2P, JC69, or pdist models are supported!\n" unless $model =~ /K2P|JC69|PDIST/i;

##Store LTR information in hash
my %scn :shared;
my $head='';
while (<List>){
	next if /^\s+$/;
	s/^\s+//;
	if (/^#/){
		$head.=$_;
		next;
		}
	my ($start, $end, $len, $ls, $le, $ll, $rs, $re, $rl, $sim, $id)=split;
	$scn{"$start..$end"}=shared_clone([split]); #store in %scn
	$scn{"$start..$end"}[9]*=0.01 if $scn{"$start..$end"}[9] ne "NA"; #convert % to decimal
	}
close List;

##protein superfamily and strand annotation
my %anno;
if (defined $ANNO){
	# listing orders that belong to LTRs and not LTRs
	my @notLTR = ("mixture", "notLTR", "TIR", "Helitron", "LINE", "SINE", "Maverick", "pararetrovirus");
	my @yesLTR = ("LTR", "DIRS");
	open ANNO, "<$ANNO" or die "ERROR: Can't read the .anno file!\n";
	while (<ANNO>){
		next if /^#/;
		s/^\s+//;
		my ($id, $order, $superfamily, $strand)=(split)[0,1,2,3];
		#id is like: Chr1:106472..118130|Chr1:106522..118080
		$id=~s/.*\|(?:\S+:)?([0-9]+\.\.[0-9]+)$/$1/;
		next unless defined $scn{$id};

		# convert order to two categories
		$order = "LTR" if grep { $_ eq $order } @yesLTR;
		$order = "notLTR" if grep { $_ eq $order } @notLTR;

		# assign values to the big table
		if (defined $scn{$id}[12]){
			$scn{$id}[12] = $strand if $scn{$id}[12] eq '?';
			} else {
			$scn{$id}[12] = $strand;
			}
		if (defined $scn{$id}[17]){
			$scn{$id}[17] = "notLTR" if $order eq "notLTR";
			} else {
			$scn{$id}[17]=$order;
			}
		if (defined $scn{$id}[18]){
			$scn{$id}[18] = $superfamily if $scn{$id}[18] eq 'unknown';
			} else {
			$scn{$id}[18]=$superfamily;
			}
		}
	close ANNO;
	}

#Store sequence information into an ordered, shared candidate list (so worker threads in
#each wave can read it by index without cloning the sequence data).
my @candidates :shared;
$/ = "\n>";
while (<FA>){
	chomp;
	s/>//g;
	my ($name, $seq) = (split /\n/, $_, 2);
	next unless defined $seq;
	next if $seq eq '';
	$seq=~s/\s+//g;
	$seq=uc $seq;
	my $id;
	if ($name =~ /\|(?:\S+:)?([0-9]+\.\.[0-9]+)$/) {
		$id = $1;
	}
	next unless defined $id and defined $scn{$id};
	push @candidates, shared_clone([$name, $seq]);
	}
close FA;
$/ = "\n";

##open scn.adj file and print out the header
open SCN, ">$List.adj" or die "ERROR: $!";
print SCN "#LTR boundary fine-grain adjustment and annotation have been performed by LTR_retriever (Shujun Ou, oushujun\@msu.edu)\n$head";

##Batched, load-balanced blastn.
## Each candidate needs (1) a self-alignment of its LTR region (boundary correction +
## divergence) and (2) two flanking-region alignments. The flanking inputs depend on the
## coordinates that the self-alignment corrects, so the blasts form two dependent stages.
## The original spawned blastn (plus a perl helper and several shell subshells) inline, one
## candidate at a time. Each blastn is cheap on CPU (~0.1s) but latency-bound (~1s wall under
## load), so total runtime was dominated by per-call overhead and limited blast concurrency.
## Here we instead run the analysis in three parallel passes:
##   wave 1 - collect each candidate's self-alignment job;   then run them all as one batch
##   wave 2 - collect each candidate's flanking jobs;         then run them all as one batch
##   wave 3 - all blasts are cached, so analyze and emit
## Each pass is multi-threaded; the batches run blastn back-to-back at full concurrency with
## no perl/shell overhead in between. Blast commands and inputs are byte-identical to the
## per-candidate version, so the output is identical; only the orchestration changes.
my %BLAST_CACHE :shared;   #job key => shared array of that job's own blastn output lines
my %PENDING :shared;       #job key => shared [query_seq, subject_seq, blast_spec] to run
my %DEFALSE :shared;       #candidate index => defalse text (collected for ordered output)
my %STATE :shared;         #candidate index => boundary-corrected state carried from wave 2 to wave 3
my $COLLECT_STAGE :shared = 0;  #blast stage currently being collected (1=self, 2=flanking)
my $EMIT :shared = 0;           #1 only in the final pass: collect defalse + update %scn
#NOTE: no -parse_deflines. We tag each pooled sequence with its numeric job index as the defline.
#With -parse_deflines, BLAST reads a bare integer defline as a GI number and can mis-attribute hits
#among near-identical pooled sequences (subtle, breaks the flanking filter); without it the integer
#is used verbatim, so the qseqid==sseqid self-hit filtering is reliable. (qseqid/sseqid are not used
#downstream anyway - only the alignment columns are.)
my $self_spec  = "-dust no -outfmt '6 qseqid sseqid sstart send slen qstart qend qlen length nident btop'";
#Each pending job (a single-sequence self-alignment, or a flank-vs-flank pair) is too small to
#justify its own blastn process - the launch overhead dwarfs the alignment. So we pool jobs into
#chunks and run ONE blastn per chunk: all chunk queries vs all chunk subjects, then keep only the
#hits a job has against ITS OWN subject (qseqid==sseqid, both tagged with the job's index). This
#amortizes the per-process cost over a whole chunk while leaving each job's result identical to
#its standalone blast. Chunk size is bounded so within-chunk cross-hits stay cheap.
my $chunk_max = 500;       #ceiling on jobs/chunk; the actual size is chosen adaptively per batch
                           #(estimate_chunk_size) from the measured cross-hit rate
#Blast inputs are written to per-worker temp files on fast local storage (RAM-backed /dev/shm if
#available, else the system temp dir) rather than the -list directory, which may be slow/networked.
#tempdir(CLEANUP=>1) removes the directory automatically on exit (including on error).
my $tmpbase = (-d "/dev/shm" && -w "/dev/shm") ? "/dev/shm" : File::Spec->tmpdir;
my $tmpdir = tempdir("LTRid.XXXXXXXX", DIR => $tmpbase, CLEANUP => 1);
my $idx_queue;             #per-wave shared queue of candidate indices
#The flanking blast (60bp pairs, -evalue 1000) reports hits whose e-values depend on the blast
#database size. Pooling many pairs into one chunk enlarges that database and would drop borderline
#hits differently depending on chunk size (i.e. on -threads). To keep each pair's reported hits
#identical to its standalone blast regardless of chunking, we fix the e-value search space with
#-searchsp, set to the value a single 60bp-vs-60bp blast uses (probed once here so it stays correct
#across blast versions / scoring).
my $flank_searchsp = 0;
{
	#probe with ONE 60bp sequence as both query and subject, so the search space matches a
	#single flank-vs-flank blast (a 1-sequence, 60bp database)
	my $pf = "$tmpdir/searchsp_probe.fa";
	if (open my $ph, '>', $pf){
		print $ph ">a\n", ("ACGT" x 15), "\n"; close $ph;
		foreach (qx(${blastplus}blastn -query $pf -subject $pf -evalue 1000 -word_size $w_size -dust no 2> /dev/null)){
			$flank_searchsp = $1 if /Effective search space used:\s*(\d+)/;
			}
		unlink $pf;
		}
}
my $flank_spec = "-evalue 1000 -word_size $w_size -dust no -outfmt 6"
	. ($flank_searchsp > 0 ? " -searchsp $flank_searchsp" : "");

#Stage 1: collect + batch-run each candidate's self-alignment
$COLLECT_STAGE = 1; $EMIT = 0; &run_wave;
&flush_blast;
#Stage 2: collect + batch-run flanking alignments (need stage-1 corrected coordinates)
$COLLECT_STAGE = 2; $EMIT = 0; &run_wave;
&flush_blast;
#Final pass: all blasts cached; analyze and collect results
$COLLECT_STAGE = 0; $EMIT = 1; &run_wave;
#($tmpdir is auto-removed at exit by File::Temp CLEANUP)

##print defalse results to STDOUT in candidate (input) order - deterministic, thread-count independent
foreach my $i (sort {$a<=>$b} keys %DEFALSE){ print $DEFALSE{$i}; }

##print out entries that could not pass initial screening criteria to scn.adj
foreach my $key (sort{$a cmp $b}(keys %scn)){
	foreach (0..$#{$scn{$key}}){
		next unless defined $scn{$key};
		next unless defined $scn{$key}[$_];
		print SCN "$scn{$key}[$_]  ";
		}
	print SCN "\n";
	}
close SCN;


##LTR structural analysis, split so each candidate's pre-work is computed once. Identifier just
##dispatches by pass: process_pre does the parse / self-align / boundary-correct / build-flanking-
##windows part (registering blast jobs); process_post finishes from the cached blasts. The boundary
##correction therefore runs once (wave 2) and its result is reused in wave 3, not recomputed.
sub Identifier() {
	while (defined(my $cand_idx = $idx_queue->dequeue())){
		if ($EMIT){                              #wave 3: finish from cached blasts + saved state
			my $st = $STATE{$cand_idx};
			&process_post($cand_idx, $st) if defined $st;
			} else {                             #waves 1-2: collect blast jobs
			my $st = &process_pre($cand_idx);
			$STATE{$cand_idx} = $st if defined $st;  #a serialized string (one shared scalar), wave-2 only
			}
		}
}

##Waves 1-2. Wave 1 stops after registering the self-alignment job; wave 2 also does the boundary
##correction and registers the two flanking jobs, then returns the state process_post needs (or
##undef if this candidate is skipped).
sub process_pre {
	my ($cand_idx) = @_;
##Structural analysis, main program
	my ($name, $seq)=(@{$candidates[$cand_idx]}[0], @{$candidates[$cand_idx]}[1]);
	my $decision="raw"; #conclusion of whether the element is a LTR
	my ($chr, $seq_start, $seq_end, $ltr_start, $ltr_end);
	($chr, $seq_start, $seq_end, $ltr_start, $ltr_end)=($1, $2, $3, $4, $5) if $name=~/^(\S+):([0-9]+)\.\.([0-9]+)\|(?:\S+:)?([0-9]+)\.\.([0-9]+)/;  #eg: Chr4:10009589..10017157|10009609..10017137 or 10.dna.chromosome.ch:100016935..100026312|100016935..100026312
	my $id="$ltr_start..$ltr_end";
	return if $id eq '';
	my @info = @{$scn{$id}};

##Coarse boundary correction - after correction, coordinates may still have 1-2 bp shifted from the real case
	my $ltr=substr $seq, $ltr_start-$seq_start, $ltr_end-$ltr_start+1;
	my $seq_len = length $ltr;
	#self-alignment of the candidate LTR region (incl. internal); query and subject are the same seq
	my @Blast=&cached_blast("$cand_idx:self", 1, $ltr, $ltr, $self_spec);
	return if $COLLECT_STAGE==1; #wave 1 only needs to collect the self-alignment job
	my @seq=(split '', $ltr);
	my $motif1="$seq[0]"."$seq[1]";
	my $motif2="$seq[-2]"."$seq[-1]";

	my ($div, $aln_len, $sim, $mismatch, $age, $cor_adj) = (0,0,1,0,0,0);
	my ($q_start, $q_end, $qlen, $s_start, $s_end, $slen, $ls, $le, $rs, $re, $ll, $rl)=(0,0,0,0,0,0,0,0,0,0,0,0);
	my ($nident, $btop, $qseqid, $sseqid) = (0, '', '', '');
	my $adjust="NO";
	$decision="false" if $#Blast==0;

	if ($#Blast>0){
	my $pair=0; #0 indicates no alignment pair seems correct, 1 indicates at least 1 alignment pair seems right
	my $aln_diff = 0;
	for (my $i=1; defined $Blast[$i+1]; $i++){
		$Blast[$i]=~s/^\s+//;
		$decision="false" if $i>8;
		last if $i>8;
		# print "$Blast[$i]\n"; #test
		# Chr1:106472..118130|Chr1:106522..118080 Chr1:106472..118130|Chr1:106522..118080 1       3085    11559   8475    11559   11559   3085    3084    2955CA129

		($qseqid, $sseqid, $s_start, $s_end, $slen, $q_start, $q_end, $qlen, $aln_len, $nident, $btop) = (split /\s+/,  $Blast[$i]); #btop=Blast trace-back operations, contains alignment info
		$cor_adj=$info[0]-1;
		($ls, $le, $rs, $re)=($info[3]-$cor_adj, $info[4]-$cor_adj, $info[6]-$cor_adj, $info[7]-$cor_adj);
		($q_start, $q_end)=($q_end, $q_start) if $q_start>$q_end;
		($s_start, $s_end)=($s_end, $s_start) if $s_start>$s_end;
		($s_start, $s_end, $q_start, $q_end)=($q_start, $q_end, $s_start, $s_end) if $s_start>$q_start;
		if ($s_start>100 or abs($q_end-length($ltr))>100){ #if LTR alignment shift from the start or end for more than 100bp, it's probably FP
			$pair=0;
			next;
			} else {
			$pair=1;
			$aln_diff = abs((($s_start-$s_end)-($q_start-$q_end))); #for paired alignments, calculate the difference of length
			last;
			}
		}
	$decision="false" if $pair==0;

	# count variants
	$btop =~ s/\d+//g; #remove all matches
	my $len_snp = length $btop;
	return unless $len_snp % 2 == 0; #expect string is even length
	
	# count transitions and transversions
	my $n_transition = 0; #A<->G; C<->T
	my $n_transversion = 0; #A<->C; A<->T; G<->C; G<->T
	my $n_indel = 0; #[AGCT] <-> -
	while ($btop =~ s/([ATCG-][ATCG-])//i){
		my $snp = $1;
		$n_transition++ if $snp =~ /(AG)|(GA)|(CT)|(TC)/i;
		$n_transversion++ if $snp =~ /(AC)|(CA)|(AT)|(TA)|(GC)|(CG)|(GT)|(TG)/i;
		$n_indel++ if $snp =~ /\-/;
		}

	# estimate evolutionary distance
	my $tot_len = $n_transition + $n_transversion + $nident; #SNP only, indel not counted
	$tot_len = $seq_len if $tot_len == 0; #EDTA issue 564
	my $raw_d = ($n_transition+$n_transversion) / $tot_len; #percent SNP
	my $JC69_d = 1; #the Jukes-Cantor model K= -3/4*ln(1-4*d/3) adjusts for non-coding sequences, d=$raw_d
	if ($raw_d < 0.66){ #highly diverged sequence could not be adjusted by the JC69 model
		$JC69_d = -3/4*log(1-4*$raw_d/3); #log=ln
		} else {
		$JC69_d = $raw_d;
		}
	my $P = $n_transition / $tot_len; #fraction of transition
	my $Q = $n_transversion / $tot_len; #fraction of transversion
	my $K2P_d = -1/2*log((1-2*$P-$Q)*sqrt(1-2*$Q)); #The Kimura 2-parameter model controls difference b/t transition and transversion rates
        
	#estimate divergence time T = K/2u, where K stands for divergence rate, and u is mutation rate (per bp per ya)
	my $raw_T = sprintf ("%.0f", $raw_d/(2*$miu));
	my $JC69_T = sprintf ("%.0f", $JC69_d/(2*$miu));
	my $K2P_T = sprintf ("%.0f", $K2P_d/(2*$miu));
	$raw_d = sprintf ("%.4f", $raw_d);
	$JC69_d = sprintf ("%.4f", $JC69_d);
	$K2P_d = sprintf ("%.4f", $K2P_d);

	# reassign value to the array
	($div, $age) = ($K2P_d, $K2P_T) if $model eq "K2P";
	($div, $age) = ($JC69_d, $JC69_T) if $model eq "JC69";
	($div, $age) = ($raw_d, $raw_T) if $model eq "PDIST";
	$info[9]=sprintf ("%.4f", 1-$div);
	$info[19]=sprintf ("%.0f", $age);
	#print "$pair\t$K2P_d, $K2P_T, $JC69_d, $JC69_T, $raw_d, $raw_T\n"; #test

	if ($s_end != $le){
		my $i=1;
		for (; $i<100; $i++){
			my $seed="$seq[$s_end-$i-1]"."$seq[$s_end-$i]";
			if ($seed=~/$motif2/i){
				$le=$s_end-$i+1;
				$adjust="3' lLTR";
				last;
				}
			}
		}
	if ($q_start != $rs){
		my $i=1;
		for (; $i<100; $i++){
			my $seed="$seq[$q_start-$i]"."$seq[$q_start-$i+1]";
			if ($seed=~/$motif1/i){
				$rs=$q_start+$i-1;
				$adjust="5' rLTR";
				last;
				}
			}
		}
	$ll=$le-$ls+1;
	$rl=$re-$rs+1;

##update element information in %info
	if (abs(abs($ll-$rl)-$aln_diff)<=$length_diff and $s_end<$q_start){
		$info[3]=$ls+$cor_adj;
		$info[4]=$le+$cor_adj;
		$info[5]=$ll;
		$info[6]=$rs+$cor_adj;
		$info[7]=$re+$cor_adj;
		$info[8]=$rl;
		} else {
		$adjust="NO";
		$ll=$rl=$ls=$le=$rs=$re="NA";
		$decision="false";
		}
	}
##Finish correcting internal boundaries
	
##Start structural analysis
	my ($ltr1_s, $ltr1_e, $ltr2_s, $ltr2_e)=@info[3,4,6,7];
	my ($up1_seq, $do1_seq, $up2_seq, $do2_seq)=('','','','');
	$up1_seq=substr $seq, 0, $ltr1_s-$seq_start+10; #start (i.e. 50bp upstream) + 10bp lLTR
	$do1_seq=substr $seq, $ltr1_e-$seq_start-9, 60; #10bp lLTR + 50bp internal
	$up2_seq=substr $seq, $ltr2_s-$seq_start-50, 60; #50bp internal + 10bp rLTR
	$do2_seq=substr $seq, $ltr2_e-$seq_start-9; #10bp rLTR + end (i.e. 50bp downstream)

##boundary missing rate control
	my $bond_miss=0;
	foreach my $bond ($up1_seq, $do1_seq, $up2_seq, $do2_seq){
		$bond_miss++ while $bond=~/[nN\-]/gi;
		$decision=~s/raw/false/ and last if $bond_miss>=$boundary_N;
		$bond_miss=0;
		}

##boundary alignment
	my ($ac, $bc, $bd, $ad)=(' ',' ',' ',' ');#alignment results between regions
	my ($LTR1_up, $LTR1_do, $LTR2_up, $LTR2_do)=(' ',' ',' ',' ');
#	----||||||||----....----||||||||----
#	 a   5'LTR   b       c   3'LTR   d
#	up60[1]   do60[1]   up60[2]   do60[2]

	$LTR1_up=">$chr:$id\[1]\\n$up1_seq";
	$LTR1_do=">$chr:$id\[1]\\n$do1_seq";
	$LTR2_up=">$chr:$id\[2]\\n$up2_seq";
	$LTR2_do=">$chr:$id\[2]\\n$do2_seq";
	&align_flanking($a_cutoff, $s_cutoff, $w_size, $boundary_ctrl, $LTR1_up, $LTR2_up, $flank_spec, "$cand_idx:ac"); #register the flanking jobs
	&align_flanking($a_cutoff, $s_cutoff, $w_size, $boundary_ctrl, $LTR1_do, $LTR2_do, $flank_spec, "$cand_idx:bd");
	#wave 2 done: hand the boundary-corrected state to process_post (which runs the flanking analysis).
	#Serialize to one string (cheap to store in the shared %STATE) - \x00 separates fields, \x01 marks
	#an undef @info slot (so defined-ness round-trips); neither byte occurs in names/coords/sequences.
	return join("\x00", $decision, $chr, $id, $adjust, $ll, $rl, $s_start, $s_end, $q_start, $q_end,
		$bond_miss, $ltr1_s, $ltr1_e, $ltr2_s, $ltr2_e, $up1_seq, $do1_seq, $up2_seq, $do2_seq,
		map { defined $_ ? $_ : "\x01" } @info);
}

##Wave 3: from the saved boundary-corrected state + the cached flanking blasts, run the flanking
##analysis, decide pass/false/truncated, and emit. (No blast or boundary correction is redone here.)
sub process_post {
	my ($cand_idx, $st) = @_;
	my @f = split /\x00/, $st, -1;  #-1 keeps trailing empty/undef fields
	my ($decision, $chr, $id, $adjust, $ll, $rl, $s_start, $s_end, $q_start, $q_end, $bond_miss, $ltr1_s, $ltr1_e, $ltr2_s, $ltr2_e, $up1_seq, $do1_seq, $up2_seq, $do2_seq) = @f[0..18];
	my @info = map { $_ eq "\x01" ? undef : $_ } @f[19..$#f];
##boundary alignment - rebuild the same flanking inputs process_pre used (the blast itself is cached)
	my $LTR1_up=">$chr:$id\[1]\\n$up1_seq";
	my $LTR1_do=">$chr:$id\[1]\\n$do1_seq";
	my $LTR2_up=">$chr:$id\[2]\\n$up2_seq";
	my $LTR2_do=">$chr:$id\[2]\\n$do2_seq";
	my $ac=&align_flanking($a_cutoff, $s_cutoff, $w_size, $boundary_ctrl, $LTR1_up, $LTR2_up, $flank_spec, "$cand_idx:ac");
	my $bd=&align_flanking($a_cutoff, $s_cutoff, $w_size, $boundary_ctrl, $LTR1_do, $LTR2_do, $flank_spec, "$cand_idx:bd");
	$ac=~s/l3[ATGCN\-?:]+\s+l4[ATCGN\-?:]+\s+(r3[ATGCN\-?:]+\s+r4[ATCGN\-?:]+\s+)HT-align:[0|1]\s+/$1/i;
	$bd=~s/(l3[ATGCN\-?:]+\s+l4[ATCGN\-?:]+\s+)r3[ATGCN\-?:]+\s+r4[ATCGN\-?:]+\s+HT-align:[0|1]\s+/$1/i;

	if ($ac=~/aligned/i or $bd=~/aligned/i){
		$decision="false";
		} else {
		$decision=~s/raw/pass/;
		}

##identify TSD
	my ($TSD_ls, $TSD_le, $TSD_rs, $TSD_re)=(0,0,0,0);
	my ($lTSD, $rTSD)=('','');
	$lTSD=substr $up1_seq, -18, 11; #8bp TSD + 2bp motif + 1bp lLTR
	$rTSD=substr $do2_seq, 7, 11; #1bp rLTR + 2bp motif + 8bp TSD
	my $TSD="NA\t..\t..";
	my $probTSD="NA";
	my $motif="NA";
	my $first_motif="NA"; #first tier: 5bp TSD + TGCA motif
	my $second_motif="NA"; #second tier: 5bp TSD + known non-TGCA motif
	my $third_motif="NA"; #third tier: <5bp TSD + TGCA motif
	my $fourth_motif="NA"; #fourth tier: 5bp TSD + unknown non-TGCA motif
	my $fifth_motif="NA"; #fifth tier: <5bp TSD + non-TGCA motif
	my $sixth_motif="NA"; #sixth tier: >5bp TSD
	foreach my $num (0..6) { #search for longest TSD. Minimum TSD-seed length: 9-6=3.
		foreach (0..(6-$num)) {
			my $len=3+$num;
			next if ($len + $_ > length $lTSD) or (length $lTSD == 0) or (length $rTSD == 0); #avoid substr sequences out of range
			my $seed="NA";
			$seed=substr $lTSD, $_, $len;
			if ($rTSD=~/$seed/i){
				my $temp_motif="NA";
				my $temp_lf_index=$_+$num+3;
				my $temp_rf_index=(index $rTSD, $seed)-2;
				$temp_lf_index=(length $lTSD)-1 if $temp_lf_index>(length $lTSD)-1;
				$temp_rf_index=0 if $temp_rf_index<0;
				$temp_motif=(substr $lTSD, $temp_lf_index, 2).(substr $rTSD, $temp_rf_index, 2);
				$probTSD=$seed if defined $seed; #assign $seed to $probTSD
				
				$first_motif="TGCA_$probTSD" if ($temp_motif=~/TGCA/i and length $probTSD == 5); #first preferred 5bp TSD + TGCA motif
				if (uc $temp_motif ne "TGCA" and length $probTSD == 5){
					foreach my $std_motif (@motif){
						$second_motif="${temp_motif}_$probTSD" if uc $temp_motif eq uc $std_motif; #identify 5bp TSD + known non-TGCA motif
						}
					$fourth_motif="${temp_motif}_$probTSD" if $second_motif eq "NA"; #identify 5bp TSD + unknown non-TGCA motif
					}

				$third_motif="${temp_motif}_$probTSD" if ($temp_motif=~/TGCA/i and (length $probTSD < 5) and (length $probTSD > (length $third_motif) -5) ) ; #get the logest TSD (<5)
				$fifth_motif="${temp_motif}_$probTSD" if (uc $temp_motif ne "TGCA" and (length $probTSD < 5) and (length $probTSD > (length $third_motif) -5) ); #get the logest TSD (<5)
				$sixth_motif="${temp_motif}_$probTSD" if (length $probTSD > 5 and (length $probTSD > (length $third_motif) -5) ); #get the logest TSD (>5)
				}
			}
		}

#Get TSD and motif with preference
	if ($first_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $first_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		} elsif ($second_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $second_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		} elsif ($third_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $third_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		} elsif ($fourth_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $fourth_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		} elsif ($fifth_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $fifth_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		} elsif ($sixth_motif ne "NA"){
		($motif, $probTSD)=($1, $2) if $sixth_motif=~/^([ATCGN]+)_([ATCGN]+)$/;
		}

	my $TSDlen=length $probTSD;

##Correction of TSDs > 5bp
	if ($TSDlen>5 and lc $decision eq "pass") {
		foreach my $std_motif (@motif){
			my ($lm, $rm)=($1, $2) if $std_motif=~/(..)(..)/;
			my $lstart=rindex $lTSD, $lm;
			my $rstart=index $rTSD, $rm;
			if ($lstart>=6 and $rstart<=3){
				my $tsd1=substr $lTSD, $lstart-5, 5;
				my $tsd2=substr $rTSD, $rstart+2, 5;
				if (uc $tsd1 eq uc $tsd2){
					$probTSD=$tsd1;
					$motif=$std_motif;
					$TSDlen=length $probTSD;
					}
				}
			}
		}
	$motif="NA" if (length $motif ne 4);

##Adjust original coordinates
	my $l_adj=index($lTSD, $probTSD)+$TSDlen-8; #8 is the ori start of lLTR
	my $r_adj=index($rTSD, $probTSD)-3; #3 is the ori end or rLTR
	$probTSD="NA" if ((abs($l_adj)+ abs($r_adj) >= 3) and $TSDlen < 4);

	if (uc $probTSD ne "NA"){
		my ($lm, $rm)=("NA", "NA");
		($lm, $rm)=($1, $2) if $motif=~/(..)(..)/;
		$TSD_ls=$ltr1_s + $l_adj + index($lTSD, $probTSD) - rindex($lTSD, $lm);
		$TSD_le=$TSD_ls + $TSDlen - 1;
		my $rindex=index($rTSD, $probTSD) + index($rTSD, $rm) - 1;
		$TSD_rs=$ltr2_e - $r_adj + index($rTSD, $probTSD) + index($rTSD, $rm) - 3;
		$TSD_re=$TSD_rs + $TSDlen - 1;
		$ltr1_s=$ltr1_s+$l_adj;
		$ltr2_s=$ltr2_s+$l_adj;
		$ltr1_e=$ltr1_e+$r_adj;
		$ltr2_e=$ltr2_e+$r_adj;

##adjust the TSD and motif sequence based on element direction
		if (defined $info[12]){
		if ($info[12] eq '-'){
			$motif=~tr/tgcaTGCA/acgtACGT/;
			$motif=reverse $motif;
			$motif="NA" if $motif=~/^TN$/i; #correct 'NA'
			$probTSD=~tr/tgcaTGCA/acgtACGT/;
			$probTSD=reverse $probTSD;
			}
			}
		$TSD="$probTSD\t$TSD_ls..$TSD_le\t$TSD_rs..$TSD_re";

##update coordinates and structural info
#($start, $end, $len, $ls, $le, $ll, $rs, $re, $rl, $sim, $id)
#start end len lLTR_str lLTR_end lLTR_len rLTR_str rLTR_end rLTR_len similarity seqid chr direction TSD lTSD rTSD motif order superfamily age(ya)
#10030396  10042892  12497  10030396  10031396  1001  10041892  10042892  1001  0.985  0  Chr1  NA  CATAC  10030391..10030395  10042893..10042897  TGCA  LTR  Gypsy  27361111
		$info[0]=$ltr1_s;
		$info[1]=$ltr2_e;
		$info[2]=$ltr2_e-$ltr1_s+1;
		$info[3]=$ltr1_s;
		$info[4]=$ltr1_e;
		$info[5]=$ltr1_e-$ltr1_s+1;
		$info[6]=$ltr2_s;
		$info[7]=$ltr2_e;
		$info[8]=$ltr2_e-$ltr2_s+1;
		$info[11]=$chr;
		$info[13]=$probTSD;
		$info[14]="$TSD_ls..$TSD_le";
		$info[15]="$TSD_rs..$TSD_re";
		$info[16]=$motif;
		}

##fill in these variables if not definded
	$info[10]="NA" unless defined $info[10]; #seqid from LTRharvest
	$info[11]="NA" unless defined $info[11]; #chr
	$info[12]="?" unless defined $info[12]; #strand
	$info[13]="NA" unless defined $info[13]; #TSD seq
	$info[14]="NA" unless defined $info[14]; #TSD left coor
	$info[15]="NA" unless defined $info[15]; #TSD right coor
	$info[16]="NA" unless defined $info[16]; #motif
	$info[17]="NA" unless defined $info[17]; #order
	$info[18]="unknown" unless defined $info[18]; #superfamily
	$info[19]="NA" unless defined $info[19]; #age (ya)
	$info[19]="0" if $info[19] eq "-0";

##TSD control, boundary control, MISC control, and reporting
	my $internal=($ltr1_e+1)."..".($ltr2_s-1);
	my $overlap=99; #initial value
	$overlap=abs($ltr1_s-$TSD_le-1)+ abs($TSD_rs-$ltr2_e-1) if (defined $TSD_le and defined $ltr1_s);
	$decision="false" if ($info[12] eq "?" and lc $info[18] eq "unknown" and uc $info[17] eq "NA" and $overlap>0); #?+unknown+NA+runin = false
	$decision="false" if ($motif !~ /TGCA/i and $overlap>0 and abs($info[8])<$minlen); #Nmtf+runin+LTR length<100=false
	$decision="false" if ($motif !~ /TGCT|TACA|TACT|TGGA|TGGT|TATA|TGTA|TGCC|TGCA/i and $info[12] eq "?" and lc $info[18] eq "unknown" and uc $info[17] eq "NA");
	$decision="false" if ($TSD_ctrl==1 and $TSD=~/NA/);
	$decision="false" if ($motif !~ /TGCA/i and length $probTSD ne 5 and $info[12] eq "?" and lc $info[18] eq "unknown" and uc $info[17] eq "NA");
	$decision="false" if uc $info[17] eq "NOTLTR";

	if ($boundary_ctrl and lc $decision eq "pass"){
		unless ($ac=~/right/i and $bd=~/left/i and $motif ne "NA"){
			$decision="truncated";
			}
		unless (uc $motif eq "TGCA" or (length $probTSD == 5 and $motif =~ /^T/i)){ #TGCA + varying length TSD = pass; Txxx + 5bp TSD = pass
			$decision="truncated";
			}
		}

	#last four variables: strand/superfamily/order/age
	my $locus = $info[12] eq '-' ? "$chr:$ltr2_e..$ltr1_s" : "$chr:$ltr1_s..$ltr2_e";
	my $defalse = "$locus\t$decision\tmotif:$motif\tTSD:$TSD\tIN:$internal\t$info[9]\t$info[12]\t$info[18]\t$info[17]\t$info[19]
	Adjust: $adjust\tlLTR: $ll\trLTR: $rl
	Alignment regions: $s_start, $s_end, $q_start, $q_end
	LTR coordinates: $ltr1_s, $ltr1_e, $ltr2_s, $ltr2_e
	TSD-LTR overlap: $overlap
	Boundary missing: $bond_miss\n\n";
	$DEFALSE{$cand_idx} = $defalse; #collected for deterministic ordered printing after join
	$scn{$id} = shared_clone([@info]);
}

##Run one analysis wave: hand every candidate index to a pool of worker threads. Each wave
##re-scans all candidates; the body decides (via $COLLECT_STAGE/$EMIT) what to do this pass.
sub run_wave {
	$idx_queue = Thread::Queue->new();
	$idx_queue->enqueue(0..$#candidates);
	$idx_queue->end();
	my @workers = map { threads->create(\&Identifier) } (1..$threads);
	$_->join() foreach @workers;
}

##Memoized blastn. Returns the cached output lines if available; otherwise, while collecting
##the matching stage, records the job (query seq + subject seq + blast parameters) for the next
##batch. Returns an empty list until the job has been run by flush_blast().
sub cached_blast {
	my ($key, $stage, $query_seq, $subject_seq, $spec) = @_;
	return @{$BLAST_CACHE{$key}} if exists $BLAST_CACHE{$key};
	if ($stage == $COLLECT_STAGE){
		lock(%PENDING);
		$PENDING{$key} = shared_clone([$query_seq, $subject_seq, $spec]) unless exists $PENDING{$key};
		}
	return ();
}

##Run all pending blast jobs, pooled into chunks. Each job is tiny (one self-alignment or one
##flank pair), so a process per job is dominated by blastn start-up. Instead we bin-pack jobs
##into chunks (balanced by sequence length) and run ONE blastn per chunk: every chunk query vs
##every chunk subject, then keep only each job's hit against its OWN subject (qseqid==sseqid,
##both tagged with the job index). This amortizes process start-up over a whole chunk while each
##job's kept hits are identical to its standalone blast. Chunks run in a pool of $threads workers.
sub flush_blast {
	return unless %PENDING;
	my @keys = keys %PENDING;
	my $M = scalar @keys;
	return unless $M;
	# assign each job a numeric index; gather seqs (uniform spec within a flush)
	my (@jkey, @jq, @js, $spec);
	for (my $j=0; $j<$M; $j++){
		my $p = $PENDING{$keys[$j]};
		$jkey[$j]=$keys[$j]; $jq[$j]=$p->[0]; $js[$j]=$p->[1]; $spec=$p->[2];
		}
	# pick a chunk size: large enough to amortize blastn start-up, small enough that the wasted
	# within-chunk cross-hits stay cheap. The cross-hit rate is data-dependent (60bp flank pairs
	# match by chance a lot; self-alignments match however many family members exist), so estimate
	# it from a quick pilot and size chunks to a roughly constant wasted-output budget per chunk.
	my $cs = &estimate_chunk_size(\@jq, \@js, $spec, $M);
	my $nchunks = int(($M + $cs - 1)/$cs);
	$nchunks = $threads if $nchunks < $threads;
	$nchunks = $M       if $nchunks > $M;
	# greedy bin-pack jobs (longest first) into the lightest chunk by total sequence length
	my @order = sort { (length($jq[$b])+length($js[$b])) <=> (length($jq[$a])+length($js[$a])) } (0..$M-1);
	my @chunk_of; my @load = (0) x $nchunks;
	foreach my $j (@order){
		my $min=0; for (my $c=1; $c<$nchunks; $c++){ $min=$c if $load[$c]<$load[$min]; }
		$chunk_of[$j]=$min; $load[$min]+=length($jq[$j])+length($js[$j]);
		}
	# build per-chunk job lists and enqueue (data travels in the shared queue, not cloned per thread)
	my @members; push @{$members[$chunk_of[$_]]}, $_ for (0..$M-1);
	my $cq = Thread::Queue->new();
	foreach my $c (0..$nchunks-1){
		next unless $members[$c] and @{$members[$c]};
		my @data = map { [$_, $jq[$_], $js[$_], $jkey[$_]] } @{$members[$c]};
		$cq->enqueue(shared_clone([$spec, \@data]));
		}
	$cq->end();
	my $nw = $nchunks < $threads ? $nchunks : $threads;
	my @workers = map { threads->create(\&blast_chunk_worker, $cq) } (1..$nw);
	$_->join() foreach @workers;
	%PENDING = ();
}

##Estimate a good chunk size for this batch. One blastn per chunk does an all-vs-all of the chunk,
##but we only keep each job's hit against its OWN subject; the rest (cross-hits) is wasted output
##whose volume is chunk_size^2 * cross_rate per chunk. Bigger chunks amortize start-up but inflate
##that waste. We run a small pilot all-vs-all to measure cross_rate, then pick chunk_size so the
##wasted output per chunk is ~$cross_budget lines (which balances start-up vs waste). Chunk size
##does NOT affect the result (each job's kept hits are identical) - only speed/memory.
sub estimate_chunk_size {
	my ($jq, $js, $spec, $M) = @_;
	return $chunk_max if $M <= 200;                 #small batch: not worth a pilot
	my $is_flank = ($spec =~ /evalue 1000/);
	my $P = $is_flank ? 128 : 64;                   #flank pairs are 60bp (cheap pilot); self is kb
	$P = $M if $P > $M;
	my $qf = "$tmpdir/pilot.q.fa"; my $sf = "$tmpdir/pilot.s.fa";
	open(my $Q, '>', $qf) or return $chunk_max;
	open(my $S, '>', $sf) or return $chunk_max;
	for (my $i=0; $i<$P; $i++){ my $j = int($i*$M/$P); print $Q ">$i\n$jq->[$j]\n"; print $S ">$i\n$js->[$j]\n"; }
	close $Q; close $S;
	my @out = qx(timeout -s KILL $timeout ${blastplus}blastn -subject $sf -query $qf $spec 2> /dev/null);
	unlink $qf, $sf;
	my $self = 0; foreach (@out){ my ($q,$s) = split /\t/, $_, 3; $self++ if defined $s and $q eq $s; }
	my $cross = scalar(@out) - $self;
	my $pairs = $P * ($P - 1);                       #ordered cross pairs in the pilot all-vs-all
	my $rate  = $pairs > 0 ? $cross / $pairs : 0;
	$rate = 1e-6 if $rate < 1e-6;
	my $cross_budget = 5000;                         #target wasted output lines per chunk
	my $cs = int(sqrt($cross_budget / $rate) + 0.5); #chunk_size^2 * rate == budget
	$cs = 16         if $cs < 16;
	$cs = $chunk_max if $cs > $chunk_max;
	return $cs;
}

##Worker: pull a chunk, write its queries and subjects (tagged by job index) to per-worker temp
##files, run one blastn (all queries vs all subjects), then split the output back to each job -
##keeping only that job's self hits (qseqid==sseqid) so the result equals its standalone blast.
sub blast_chunk_worker {
	my ($cq) = @_;
	local $/ = "\n";
	my $tid = threads->tid();
	my $qfile = "$tmpdir/c$tid.q.fa";
	my $sfile = "$tmpdir/c$tid.s.fa";
	while (defined(my $chunk = $cq->dequeue())){
		my ($spec, $data) = @$chunk;
		open(my $qfh, '>', $qfile) or die "ERROR: $qfile: $!\n";
		open(my $sfh, '>', $sfile) or die "ERROR: $sfile: $!\n";
		foreach my $d (@$data){ print $qfh ">$d->[0]\n$d->[1]\n"; print $sfh ">$d->[0]\n$d->[2]\n"; }
		close $qfh; close $sfh;
		my $cmd = "timeout -s KILL $timeout ${blastplus}blastn -subject $sfile -query $qfile $spec 2> /dev/null";
		my @out=();
		for (my $try=0; $try<10; $try++){ #retry guards against transient blast failures
			@out=qx($cmd);
			last if $? == 0;
			}
		# keep each job's hits against its own subject (qseqid==sseqid), in blastn output order
		my %byjob;
		foreach my $line (@out){
			my ($q, $s) = split /\t/, $line, 3;
			next unless defined $s and $q eq $s;
			push @{$byjob{$q}}, $line;
			}
		lock(%BLAST_CACHE);
		foreach my $d (@$data){ $BLAST_CACHE{$d->[3]} = shared_clone($byjob{$d->[0]} // []); }
		}
}

##Inlined from bin/align_flanking.pl: overall aligned length + identity of two flanking
##sequences. Logic is verbatim; only the blastn call is routed through the batched cache, and
##the result is returned instead of printed. $File3 is the subject, $File4 the query (the two
##records are joined by a literal "\n" by the caller, hence the /\\n/ split).
sub align_flanking {
	my ($a_cutoff, $s_cutoff, $w_size, $boundary_ctrl, $File3, $File4, $spec, $key) = @_;
	my ($left_align, $right_align) = ('None', 'None');
	my $boundary_aln = "NA";
	$File3=~s/^\s+//;
	$File4=~s/^\s+//;
	my $seq3=(split /\\n/, $File3)[1];
	my $seq4=(split /\\n/, $File4)[1];
	my $seq_l = (length($seq3) <= length($seq4)) ? length($seq3) : length($seq4);
	chomp ($seq3, $seq4);
	my $left3=substr $seq3, 0, 10;
	my $left4=substr $seq4, 0, 10;
	my $right3=substr $seq3, -10;
	my $right4=substr $seq4, -10;
	my $bond="l3:$left3 l4:$left4 r3:$right3 r4:$right4";
	my @left3=split //, $left3;
	my @left4=split //, $left4;
	my @right3=split //, $right3;
	my @right4=split //, $right4;
	my ($i, $j, $a, $b, $align_cutoff)=(0, 0, 0, 0, 7);
	foreach my $base (@left3){ $a++ if $base eq $left4[$i]; $i++; }
	foreach my $base (@right3){ $b++ if $base eq $right4[$j]; $j++; }
	$left_align="left" if $a>=$align_cutoff;
	$right_align="right" if $b>=$align_cutoff;
	$boundary_aln="$left_align,"."$right_align" if $boundary_ctrl==1;

	## Align the flanking sequence (batched blastn; subject=File3/seq3, query=File4/seq4)
	my @Blast=&cached_blast($key, 2, $seq4, $seq3, $spec);

	my (%q_bank, %s_bank, $sim, $q_start, $q_end, $s_start, $s_end, $q_index, $s_index);
	my ($length, $similarity, $m, $q_sim, $n, $s_sim, $s_e_align);
	$sim=$q_start=$q_end=$s_start=$s_end=$similarity=$q_sim=$s_sim=$s_e_align=0.00;
	$m=$n=$length=1;
	$q_index=$s_index='';

	if (@Blast){
		foreach (@Blast){
			s/^\s+//;
			($sim, $q_start, $q_end, $s_start, $s_end)=(split)[2,6,7,8,9];
			$q_index="$q_index"."$q_start..$q_end;";
			$s_index="$s_index"."$s_start..$s_end;";
			($q_start, $q_end)=($q_end, $q_start) if $q_start>$q_end;
			($s_start, $s_end)=($s_end, $s_start) if $s_start>$s_end;
			$s_e_align=1 if (10-$q_start>=4 && $s_end-(length($seq3)-10)>=4); #head of query aligns with tail of subject
			for (my $i=$q_start; $i<=$q_end; $i++){ $q_bank{"q_$i"}=$sim unless exists $q_bank{"q_$i"}; }
			for (my $j=$s_start; $j<=$s_end; $j++){ $s_bank{"s_$j"}=$sim unless exists $s_bank{"s_$j"}; }
			}
		while ((my $k2, my $value)=each (%q_bank)){ $m++; $q_sim+=$value; }
		while ((my $k2, my $value)=each (%s_bank)){ $n++; $s_sim+=$value; }
		$length=($n+$m)/(2*$seq_l);
		$similarity=($q_sim+$s_sim)/($m+$n);
		$q_sim=sprintf("%.2f", $q_sim/$m);
		$s_sim=sprintf("%.2f", $s_sim/$n);
		} else {
		$length=0;
		}

	my $result = ($length>=$a_cutoff && $similarity>=$s_cutoff) ? "aligned" : "not match";
	return "$result\tBoundary-align:$boundary_aln\t$bond\tHT-align:$s_e_align\tm:$m\tqsim:$q_sim\tn:$n\tssim:$s_sim\tqindex:$q_index\tsindex:$s_index";
}

