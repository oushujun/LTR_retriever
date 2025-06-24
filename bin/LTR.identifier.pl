#!/usr/bin/env perl
use warnings;
use strict;
use threads;
use Thread::Queue;
use threads::shared;
use File::Basename;

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
		$id=~s/.*\|.*:([0-9]+\.\.[0-9]+)$/$1/;
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

#Store sequence information and push candidates into multithreading queue
my $queue=Thread::Queue->new();
$/ = "\n>";
while (<FA>){
	chomp;
	s/>//g;
	my ($name, $seq) = (split /\n/, $_, 2);
	next unless defined $seq;
	next if $seq eq '';
	$seq=~s/\s+//g;
	$seq=uc $seq;
	my $id = $1 if $name =~ /.*:(\S+)/;
	next unless defined $scn{$id};
	$queue->enqueue([$name, $seq]);
	}
close FA;
$/ = "\n";
$queue -> end();

##open scn.adj file and print out the header
open SCN, ">$List.adj" or die "ERROR: $!";
print SCN "#LTR boundary fine-grain adjustment and annotation have been performed by LTR_retriever (Shujun Ou, oushujun\@msu.edu)\n$head";

##Run the work queue with multiple threads
foreach (1..$threads){
	threads -> create(\&Identifier);
	}
foreach (threads -> list()){
	$_ -> join();
	}

##print out entries that could not pass initial screening criteria to scn.adj
foreach my $key (sort{$a cmp $b}(keys %scn)){
	foreach (0..$#{$scn{$key}}){
		next unless defined $scn{$key};
		print SCN "$scn{$key}[$_]  ";
		}
	print SCN "\n";
	}
close SCN;


##subrotine for LTR structural analyses
sub Identifier() {
	while (defined($_ = $queue->dequeue())){
##Structural analysis, main program
	my ($name, $seq)=(@{$_}[0], @{$_}[1]);
	my $decision="raw"; #conclusion of whether the element is a LTR
	my ($chr, $seq_start, $seq_end, $ltr_start, $ltr_end);
	($chr, $seq_start, $seq_end, $ltr_start, $ltr_end)=($1, $2, $3, $5, $6) if $name=~/^(\S+):([0-9]+)..([0-9]+)\|(\S+):([0-9]+)..([0-9]+)/;  #eg: Chr4:10009589..10017157|Chr4:10009609..10017137 or 10.dna.chromosome.ch:100016935..100026312|10.dna.chromosome.ch:100016935..100026312
	my $id="$ltr_start..$ltr_end";
	next if $id eq '';
	my @info = @{$scn{$id}};

##Coarse boundary correction - after correction, coordinates may still have 1-2 bp shifted from the real case
	my $ltr=substr $seq, $ltr_start-$seq_start, $ltr_end-$ltr_start+1;
	my $seq_len = length $ltr;
	my @seq=(split '', $ltr);
	my $motif1="$seq[0]"."$seq[1]";
	my $motif2="$seq[-2]"."$seq[-1]";
	my $candidate_seq=">$name\\n$ltr"; #just candidate LTR including internal region, no extended sequences
	my $exec="timeout -s KILL $timeout ${blastplus}blastn -subject <(echo -e \"$candidate_seq\") -query <(echo -e \"$candidate_seq\") -dust no -outfmt \"6 qseqid sseqid sstart send slen qstart qend qlen length nident btop\" -parse_deflines";
	my @Blast=();
	for (my $try=0; $try<10; $try++){ #it's possible that sequence wrote in memory is rewritten by other programs and caused blast error, this step will try 10 times to guarantee the blast is run correctly
		@Blast=qx(bash -c '$exec' 2> /dev/null) if defined $ltr;
		last if $? == 0;
		}

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
	next unless $len_snp % 2 == 0; #expect string is even length
	
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
	$ac=`perl $script_path/align_flanking.pl $a_cutoff $s_cutoff $w_size $boundary_ctrl \"$LTR1_up\" \"$LTR2_up\" $blastplus`;
	$bd=`perl $script_path/align_flanking.pl $a_cutoff $s_cutoff $w_size $boundary_ctrl \"$LTR1_do\" \"$LTR2_do\" $blastplus`;
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
	print $defalse;
	$scn{$id} = shared_clone([@info]);
	}
}

