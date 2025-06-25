#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX ":sys_wait_h";

# Shujun Ou (shujun.ou.1@gmail.com) and ChatGPT o4-mini
# 06/24/2025
# Description: Overcome memory bottleneck on Repeatmasker jobs with large library and large target sequence by:
# 	Split sequence into small chunks and execute Repeatmasker on them.
# 	Skip chunks if the .masked file exists for that chunk.
# 	Return the merged .masked file at the end. 

# defaults
my ($input, $lib, $jobs, $cpus, $chunk_size, $out_prefix, $repeatmasker, $blastplus);
$jobs       = 2; # assume memory consumption of 40GB per job
$cpus       = 4;
$chunk_size = 1_000_000;
$out_prefix = 'combined';
$blastplus  = 'makeblastdb';  # assume in PATH

# parse args
GetOptions(
    'input=s'        => \$input,
    'lib=s'          => \$lib,
    'jobs=i'         => \$jobs,
    'cpus=i'         => \$cpus,
    'chunk_size=i'   => \$chunk_size,
    'out_prefix=s'   => \$out_prefix,
    'repeatmasker=s' => \$repeatmasker,
    'blastplus=s'    => \$blastplus,
) or usage();
usage() unless $input && $lib;

# locate RepeatMasker
# allow override via -repeatmasker / RepeatMasker path
if (!$repeatmasker) {
    chomp($repeatmasker = `which RepeatMasker 2>/dev/null`);
}
$repeatmasker =~ s/RepeatMasker\n?$//;      # strip trailing program name
$repeatmasker .= '/' if $repeatmasker && $repeatmasker !~ /\/$/;
die "Error: RepeatMasker not found at ${repeatmasker}RepeatMasker\n"
    unless -x "${repeatmasker}RepeatMasker";

# build BLAST DB for your library once
(my $lib_base = $lib) =~ s/\.(fa|fasta|lib)$//i;
unless (-e "$lib_base.nsq") {
    warn "Building BLAST DB for library '$lib' …\n";
    system($blastplus, "-in", $lib, "-dbtype", "nucl") == 0
      or die "Error: makeblastdb failed on $lib\n";
}

# prepare workdir
my $workdir = "RM_split";
mkdir $workdir or warn "Can't mkdir $workdir: $!";

# SPLIT FASTA into whole-record chunks
open my $IN, '<', $input or die "Can't read $input: $!";
my $chunk_idx  = 0;
my $base_count = 0;
my $OUT;
while (my $line = <$IN>) {
    if ($line =~ /^>/) {
        # header: new chunk if size reached
        if (!defined $OUT or $base_count >= $chunk_size) {
            close $OUT if defined $OUT;
            $chunk_idx++;
            open $OUT, '>', sprintf("%s/chunk_%03d.fa", $workdir, $chunk_idx)
              or die "Can't write chunk: $!";
            $base_count = 0;
        }
        print $OUT $line;
    } else {
        my $seq = $line; chomp $seq;
        $seq =~ s/\s+//g;
        $base_count += length($seq);
        print $OUT "$seq\n";
    }
}
close $IN;
close $OUT if defined $OUT;

# collect chunks
opendir my $DH, $workdir or die $!;
my @chunks = sort grep { /^chunk_\d+\.fa$/ } readdir $DH;
closedir $DH;

# fork management
my %kids;
sub reap { while((my $p = waitpid(-1, WNOHANG))>0){ delete $kids{$p} } }

for my $chunk (@chunks) {
    # skip if this chunk has already been masked
    my $masked = "$workdir/" . $chunk;
    $masked =~ s/\.fa$/.fa.masked/;
    if ( -e $masked ) {
        warn "Skipping $chunk (found $masked)\n";
        next;
    }

    # throttle to $jobs
    while (keys %kids >= $jobs) {
        reap();
        sleep 1;
    }
    my $pid = fork();
    die "fork failed: $!" unless defined $pid;
    if ($pid == 0) {
        # child
        exec(
            "${repeatmasker}RepeatMasker", "-e","ncbi", "-q", "-no_is", "-norna", "-div","40", "-cutoff","225",
            "-pa",$cpus,
            "-lib",$lib,
            "$workdir/$chunk"
        ) or die "exec failed: $!";
    } else {
        $kids{$pid} = 1;
    }
}

# wait for all
1 while keys %kids and reap(), sleep 1;

# MERGE masked outputs
my $merged = "$out_prefix.masked";
open my $MER, '>', $merged or die "Can't write $merged: $!";
for my $chunk (@chunks) {
    (my $m = $chunk) =~ s/\.fa$/.fa.masked/;
    open my $INM, '<', "$workdir/$m" or die "Can't open $m: $!";
    print $MER $_ while <$INM>;
    close $INM;
}
close $MER;

`rm $lib.* 2>/dev/null`;
print "All done! Combined masked file: $merged\n";
exit 0;

sub usage {
    die <<"EOF";
Usage: $0 --input FASTA --lib LIB [--jobs N] [--cpus M] [--chunk_size B] \\
          [--out_prefix P] [--repeatmasker PATH] [--blastplus PATH]
Options:
  --input       target FASTA to mask
  --lib         RepeatMasker library
  --jobs   N    max parallel jobs (default $jobs)
  --cpus   M    CPUs per job          (default $cpus)
  --chunk_size B bases per chunk      (default $chunk_size)
  --out_prefix P prefix for merged    (default '$out_prefix')
  --repeatmasker PATH  path to RM dir (default from `which RepeatMasker`)
  --blastplus PATH     path to makeblastdb
EOF
}

