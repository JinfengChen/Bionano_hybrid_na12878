#!/usr/bin/perl -w

use strict;
use warnings;

#my $ngsDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid";	# use v1 (filtered) NGS in this case
my $ngsDir = $ARGV[0];	# use v1 (filtered) NGS in this case
my $ngsFile = "used.newSeq.cmap";

#my $superScafDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/write_up/NatureMethods/script_release_oldWays/v2/output_results/mergeNGS_BN";
my $superScafDir = $ARGV[1];
my $superScafFile = "step2.hybrid.cmap";

my $outDir = "$ngsDir/alignref_noRepeatMask";	mkdir $outDir if (! -e $outDir);
my $outFilePrefix = "selectedNewV1Seq.vs.v2Hybrid";

#my $cmd = "/home/users/apang/bin/current/RefAligner ";
my $cmd = "$ARGV[2] ";
$cmd .= "-ref $superScafDir/$superScafFile ";
$cmd .= "-i $ngsDir/$ngsFile ";
$cmd .= "-o $outDir/$outFilePrefix ";
$cmd .= "-f -T 1e-10 -A 5 -biaswt 0 -FP 0.5 -FN 0.05 -sf 0.2 -sd 0.10 -M 1 -BestRef 1 -outlier 1e-4 -endoutlier 1e-3 -deltaX 9 -deltaY 9 -xmapchim 14 -maxmem 150 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 -hash -hashdelta 50 -mres 0 -rres 0.002 -res 0 -extend 0 ";

print STDERR "command = $cmd;\n";
my $systemOutput = system($cmd);
print STDERR "systemOutput=$systemOutput;";
