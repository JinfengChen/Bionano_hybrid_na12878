#!/usr/bin/perl -w

use strict;
use warnings;

my $v2SeqFile = "input_data/Citrus_v1_quiver_round2.fasta";
my $v1SeqFile = "input_data/asm.scf.fasta";
my $bngFile = "input_data/asm.scf_BspQI_0Kb_0labels.cmap";
mkdir "v1" if (! -e "v1");
mkdir "v2" if (! -e "v2");

### v1
my $cmd = "/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks hybridScaffold_v1.pl ";
$cmd .= "-n $v1SeqFile ";
$cmd .= "-b $bngFile ";
$cmd .= "-c xml/hybridScaffold_config_v1.xml ";
$cmd .= "-o v1/output_results ";
$cmd .= "-f ";
system($cmd);
print "done v1\n";

### merge results
$cmd = "./RefAligner -merge ";
$cmd .= "-i v1/output_results/mergeNGS_BN/step2.hybrid.cmap ";
$cmd .= "-i v1/output_results/mergeNGS_BN/step1.BN.naive.cmap ";
$cmd .= "-o v2/v1.hybrid.BN.naive ";
$cmd .= "-f ";
system($cmd);
print "done merge\n";

### v2
### trust the BN/hybrid
$cmd = "/rhome/cjinfeng/BigData/software/perl5.10.1/perl-5.10.1/perl -Mforks hybridScaffold_v2.pl ";
$cmd .= "-b v2/v1.hybrid.BN.naive.cmap ";
$cmd .= "-n $v2SeqFile ";
$cmd .= "-c xml/hybridScaffold_config_v2.xml ";
$cmd .= "-o v2/output_results ";
$cmd .= "-B -f ";
system($cmd);
print "done v2\n";
