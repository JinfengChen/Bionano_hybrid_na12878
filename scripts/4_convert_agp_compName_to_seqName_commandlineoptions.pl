#!/usr/bin/perl -w

use strict;
use warnings;
# this is for submission to ncbi, so that we convert the cmap id of the contigs back to fasta id so that there is a relation between the sequence submitted and the v2 scaffold's components

#my $homeDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015";
#my $homeDir = $ARGV[0];
my $agpDir = $ARGV[0];
my $agpFile = $agpDir . $ARGV[1];
#my $agpDir = "$homeDir/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/agp";
#my $agpFile = "selectedNewV1Seq.vs.v2Hybrid.agp";

#my $keyDir = "$homeDir/new_pac_bio_seq";
#my $keyFile = "consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels_key.txt";

my $keyDir = $ARGV[2];
my $keyFile = $ARGV[3];

#my $agpFastaIdDir = "$agpDir";
#my $agpFastaIdFile = "selectedNewV1Seq.vs.v2Hybrid.fastaId.agp";
my $agpFastaIdDir = "$agpDir";
my $agpFastaIdFile = $ARGV[4];
my $translationRef = getTranslation($keyDir, $keyFile);

# make the translation
open(IN, "$agpFile") or die "Cannot open $agpFile: $!\n";
open(OUT, ">","$agpFastaIdDir/$agpFastaIdFile") or die "Cannot write to $agpFastaIdDir/$agpFastaIdFile: $!\n";
my $count = 0;	my $outCount = 0;
print STDERR "reading in $agpFile and writing to $agpFastaIdFile\n";
while (my $line = <IN>)	{
	chomp $line;
$count += 1;
	if ($line =~ /^#/)	{
		print OUT "$line\n";
$outCount += 1;
	} else	{
		my @content = split(/\t/, $line);
		if ($content[6] =~ /^scaffold$/i)	{
			# inter contig lines
			print OUT "$line\n";
$outCount += 1;
		} else	{
			# contig lines
			my $compId = $content[5];
			die "Cannot find translation for componentId=$compId\n" if (! exists $translationRef->{$compId});
			print OUT join("\t", @content[0..4])."\t$translationRef->{$compId}\t".join("\t", @content[6..$#content])."\n";
$outCount += 1;
		} # if 
	} # if line
} # while line
print STDERR "\tread in $count records; written out $outCount records\n";
close IN;	close OUT;



sub getTranslation	{
	my ($dir, $file) = @_;
	my %translation = ();
	open(IN, "$dir/$file") or die "getTranslation: cannot open file $dir/$file: $!\n";
	my $skip = <IN>;
my $count = 0;	print STDERR "getTranslation: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		my @content = split(/\t/, $line);
		my ($cmapId, $fastaId) = @content[0..1];
		$translation{$cmapId} = $fastaId;
$count += 1;
	} # while line
	close IN;
print STDERR "\tread in $count records\n";
	return \%translation;
} # getTranslation
