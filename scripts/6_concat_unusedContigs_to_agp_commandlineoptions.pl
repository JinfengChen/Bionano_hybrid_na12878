#!/usr/bin/perl -w

use strict;
use warnings;

#my $agpDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/agp";
#my $agpFile = "selectedNewV1Seq.vs.v2Hybrid.fastaId.twoPrimaryContigScaffolds.agp";
my $agpDir = $ARGV[0];
my $agpFile = $ARGV[1];

#my $fastaDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/new_pac_bio_seq";
#my $fastaFile = "consensus_diploid_for_andy_forncbi_filtered2b.fa";
my $fastaDir = $ARGV[2];
my $fastaFile = $ARGV[3];

my $outDir = "$agpDir";
my $outFile = $agpFile;	$outFile =~ s/agp$/notUsed.agp/;

my $fastaInfoRef = getFastaInfo($fastaDir, $fastaFile);
my @agp = ();	my $agpRef = \@agp;
($fastaInfoRef, $agpRef) = getAgpInfo($agpDir, $agpFile, $agpRef, $fastaInfoRef);
$agpRef = addUnusedSeq($agpRef, $fastaInfoRef);
printAgpInfo($outDir, $outFile, $agpRef);

sub printAgpInfo	{
	my ($dir, $file, $agpRef) = @_;
	open(OUT, ">$dir/$file") or die "printAgpInfo: writing to $dir/$file: $!\n";	
my $count = 0;	print STDERR "printAgpInfo: writing to $file\n";
	for (my $i = 0; $i < scalar(@$agpRef); $i += 1)	{
		my $theLine = $agpRef->[$i];
		print OUT "$theLine\n";
$count += 1 if ($theLine !~ /^#/);
	} # for i
print STDERR "\twritten $count records\n";
	close OUT;
} # printAgpInfo

sub addUnusedSeq	{
	my ($agpRef, $fastaInfoRef) = @_;
my $count = 0;
	foreach my $seqName (sort keys %$fastaInfoRef)	{
		my $fRef = $fastaInfoRef->{$seqName};
		if ($fRef->{usedInAgp} == 0)	{
			# append to agpRef, those sequence that were not used in the scaffolding process
			# must distinguish object name from component name
			my $objName = $seqName."_obj";
			my $newLine = join("\t", $objName, 1, $fRef->{size}, 1, "W", $seqName, 1, $fRef->{size}, "+");
			push(@$agpRef, $newLine);	
$count += 1;
		} # if fRef
	} # foreach seqName
print STDERR "addUnusedSeq: added $count unused sequences\n";
	return $agpRef;
} # addUnusedSeq

sub getAgpInfo	{
	my ($dir, $file, $agpRef, $fastaInfoRef) = @_;
	open(IN, "$dir/$file") or die "getAgpInfo: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getAgpInfo: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header lines
			push(@$agpRef, $line);
			next;
		} # if line
		my @content = split(/\t/, $line);
		my ($compType, $compId) = ($content[4], $content[5]);
		if ($compType =~ /W/i)	{
			# this is a component line, update the fastaInfoRef 
			die "getAgpInfo: the component sequence=$compId does NOT have an entry in the original fasta file\n" if (! exists $fastaInfoRef->{$compId});
			$fastaInfoRef->{$compId}{usedInAgp} = 1;
$count += 1;
		} # if compType
		push(@$agpRef, $line);
	} # while line
print STDERR "\tread in $count records\n";
	close IN;	
	return ($fastaInfoRef, $agpRef);
} # getAgpInfo

sub getFastaInfo	{
	my ($dir, $file) = @_;
	my %fastaInfo = ();	# fastaInfo{seqName} = {size, usedInAgp}
	open(IN, "$dir/$file") or die "Cannot open getFastaInfo: $!\n";
my $count = 0;	print STDERR "getFastaInfo: reading in $file\n";
	my ($curName, $curSeq) = ("", "");
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^>(.+)/)	{
			# if a header line
			if ($curName ne "" && $curSeq ne "")	{
				$fastaInfo{$curName} = {size => length($curSeq), usedInAgp => 0};
				($curName, $curSeq) = ("", "");
$count += 1;
			} # if curName
			$curName = $1;
		} else	{
			$curSeq .= $line if ($curName ne "");
		} # if line
	} # while line
	# last sequence
	if ($curName ne "" && $curSeq ne "")	{
		$fastaInfo{$curName} = {size => length($curSeq), usedInAgp => 0};
$count += 1;
	} # curName
print STDERR "\tgetFastaInfo: read in $count records\n";
	close IN;
	return \%fastaInfo;
} # getFastaInfo
