#!/usr/bin/perl -w

use strict;
use warnings;

#my $agpDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015/align_usedQuiveredCelera_v2Hybrid/alignref_noRepeatMask/agp";
#my $agpFile = "selectedNewV1Seq.vs.v2Hybrid.fastaId.twoPrimaryContigScaffolds.notUsed.agp";
my $agpDir = $ARGV[0];
my $agpFile = $ARGV[1];

my $outDir = "$agpDir";
my $outAgpFile = $agpFile;	$outAgpFile =~ s/agp$/trimHeadTailGap.agp/;
my $outTranslationFile = $agpFile;	$outTranslationFile =~ s/agp$/trimHeadTailGap.coordTable/;

my ($headerRef, $agpInfoRef) = getAgpInfo($agpDir, $agpFile);
$agpInfoRef = findHeadTailGap($agpInfoRef);
printNewAgp($outDir, $outAgpFile, $headerRef, $agpInfoRef);
printTranslation($outDir, $outTranslationFile, $agpInfoRef);

sub printTranslation	{
	my ($dir, $file, $agpInfoRef) = @_;
	open(OUT, ">$dir/$file") or die "printTranslation: cannot write to $dir/$file: $!\n";
my $count = 0;	print STDERR "printTranslation: writing to $file\n";
	print OUT "ObjId\tMap_coord\tAgp_coord\n";
	foreach my $objName (sort keys %$agpInfoRef)	{
		#$agpInfo{$objId} = {headGapFlag => 0, tailGapFlag => 0, oriAgpHeadObjStart => 1, newAgpHeadObjStart => 1, offSet => 0}; # initialize
		my $hRef = $agpInfoRef->{$objName};
		print OUT join("\t", $objName, $hRef->{oriAgpHeadObjStart}, $hRef->{newAgpHeadObjStart})."\n";	
$count += 1;
	} # foreach objName
print STDERR "\twritten $count records\n";
	close OUT;
} # printTranslation

sub printNewAgp	{
	my ($dir, $file, $headerRef, $agpInfoRef) = @_;
	open(OUT, ">$dir/$file") or die "printNewAgp: cannot write to $dir/$file: $!\n";	
my $count = 0;	print STDERR "printNewAgp: writing to $file\n";
	foreach my $header (@$headerRef)	{
		print OUT "$header\n";
	} # foreach header
	foreach my $objName (sort keys %$agpInfoRef)	{
		my $hRef = $agpInfoRef->{$objName};
		my $startIndex = ($hRef->{headGapFlag} == 1) ? (1) : (0);
		my $endIndex = ($hRef->{tailGapFlag} == 1) ? ($#{$hRef->{content}} - 1) : ($#{$hRef->{content}});
		for (my $i = $startIndex; $i <= $endIndex; $i += 1)	{
			my $iRef = $hRef->{content}[$i];
			my ($newObjStart, $newObjEnd) = ($iRef->{objStart} - $hRef->{offSet}, $iRef->{objEnd} - $hRef->{offSet});
			my $newPartNum = ($hRef->{headGapFlag} == 1) ? ($iRef->{partNum} - 1) : ($iRef->{partNum});
			print OUT join("\t", $objName, $newObjStart, $newObjEnd, $newPartNum, $iRef->{compType}, $iRef->{compId}, $iRef->{compStart}, $iRef->{compEnd}, $iRef->{orientation})."\n";
$count += 1;
		} # for i
	} # foreach objName
print STDERR "\twritten $count records\n";
	close OUT;
} # printNewAgp

sub findHeadTailGap	{
	my ($agpInfoRef) = @_;
	foreach my $objName (keys %$agpInfoRef)	{
		if ($agpInfoRef->{$objName}{content}[0]{compType} =~ /N/i)	{
			# if the first segment is a gap
			# there has to be more than one segment for this to be a valid hybrid scaffold
			die "findHeadTailGap: $objName object has a head gap, but there is only 1 segment in this object\n" if (scalar(@{$agpInfoRef->{$objName}{content}}) == 1);
			# the second segment must be sequence
			die "findHeadTailGap: $objName object has a head gap, but the next segment is not a sequence\n" if ($agpInfoRef->{$objName}{content}[1]{compType} !~ /W/i); 
			$agpInfoRef->{$objName}{headGapFlag} = 1;
			$agpInfoRef->{$objName}{oriAgpHeadObjStart} = $agpInfoRef->{$objName}{content}[1]{objStart};	
			$agpInfoRef->{$objName}{offSet} = $agpInfoRef->{$objName}{content}[0]{objEnd};			# always subtract this number
		} # if agpInfoRef
		if ($agpInfoRef->{$objName}{content}[$#{$agpInfoRef->{$objName}{content}}]{compType} =~ /N/i)	{
			# if the last segment is a gap
			# there has to be more than one segment for this to be a valid hybrid scaffold
			die "findHeadTailGap: $objName object has a tail gap, but there is only 1 segment in this object\n" if (scalar(@{$agpInfoRef->{$objName}{content}}) == 1);
			$agpInfoRef->{$objName}{tailGapFlag} = 1;
		} # if agpInfoRef

		if ($agpInfoRef->{$objName}{headGapFlag} == 1 && $agpInfoRef->{$objName}{tailGapFlag} == 1)	{
			# if there is a head gap and a tail gap, then there must be more two segments
			die "findHeadTailGap: $objName has head and tail gap, but there are no more than 2 segment in this object\n" if (scalar(@{$agpInfoRef->{$objName}{content}}) <= 2);
		} # if agpInfoRef
	} # foreach objName
	return $agpInfoRef;
} # findHeadTailGap

sub getAgpInfo	{
	my ($dir, $file) = @_;
	my %agpInfo = ();	my @header = ();
	open(IN, "$dir/$file") or die "getAgpInfo: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getAgpInfo: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			push(@header, $line);
			next;
		} # if line

		my ($objId, $objStart, $objEnd, $partNum, $compType, $compId, $compStart, $compEnd, $orientation) = split(/\t/, $line);
		$agpInfo{$objId} = {headGapFlag => 0, tailGapFlag => 0, oriAgpHeadObjStart => 1, newAgpHeadObjStart => 1, offSet => 0, content => []} if (! exists $agpInfo{$objId});	# initialize
		push(@{$agpInfo{$objId}{content}}, {objStart => $objStart, objEnd => $objEnd, partNum => $partNum, compType => $compType, compId => $compId, compStart => $compStart, compEnd => $compEnd, orientation => $orientation});
$count += 1;		
	} # while line
print STDERR "There are ".(scalar(@{$agpInfo{"Super-Scaffold_2"}{content}}))." segments in Super-Scaffold_2\n";
print STDERR "\tread in $count records\n";	
	close IN;
	return (\@header, \%agpInfo);
} # getAgpInfo
