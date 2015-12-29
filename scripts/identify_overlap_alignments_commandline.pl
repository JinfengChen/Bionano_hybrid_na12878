#!/usr/bin/perl -w

use strict;
use warnings;

# identify cases when two sequence contigs overlaps (a - alignment touch, b - alignment touch;sequence overlap, c - alignment overlap, d - alignment do not touch, but sequence overlap)

my $homeDir = $ARGV[0]; 
my $inDir = "$homeDir";
my $xmapDir = $ARGV[1]; 
my $xmapFile = $ARGV[2];
#my $xmapFile = $ARGV[2]; # xmap of alignments you are trying to finding overlaps for # 
#my $xmapFile = "selectedNewV1Seq.vs.v2Hybrid.xmap";	# no chimera, no not used NGS
my $ngsCmapFile = "selectedNewV1Seq.vs.v2Hybrid_q.cmap";
#my $ngsCmapFile = $ARGV[3]; # cmp of query matching above xmap

my $fastaDir = "$homeDir/fasta/";	mkdir $fastaDir if (! -e $fastaDir);
my $outDir = "$homeDir/fasta/ngs_overlap_status";	mkdir $outDir if (! -e $outDir);
#my $outDir = $ARGV[2];	mkdir $outDir if (! -e $outDir);
my $outFile = "adjacent_ngs_overlap_status.txt";


my $xmapRef = getXmap($xmapDir, $xmapFile);
$xmapRef = sortByCoord($xmapRef);

my $cmapRef = getCmap($inDir, $ngsCmapFile);

my $overlapRef = getOverlap($xmapRef, $cmapRef);

printOverlap($outDir, $outFile, $overlapRef);

sub printOverlap	{
	my ($dir, $file, $overlapRef) = @_;
	open(OUT, ">$dir/$file") or die "printOverlap: cannot write to $dir/$file: $!\n";
my $count = 0;	print STDERR "printOverlap: writing to $file\n";
	print OUT "anchorId\tqId1\tstart1\tend1\tqStart1\tqEnd1\torientation1\t";
	print OUT "qId2\tstart2\tend2\tqStart2\tqEnd2\torientation2\t";
	print OUT "overlapType\n";
	foreach my $id (sort numeric keys %$overlapRef)	{
		for (my $i = 0; $i < scalar(@{$overlapRef->{$id}}); $i += 1)	{
			my $aRef = $overlapRef->{$id}[$i];
			print OUT join("\t", $id, $aRef->{curQId}, $aRef->{curStart}, $aRef->{curEnd}, $aRef->{curQStart}, $aRef->{curQEnd}, $aRef->{curOrientation})."\t";
			print OUT join("\t", $aRef->{nextQId}, $aRef->{nextStart}, $aRef->{nextEnd}, $aRef->{nextQStart}, $aRef->{nextQEnd}, $aRef->{nextOrientation})."\t";
			print OUT join("\t", $aRef->{type})."\n";
$count += 1;
		} # for i
	} # foreach id
print STDERR "\t written $count records\n";
	close OUT;
} # printOverlap

sub numeric	{$a	<=>	$b}

sub getOverlap	{
	my ($xmapRef, $cmapRef) = @_;
	my %overlap = ();

	foreach my $id (keys %$xmapRef)	{
		for (my $i = 0; $i < scalar(@{$xmapRef->{$id}}) - 1; $i += 1)	{
			my $curContigRef = $xmapRef->{$id}[$i];
			my $nextContigRef = $xmapRef->{$id}[$i + 1];
			
			if ($curContigRef->{end} < $nextContigRef->{start})	{
				# if two alignments do not touch
				# check if there are unaligned sequence that overlap
				die "getOverlap 1: cannot find cmap size information for qId = $curContigRef->{qId}\n" if (! exists $cmapRef->{$curContigRef->{qId}});
				die "getOverlap 2: cannot find cmap size information for qId = $nextContigRef->{qId}\n" if (! exists $cmapRef->{$nextContigRef->{qId}});
				my $curUnalignSize = ($curContigRef->{orientation} =~ /^\+$/) ? ($cmapRef->{$curContigRef->{qId}} - $curContigRef->{qEnd}) : ($curContigRef->{qEnd} - 1);
				my $nextUnalignSize = ($nextContigRef->{orientation} =~ /^\+$/) ? ($nextContigRef->{qStart} - 1) : ($cmapRef->{$nextContigRef->{qId}} - $nextContigRef->{qStart});				
				
				# now see if the unaligned portion of the query sequences would potentially overlap
				if ($curUnalignSize + $nextUnalignSize > $nextContigRef->{start} - $curContigRef->{end})	{
					# unaligned portions are larger than the distance between the aligned portions, thus overlap
					push(@{$overlap{$id}}, {type => "alignment-untouch;sequence-overlap", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
				} elsif ($curUnalignSize + $nextUnalignSize == $nextContigRef->{start} - $curContigRef->{end})	{
					# unaligned portions size is equal to the distance between the aligned portions, thus touch
					push(@{$overlap{$id}}, {type => "alignment-untouch;sequence-touch", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
				} else	{
print STDERR "curUnalignSize=$curUnalignSize;nextUnalignSize=$nextUnalignSize;distance=".($nextContigRef->{start} - $curContigRef->{end}).";\n" if ($curContigRef->{qId} == 10596 && $nextContigRef->{qId} == 10594);
					# no overlapping sequence
					push(@{$overlap{$id}}, {type => "alignment-untouch;sequence-untouch", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
						
				} # if unaligned portion can potentially overlap

			} elsif ($curContigRef->{end} == $nextContigRef->{start})	{
				# if two alignments touch
				die "getOverlap 3: cannot find cmap size information for qId = $curContigRef->{qId}\n" if (! exists $cmapRef->{$curContigRef->{qId}});
				die "getOverlap 4: cannot find cmap size information for qId = $nextContigRef->{qId}\n" if (! exists $cmapRef->{$nextContigRef->{qId}});
				my $curUnalignSize = ($curContigRef->{orientation} =~ /^\+$/) ? ($cmapRef->{$curContigRef->{qId}} - $curContigRef->{qEnd}) : ($curContigRef->{qEnd} - 1);
				my $nextUnalignSize = ($nextContigRef->{orientation} =~ /^\+$/) ? ($nextContigRef->{qStart} - 1) : ($cmapRef->{$nextContigRef->{qId}} - $nextContigRef->{qStart});				
				
				# now see if the unaligned portion of the query sequences would potentially overlap
				if ($curUnalignSize > 0 || $nextUnalignSize > 0)	{
					# additional sequence sticks out beyond the last label
					push(@{$overlap{$id}}, {type => "alignment-touch;sequence-overlap", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
				} else	{
					# no additional sequence still sticks out
					push(@{$overlap{$id}}, {type => "alignment-touch;sequence-overlap", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
				} # if curUnalignSize
			} else	{
				# two alignments overlap
					push(@{$overlap{$id}}, {type => "alignment-overlap;sequence-overlap", 
						curQId => $curContigRef->{qId}, curQStart => $curContigRef->{qStart}, curQEnd => $curContigRef->{qEnd}, curStart => $curContigRef->{start}, curEnd => $curContigRef->{end}, curOrientation => $curContigRef->{orientation}, 
						nextQId => $nextContigRef->{qId}, nextQStart => $nextContigRef->{qStart}, nextQEnd => $nextContigRef->{qEnd}, nextStart => $nextContigRef->{start}, nextEnd => $nextContigRef->{end}, nextOrientation => $nextContigRef->{orientation}});
				
			} # if 
		} # for i
	} # foreach id
	return \%overlap;
} # getOverlap
	
sub getCmap	{
	my ($dir, $file) = @_;
	my %cmap = ();

	open(IN, "$dir/$file") or die "getCmap: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getCmap: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r//g;
		next if ($line =~ /^#/);
		my @content = split(/\t/, $line);
		my ($qId, $qSize) = @content[0..1];
		$cmap{$qId} = $qSize if (! exists $cmap{$qId});
$count += 1;
	} # while line
print STDERR "\tread in $count records\n";
	close IN;
	return \%cmap;
} # 

sub sortByCoord	{
	my ($dataRef) = @_;
	foreach my $id (keys %$dataRef)	{
		@{$dataRef->{$id}}	=	sort	{
			$a->{start}	<=> $b->{start}
		} @{$dataRef->{$id}}
	} # foreach id
	return $dataRef;
} # sortByCoord

sub getXmap	{
	my ($dir, $file) = @_;
	my %xmap = ();
	open(IN, "$dir/$file") or die "getXmap: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getXmap: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r//g;	
		next if ($line =~ /^#/);
		my @content = split(/\t/, $line);
		my ($qId, $chrom, $qStart, $qEnd, $start, $end, $orientation) = ($content[1], $content[2], int($content[3]), int($content[4]), int($content[5]), int($content[6]), $content[7]);	
		push(@{$xmap{$chrom}}, {start => $start, end => $end, qId => $qId, qStart => $qStart, qEnd => $qEnd, orientation => $orientation});
$count += 1;
	} # while line
print STDERR "\tread in $count records\n";
	close IN;
	return \%xmap;
} # getXmap
