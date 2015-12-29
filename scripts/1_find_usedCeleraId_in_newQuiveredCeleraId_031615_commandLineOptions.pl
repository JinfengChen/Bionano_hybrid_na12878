#!/usr/bin/perl -w

use strict;
use warnings;

sleep(5);

# this is to enable the translation of the Celera contigs constituted to hybrid scaffolding to new Quivered Celera contigs

# thie list of v1 ngs contigs that actually participated in scaffolding
#my $v1UsedNgsDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/write_up/NatureMethods/script_release_oldWays/v1/output_results/mergeNGS_BN";

#my $v1UsedNgsDir = "/sc/orga/scratch/pendlm02/testBioNano/oldWays/script_release_oldWays/v1/output_results/mergeNGS_BN";
my $v1UsedNgsDir = $ARGV[0]; # the directory of ngs contigs used in the original pipeline

my $v1UsedNgsFile = "ngsContigs.txt"; # file in the v1 direotyr indicating ngs contigs

# the translation of old Celera contigs to cmap (assuming the format is fa2cmap format)
#my $oldKeyDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/write_up/NatureMethods/script_release_oldWays/v1/output_results"
#my $oldKeyDir = "/sc/orga/scratch/pendlm02/testBioNano/oldWays/script_release_oldWays/v1/output_results";
my $oldKeyDir = $ARGV[1]; # the directory in v1 of ngs key file used in the original pipeline
my $oldKeyFile = $ARGV[2]; # the key file in the oldKeyDir ... "mar3_NA12878.scf_BspQI_0Kb_0labels_key.txt"

###This step requires running of fa2cmap on the original fasta
# the translation of quivered Celera contigs (OR illumina corrected OR alternatively updated) to cmap
#my $newHomeDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015";
my $newHomeDir = $ARGV[3]; # this is where new data is AND where files will be written to "/sc/orga/scratch/pendlm02/testBioNano/../rerun_coassembly_050814/v2_output_T10/output_results/for_ncbi_031015";
my $newKeySubDir = $ARGV[4]; # subdirectory in output of new key files "$newHomeDir/new_pac_bio_seq";
my $newKeyDir = "$newHomeDir/$newKeySubDir";
my $newKeyFile = $ARGV[5];  # $"consensus_diploid_for_andy_forncbi_filtered2b_BspQI_0Kb_0labels_key.txt";

my $outDir = "$newHomeDir/align_usedQuiveredCelera_v2Hybrid"; 

mkdir $outDir if (! -e $outDir);
my $usedFile = "used.translation.fasta.cmap.txt";
my $newUsedFile = "notUsed.translation.fasta.cmap.txt";

my $usedCeleraOldIdsRef = getUsedCeleraOldIds($v1UsedNgsDir, $v1UsedNgsFile);
my ($usedOldKeysRef, $notUsedOldKeysRef) = getOldKeys($oldKeyDir, $oldKeyFile, $usedCeleraOldIdsRef);
($usedOldKeysRef, $notUsedOldKeysRef) = getNewKeys($newKeyDir, $newKeyFile, $usedOldKeysRef, $notUsedOldKeysRef);
printOutput($outDir, $usedFile, $newUsedFile, $usedOldKeysRef, $notUsedOldKeysRef);

sub printOutput	{
	my ($dir, $usedFile, $notUsedFile, $usedOldKeysRef, $notUsedOldKeysRef) = @_;
	open(OUT, ">$dir/$usedFile") or die "printOutput: cannot write to $dir/$usedFile: $!\n";
	print OUT "oldCmapId\tfastaId\tnewCmapId\toldSeqLength\tnewSeqLength\n";
my $count = 0;
	open(OUT2, ">$dir/$notUsedFile") or die "printOutput: cannot write to $dir/$notUsedFile: $!\n";
	print OUT2 "oldCmapId\tfastaId\tnewCmapId\n";
my $count2 = 0;
	foreach my $fastaId (keys %$usedOldKeysRef)	{
		my $aRef = $usedOldKeysRef->{$fastaId};
		print OUT "$aRef->{oldCmapId}\t$fastaId\t$aRef->{newCmapId}\t$aRef->{oldSize}\t$aRef->{newSize}\n";
$count += 1;
	} # foreach fastaId

	foreach my $fastaId (keys %$notUsedOldKeysRef)	{
		my $aRef = $notUsedOldKeysRef->{$fastaId};
		print OUT2 "$aRef->{oldCmapId}\t$fastaId\t$aRef->{newCmapId}\t$aRef->{oldSize}\t$aRef->{newSize}\n";
$count2 += 1;
	} # foreach fastaId
print STDERR "\tprinted $count used records and $count2 not used records\n";
	close OUT;	close OUT2;
} # printOutput

sub getNewKeys	{
	my ($dir, $file, $usedOldKeysRef, $notUsedOldKeysRef) = @_;
	open(IN, "$dir/$file") or die "getNewKeys: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getNewKeys: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		next if ($line !~ /^\d+/);
		my @content = split(/\t/, $line);
		my ($cmapId, $fastaId, $size) = @content;
		if (exists $usedOldKeysRef->{$fastaId})	{
			# used ones
			$usedOldKeysRef->{$fastaId}{newCmapId} = $cmapId;
			$usedOldKeysRef->{$fastaId}{newSize} = $size;
$count += 1;
		} elsif (exists $notUsedOldKeysRef->{$fastaId})	{
			# not used ones
			$notUsedOldKeysRef->{$fastaId}{newCmapId} = $cmapId;
			$notUsedOldKeysRef->{$fastaId}{newSize} = $size;
		} else	{
			# add to the not used ones
			$notUsedOldKeysRef->{$fastaId} = {oldCmapId => -1, newCmapId => $cmapId, oldSize => -1, newSize => $size};
		} # if exists
	} # while line
print STDERR "\tread in $count used records (should be the same as above)\n";	
	close IN;
	return ($usedOldKeysRef, $notUsedOldKeysRef);
} # getNewKeys

# takes in Heng's translation file format (not in this version)
#sub getOldKeys	{
#	my ($dir, $file, $usedCeleraOldIdsRef) = @_;
#	my %usedOldKeys = ();	my %notUsedOldKeys = ();
#	open(IN, "$dir/$file") or die "Cannot open $dir/$file: $!\n";
#my $count = 0;	print STDERR "getOldKeys: reading in $file\n";
#	my $skip = <IN>;
#	while (my $line = <IN>)	{	
#		chomp $line;
#		$line =~ s/\r//g;
#		$line =~ s/"//g;
#		$line =~ s/\s+/\t/g;
#		my @content = split(/\t/, $line);
#		my ($size, $cmapId, $fastaId) = @content[1..3];
#		if (exists $usedCeleraOldIdsRef->{$cmapId})	{
#			$usedOldKeys{$fastaId} = {oldCmapId => $cmapId, newCmapId => -1, oldSize => $size, newSize => -1};
#$count += 1;
#		} else	{
#			$notUsedOldKeys{$fastaId} = {oldCmapId => $cmapId, newCmapId => -1, oldSize => $size, newSize => -1};
#		} # if exists
#	} # while line
#print STDERR "\tread in $count used entries (should be the same as above count)\n";
#	close IN;
#	return (\%usedOldKeys, \%notUsedOldKeys);
#} # getOldKeys

# takes in fa2cmap translation key file format (in this version)
sub getOldKeys	{
	my ($dir, $file, $usedCeleraOldIdsRef) = @_;
	my %usedOldKeys = ();	my %notUsedOldKeys = ();
	open(IN, "$dir/$file") or die "Cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getOldKeys: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		next if ($line !~ /^\d+/);
		$line =~ s/\r//g;
		$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		my ($cmapId, $fastaId, $size) = @content;
		if (exists $usedCeleraOldIdsRef->{$cmapId})	{
			$usedOldKeys{$fastaId} = {oldCmapId => $cmapId, newCmapId => -1, oldSize => $size, newSize => -1};
$count += 1;
		} else	{
			$notUsedOldKeys{$fastaId} = {oldCmapId => $cmapId, newCmapId => -1, oldSize => $size, newSize => -1};
		} # if exists
	} # while line
print STDERR "\tread in $count used entries\n";
	close IN;
	return (\%usedOldKeys, \%notUsedOldKeys);
} # getOldKeys

sub getUsedCeleraOldIds	{
	my ($dir, $file) = @_;
	my %usedCeleraOldIds = ();
	open(IN, "$dir/$file") or die "getUsedCeleraOldIds: cannot open $dir/$file: $!\n";
my $count = 0;	print STDERR "getUsedCeleraOldIds: reading in $file\n";
	my $skip = <IN>;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		$usedCeleraOldIds{$line} = 1;
$count += 1;
	} # while line
print STDERR "\tread in $count records\n";
	close IN;
	return \%usedCeleraOldIds;
} # getUsedCeleraOldIds

