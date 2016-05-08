#!/Perl/bin/perl 
use Bio::SeqIO;

#################################################################################
#										#
#	Script name:	AMY-tree_v2.0.pl					#
#	Author: 	Anneleen Van Geystelen					#
#	Version:	v1.0							#
#	Created:	02/04/2011						#
#	Last Modified:	17/01/2014						#
#										#
#	Description:	This Perl script will define the haplogroup of a sample #
#			based on the called SNPs of a whole genome sequencing 	#
#			experiment. 					 	#
#										#
#										#
#	INPUT:									#
#	    variable:								#
# 	    fixed: 	path of input file					#
#			path of output directory				#
#			path of tree file					#
#			path of conversion file					#
#			path of hg18 or hg19 Y-chromosome fasta file		#
#			path of status hg18 or hg19 Y-chromosome file		#
#			quality control filePartialOption			#
#			version (i.e. hg18 or hg19)				#
#			file with defined regions if option PARTIAL is chosen	#
#										#
#	OUTPUT: 	output file						#
#										#
# 	PERL MODULES:								#
#										#
#################################################################################

#################################################################################
#				VARIABLES					#
#################################################################################

if ($#ARGV < 7) {
	print "Please, give valid arguments\n"; die;
}


# Set path of called SNP file.
my $inFile = $ARGV[0];
if (!(-f $inFile)) { die "$inFile does not exist #0\n"; }

# Set path of directory for output files. 
my $outDir = $ARGV[1];
if (!(-d $outDir)) { die "$outDir does not exist #1\n"; }

# Set path of the tree file. 
my $treeFile = $ARGV[2];
if (!(-f $treeFile)) { die "$treeFile does not exist #2\n"; }

# Set path of file to the conversion of the SNPs to their position, mutant and ancestral. 
my $convFile = $ARGV[3];
if (!(-f $convFile)) { die "$convFile does not exist #3\n"; }

# Set path of the Y-chromosome fasta file hg18/19.
my $hgFile = $ARGV[4];
if (!(-f $hgFile)) { die "$hgFile does not exist #4\n"; }

# Set path of the Y-chromosome fasta file hg18/19.
my $stathgFile = $ARGV[5];
if (!(-f $stathgFile)) { die "$stathgFile does not exist #5\n"; }

# Set path of quality control file. 
my $qcFile = $ARGV[6];
if (!(-f $qcFile)) { die "$qcFile does not exist #6\n"; }

# Set version variable. 
my $version = $ARGV[7];


# Set version variable. 
my $filePartialOption = $ARGV[8];
if ($filePartialOption ne '') {
    if (!(-f $filePartialOption)) { die "$filePartialOption does not exist #8\n"; }
}


# Set cutoff value for high quality samples. 
$cutoffHighQ = 0.85;

# Set cutoff value for low quality samples. 
$cutoffLowQ = 0.05;

#################################################################################
#				END OF VARIABLES				#
#################################################################################

#################################################################################
#				SUBROUTINES					#
#################################################################################

# Function to trim leading and trailing underscores.
sub trimUnderscore {
	my $string = shift;
	while (substr($string, 0, 1) eq '_') {
		$string = substr($string, 1);
	}
	while (substr($string, -1) eq '_') {
		$string = substr($string, 0,-1);
	}
# 	$string =~ s/^_+//;
# 	$string =~ s/_+$//;
	return $string;
}

# Function to trim leading and trailing underscores.
sub trimStar {
	my $string = shift;
	$string =~ s/^\*+//;
	$string =~ s/\*+$//;
	return $string;
}

# Function to trim leading and trailing spaces.
sub trimSpace {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# Function to check if element is in array.
sub isInArray {
	my ($element, @array) = @_;
	my %array;
	foreach $elem (@array) {
		$array{$elem} = 1;
	}
	my $ret = $array{$element};
}


# Function to get status of node. 
sub getMetric {
	my @results = @_; # e.g. 1 -1 0 -1 1 1 1 
	my $ones = 0;
	my $total = 0;
	foreach $result (@results) {
		if ($result ne '') {
			$total++;
			if ($result == 1) {
				$ones++;
			} 
		}
	}
	my $metric ;
	if ($total != 0) {
		$metric = $ones/$total;
	} 

	return $metric;
}

# Function to get the difference between 2 arrays. 
sub getDiffArrays {
	my ($a1, $a2) = @_;
	my @a1 = split('_', $a1);
	my @a2 = split('_', $a2);
	my %diff1; my %diff2;
	foreach $elem2 (@a2) {
		$diff2{$elem2} = 1;
	}
	foreach $elem1 (@a1) {
		if (!($diff2{$elem1} == 1)) {
			$diff1{$elem1} = 1;
		}
	}
	# %diff1 contains elements from '@a1' that are not in '@a2'
	@k = keys %diff1;
	return @k;
}

# Function to get the difference between 2 arrays. 
sub getDiffArraysSNPs {
	my (@conv, @tree) = @_;
	@conv = uniqueArray(@conv);
	@tree = uniqueArray(@tree);
	my %tree;
	foreach $elem (@tree) {
		$elem = trimSpace($elem);
		$tree{$elem} = 1;
	}
	my @res;
	foreach $elem (@conv) {
		$elem = trimSpace($elem);
		if ($tree{$elem} eq '') {
			push(@res, $elem);
		}
	}
	# @res contains elements from '@conv' that are not in '@tree'
	return @res;
}

# Function to make an array unique.
sub uniqueArray {
    my %seen = ();
    my @r = ();
    foreach my $a (@_) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
    }
    return @r;
}

# Function to arrange layout of output file. 
sub layout {
	my ($string) = @_;
	$string =~ s/^ //;
	$string =~ s/ $//;
	my $newString;
	if ($string =~ m/ or / && substr($string, -1) eq '*') {
		my @string = split(/ or /, substr($string, 0, -1));
		foreach $elem (@string) {
			$newString = $newString.$elem.'* or ';
		}
		$newString = substr($newString, 0, -4);
	} else {
		$newString = $string
	}
	return $newString;
}

# Function to get the alternative name. 
sub altName {
	my ($string) = @_;
	$string =~ s/^ //;
	$string =~ s/ $//;
	my $stringOri = $string;
	my $newString = $altnames{$string};
	while($newString eq '') {
		my $stringTemp = reverse $string;
		$stringTemp =~ s/\*//;
		$string = reverse $stringTemp;
		$newString = $altnames{$string};
	}
	my $ret;
		$ret = $stringOri.' ['.$newString.']';
	return $ret;
}

# Function to run SHELLY.
sub SHELLY {
	($analysisFile, $statSamFile, $treeFile, $statRefFile, $qualFile, $statusOutFile, $manualHaplo) = @_;
	if (!(-f $analysisFile)) { die "$analysisFile does not exist\n"; }
	if (!(-f $statSamFile)) { die "$statSamFile does not exist\n"; }
	if (!(-f $treeFile)) { die "$treeFile does not exist\n"; }
	if (!(-f $statRefFile)) { die "$statRefFile does not exist\n"; }
	###################### Create hash of status file of sample. ####################
	my %name2statSam;
	my %name2origin; 
	open (STATSAM, "< $statSamFile") or die "could not open $statSamFile\n";
	while(<STATSAM>) {
		my $line = $_; 
		chomp($line); 
		my ($name, $status, $origin) = split(/\t/, $line);
		$name2statSam{$name} = $status;
		$name2origin{$name} = $origin;
	}
	close STATSAM;
	#################### Create hash of status file of reference. ###################
	my %name2statRef;
	open (STATREF, "< $statRefFile") or die "could not open $statRefFile\n";
	while(<STATREF>) {
		my $line = $_; 
		chomp($line); 
		my ($name, $status) = split(/\t/, $line);
		$name2statRef{$name} = $status;
	}
	close STATREF;
	########################## Create array of analysis file. ########################
	my @haplogroups; 
	open (ANALYSIS, "< $analysisFile") or die "could not open $analysisFile\n";
	while(<ANALYSIS>) {
		my $line = $_; 
		chomp($line); 
		if (substr($line, 0, 1) eq '>') {
			my $haplogroup = substr($line, 2);
			if ($haplogroup ne 'Root [Root]' && $haplogroup ne 'Root* [Root]' ) {
				push(@haplogroups, $haplogroup); 
			}
		}
	}
	close ANALYSIS; 
	############################ Create hash of tree file. ##########################
	my %node2parent;
	my %node2snp;
	my %parent2nodes;
	my %SNPsInTree;
	open (TREE, "< $treeFile") or die " Could not open tree file: $treeFile\n";
	while (<TREE>) {
		my $line = $_; 
		$line = substr($line, 0,-1);
		my @splittedLine = split(/\t/, $line); 
		my $node = $splittedLine[0];
		my $parent = $splittedLine[2];
		my @SNPs = @splittedLine[3..$#splittedLine];
		my $SNPs;
		foreach $elem (@SNPs) {
			if ($elem ne '' && $elem ne '-' && $elem ne ' ') {
				$elem =~ s/\r//;
				$SNPs = $SNPs.$elem."_";
				$SNPsInTree{$elem} = 1;
			}
		}
		$SNPs = trimUnderscore($SNPs);
		$node2parent{$node} = $parent;
		$node2snp{$node} = $SNPs;
		my $children = $parent2nodes{$parent};
		$children = $children.$node."_";
		$parent2nodes{$parent} = $children;
	} 
	close TREE;
	foreach $parent (sort keys %parent2nodes) {
		my $nodes = $parent2nodes{$parent};
		$nodes = trimUnderscore($nodes);
		$parent2nodes{$parent} = $nodes;
	}
	############################## Calculate quality. ##############################
	my $haplogroup;
	if ($#haplogroups == 0) {
		$haplogroup = $haplogroups[0];
	} 
	if ($manualHaplo ne '') {
		$haplogroup = $manualHaplo;
	}
	open (STATOUT, "> $statusOutFile");
	open (OUT, "> $qualFile");
	if ($haplogroup ne '' && $haplogroup ne 'Root* [Root]' ) {
		# get SNPs on its path
		$haplogroup  =~ s/\[.*\]//;
		$haplogroup = substr($haplogroup, 0, -1);
		my $name = $haplogroup;
		my %SNPsInPath;
		my %SNPsToIgnore;
		if ( $name =~ m/\*/ ) {
			while($name =~ m/\*/) {
				$name =~ s/\*//;
			}
			if (substr($name, -1) eq ' ') {
				$name = substr($name, 0, -1);
			}
			my $parent = $node2parent{$name};
			my @children = split(/\_/, $parent2nodes{$name});
			# To ignore
			my @childrenToIgnore;
			while (@children) {
				my $child = $children[0];
				if ($child ne $name) {
					push(@childrenToIgnore, $child);
					my @grandchildren = split(/\_/, $parent2nodes{$child}); 
					foreach $grandchild (@grandchildren) {
						push(@children, $grandchild);
					}
				}
				shift(@children);
			}
			@children = split(/\_/, $parent2nodes{$parent});
			# To ignore
			while (@children) {
				my $child = $children[0];
				if ($child ne $name) {
					push(@childrenToIgnore, $child);
					my @grandchildren = split(/\_/, $parent2nodes{$child}); 
					foreach $grandchild (@grandchildren) {
						push(@children, $grandchild);
					}
				}
				shift(@children);
			}
			foreach $childToIgnore (@childrenToIgnore) {
				my @SNPs = split(/\_/, $node2snp{$childToIgnore}); 
				foreach $SNP (@SNPs) {
					$SNPsToIgnore{$SNP} = 1;
				}
			}
			# In Path
			my @SNPs = split(/\_/, $node2snp{$name});
			foreach $SNP (@SNPs) {
				$SNPsInPath{$SNP} = 1;
			}
			while($parent ne 'Root') {
				my @SNPs = split(/\_/, $node2snp{$parent});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}
				$parent = $node2parent{$parent};
			}
		} else {
			# In Path
			while(substr($name, -1) eq ' ') {
				$name = substr($name, 0, -1);
			}
			my $parent = $node2parent{$name};
			my @SNPs = split(/\_/, $node2snp{$name});
			foreach $SNP (@SNPs) {
				$SNPsInPath{$SNP} = 1;
			}
			while($parent ne 'Root') {
				my @SNPs = split(/\_/, $node2snp{$parent});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}
				$parent = $node2parent{$parent};
			}
		}
		# Create hash with expected called SNPs. 
		my %CallExp;
		my %NotCallExp;
		foreach $SNP (sort keys %SNPsInTree) {
			if (!($SNPsToIgnore{$SNP})) {
				my $statRef = $name2statRef{$SNP};
				if ($SNPsInPath{$SNP}) {
					if ($statRef != 1) {
						$CallExp{$SNP} = 1;
					} else {
						$NotCallExp{$SNP} = 1;
					}
				} else {
					if ($statRef == 1) {
						$CallExp{$SNP} = 1;
					} else {
						$NotCallExp{$SNP} = 1;
					}
				}
			}
		}
		# Get qualities. 
		my $CallInExp = 0; 
		my $CallInNotExp = 0;
		my $Correct = 0;
		my $NotCorrect = 0;
		foreach $SNP (sort keys %SNPsInTree) {
			if (!($SNPsToIgnore{$SNP})) {
				my $status = $name2statSam{$SNP};
				if ($SNPsInPath{$SNP}) {
					if ($status == 1) {
						$Correct++;
						print STATOUT $SNP."\tTP\n";
					} else {
						$NotCorrect++;
						print STATOUT $SNP."\tFN\n";
					} 
				} else {
					if ($status == 0 || $status == -1) {
						$Correct++;
						print STATOUT $SNP."\tTN\n";
					} else {
						$NotCorrect++;
						print STATOUT $SNP."\tFP\n";
					} 
				}
				my $origin = $name2origin{$SNP};
				if ($CallExp{$SNP} && $origin eq 'called') {
					$CallInExp++;
				} 

				if ($NotCallExp{$SNP} && $origin eq 'called') {
					$CallInNotExp++;
				} 
			}
		}
		my $score1 = $Correct/($Correct + $NotCorrect)*100;
		print OUT "ACC: $score1\n";
		my @CallExp;
		foreach $key (sort keys %CallExp) {
			push(@CallExp, $key);
		}	
		my $CallExp = $#CallExp + 1;
				
		my @NotCallExp;
		foreach $key (keys %NotCallExp) {
			push(@NotCallExp, $key);
		}
		my $NotCallExp = $#NotCallExp + 1;
		my $TP = $CallInExp;
		my $FP = $CallInNotExp;
		my $TN = $NotCallExp - $CallInNotExp;
		my $FN = $CallExp - $CallInExp;
		my $score2 = $TP/($TP+$FN)*100;
		print OUT "SEN: $score2\n";
		my $score3 = $TN/($FP+$TN)*100;
		print OUT "SPE: $score3\n";
		my $score5 = $TP/($TP+$FP)*100;
		print OUT "PRE: $score5\n";
		my $score6 = $TP/($TP+$FN)*100;
		print OUT "REC: $score6\n";
		my $score7 = ($score5 + $score6)/2;
		print OUT "Fsc: $score7\n";
		my $score4 = (($TP*$TN)-($FP*$FN))/sqrt(($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))*100;
		print OUT "MCC: $score4\n";
	} else {
		if (length(@haplogroups) > 1) {
			$qualFile = substr($qualFile, 0, -4);
			$statusOutFile = substr($statusOutFile, 0, -4);
			foreach $haplogroup (@haplogroups) {
				my ($name) = split(/ /,$haplogroup);
				$name =~ s/\*//;
				$qualHgFile = $qualFile. "_". $name.".txt";
				$statusOutHgFile = $statusOutFile . "_". $name.".txt";
				system("perl SHELLY.pl $analysisFile $statSamFile $treeFile $statRefFile $qualHgFile $statusOutHgFile $haplogroup");
			}
		} else {
			print OUT "When there is more than one haplogroup, SHELLY cannot determine the quality\n";
		}
	}
	close OUT;
	close STATOUT;
}

# Function to run SHELLY_v2.0.
sub SHELLY2 {
	my ($analysisFile, $statSamFile, $treeFile, $statRefFile, $regions, $convFile, $version, $qualFile, $statusOutFile, $manualHaplo) = @_;
	if (!(-f $analysisFile)) { die "$analysisFile does not exist\n"; }
	if (!(-f $statSamFile)) { die "$statSamFile does not exist\n"; }
	if (!(-f $treeFile)) { die "$treeFile does not exist\n"; }
	if (!(-f $statRefFile)) { die "$statRefFile does not exist\n"; }
	if (!(-f $regions)) { die "$regions does not exist\n"; }
	if (!(-f $convFile)) { die "$convFile does not exist\n"; }

	my %name2anc;
	my $flag = 0; 
	open (CONV, "< $convFile") or die " Could not open conversion file: $convFile\n";
	while (<CONV>) {
		my $line = $_; 
		chomp($line);
		if ($flag == 1) {
			my ($name, $pos36, $pos37, $mutation, $type, $ignore) = split(/\t/, $line); 
			if (trimSpace($ignore) eq 'no') {
				$mutation = trimSpace($mutation);
				my ($anc, $mut) = split(/->/, $mutation);
				$name = trimSpace($name);
				if ($version eq 'hg18') {
					push(@SNPsInConv, $name);
					$pos36 = trimSpace($pos36);
					$name2anc{$name} = $pos36."_".$anc;
				} elsif ($version eq 'hg19') {
					push(@SNPsInConv, $name);
					$pos37 = trimSpace($pos37);
					$name2anc{$name} = $pos37."_".$anc;
				}
			}
		} else {
			$flag = 1;
		}
	} 
	close CONV;
	###################### Create hash of status file of sample. ####################
	my %name2statSam;
	my %name2origin; 
	open (STATSAM, "< $statSamFile") or die "could not open $statSamFile\n";
	while(<STATSAM>) {
		my $line = $_; 
		chomp($line); 
		my ($name, $status, $origin) = split(/\t/, $line);
		$name2statSam{$name} = $status;
		$name2origin{$name} = $origin;
	}
	close STATSAM;
	#################### Create hash of status file of reference. ###################
	my %name2statRef;
	open (STATREF, "< $statRefFile") or die "could not open $statRefFile\n";
	while(<STATREF>) {
		my $line = $_; 
		chomp($line); 
		my ($name, $status) = split(/\t/, $line);
		$name2statRef{$name} = $status;
	}
	close STATREF;
	########################## Create array of analysis file. ########################
	my @haplogroups; 
	open (ANALYSIS, "< $analysisFile") or die "could not open $analysisFile\n";
	while(<ANALYSIS>) {
		my $line = $_; 
		chomp($line); 
		if (substr($line, 0, 1) eq '>') {
			my $haplogroup = substr($line, 2);
			if ($haplogroup ne 'Root [Root]') {
				push(@haplogroups, $haplogroup); 
			}
		}
	}
	close ANALYSIS; 
	############################ Create hash of tree file. ##########################
	my %node2parent;
	my %node2snp;
	my %parent2nodes;
	my %SNPsInTree;
	open (TREE, "< $treeFile") or die " Could not open tree file: $treeFile\n";
	while (<TREE>) {
		my $line = $_; 
		$line = substr($line, 0,-1);
		my @splittedLine = split(/\t/, $line); 
		my $node = $splittedLine[0];
		my $parent = $splittedLine[2];
		my @SNPs = @splittedLine[3..$#splittedLine];
		my $SNPs;
		foreach $elem (@SNPs) {
			if ($elem ne '' && $elem ne '-' && $elem ne ' ') {
				$elem =~ s/\r//;
				$SNPs = $SNPs.$elem."_";
				$SNPsInTree{$elem} = 1;
			}
		}
		$SNPs = trimUnderscore($SNPs);
		$node2parent{$node} = $parent;

		$node2snp{$node} = $SNPs;

		my $children = $parent2nodes{$parent};
		$children = $children.$node."_";
		$parent2nodes{$parent} = $children;
	} 
	close TREE;
	foreach $parent (sort keys %parent2nodes) {
		my $nodes = $parent2nodes{$parent};
		$nodes = trimUnderscore($nodes);
		$parent2nodes{$parent} = $nodes;
	}
	##################### Create hash of regions in regions file. ###################
	my %partial;
	open (IN, "< $regions") ;
	while(<IN>) {
		my $line = $_;
		chomp($line);
		my (undef, $start,$end) = split(/\t/, $line);
		$partial{$start} = $end;
	}
	close IN;
	my %inRegion;
	foreach $SNP (sort keys %SNPsInTree) {
		my ($position, undef) = split(/_/,$name2anc{$SNP});
		my $inRegion = 0;

		foreach $start (sort keys %partial) {
			my $end = $partial{$start};
			if ($position >= $start && $position <= $end) {
				$inRegion = 1;
			}
		}
		if ($inRegion == 1) {
			$inRegion{$SNP} = 1;
		} 
	}
	############################## Calculate quality. ##############################
	my $haplogroup;
	if ($#haplogroups == 0) {
		$haplogroup = $haplogroups[0];
	} 
	if ($manualHaplo ne '') {
		$haplogroup = $manualHaplo;
	}
	open (STATOUT, "> $statusOutFile");
	open (OUT, "> $qualFile");
	if ($haplogroup ne '') {
		# get SNPs on its path
		$haplogroup  =~ s/\[.*\]//;
		$haplogroup = substr($haplogroup, 0, -1);
		my $name = $haplogroup;
		my %SNPsInPath;
		my %SNPsToIgnore;
		if ( $name =~ m/\*/ ) {
			while($name =~ m/\*/) {
				$name =~ s/\*//;
			}
			if (substr($name, -1) eq ' ') {
				$name = substr($name, 0, -1);
			}
			my $parent = $node2parent{$name};
			my @children = split(/\_/, $parent2nodes{$name});
			# To ignore
			my @childrenToIgnore;
			while (@children) {
				my $child = $children[0];
				if ($child ne $name) {
					push(@childrenToIgnore, $child);
					my @grandchildren = split(/\_/, $parent2nodes{$child}); 
					foreach $grandchild (@grandchildren) {
						push(@children, $grandchild);
					}
				}
				shift(@children);
			}
			@children = split(/\_/, $parent2nodes{$parent});
			# To ignore
			while (@children) {
				my $child = $children[0];
				if ($child ne $name) {
					push(@childrenToIgnore, $child);
					my @grandchildren = split(/\_/, $parent2nodes{$child}); 
					foreach $grandchild (@grandchildren) {
						push(@children, $grandchild);
					}
				}
				shift(@children);
			}
			foreach $childToIgnore (@childrenToIgnore) {
				my @SNPs = split(/\_/, $node2snp{$childToIgnore}); 
				foreach $SNP (@SNPs) {
					$SNPsToIgnore{$SNP} = 1;
				}
			}
			# In Path
			my @SNPs = split(/\_/, $node2snp{$name});
			foreach $SNP (@SNPs) {
				$SNPsInPath{$SNP} = 1;
			}
			while($parent ne 'Root') {
				my @SNPs = split(/\_/, $node2snp{$parent});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}
				$parent = $node2parent{$parent};
			}
		} else {
			# In Path
			while(substr($name, -1) eq ' ') {
				$name = substr($name, 0, -1);
			}
			my $parent = $node2parent{$name};
			my @SNPs = split(/\_/, $node2snp{$name});
			foreach $SNP (@SNPs) {
				$SNPsInPath{$SNP} = 1;
			}
			while($parent ne 'Root') {
				my @SNPs = split(/\_/, $node2snp{$parent});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}
				$parent = $node2parent{$parent};
			}
		}
		# Create hash with expected called SNPs. 
		my %CallExp;
		my %NotCallExp;
		foreach $SNP (sort keys %inRegion) {
			if (!($SNPsToIgnore{$SNP})) {
				my $statRef = $name2statRef{$SNP};
				if ($SNPsInPath{$SNP}) {
					if ($statRef != 1) {
						$CallExp{$SNP} = 1;
					} else {
						$NotCallExp{$SNP} = 1;
					}
				} else {
					if ($statRef == 1) {
						$CallExp{$SNP} = 1;
					} else {
						$NotCallExp{$SNP} = 1;
					}
				}
			}
		}
		# Get qualities. 
		my $CallInExp = 0; 
		my $CallInNotExp = 0;
		my $Correct = 0;
		my $NotCorrect = 0;
		foreach $SNP (sort keys %inRegion) {
			if (!($SNPsToIgnore{$SNP})) {
				my $status = $name2statSam{$SNP};
				if ($SNPsInPath{$SNP}) {
					if ($status == 1) {
						$Correct++;
						print STATOUT $SNP."\tTP\n";
					} else {
						$NotCorrect++;
						print STATOUT $SNP."\tFN\n";
					} 
				} else {
					if ($status == 0 || $status == -1) {
						$Correct++;
						print STATOUT $SNP."\tTN\n";
					} else {
						$NotCorrect++;
						print STATOUT $SNP."\tFP\n";
					} 
				}
				my $origin = $name2origin{$SNP};
				if ($CallExp{$SNP} && $origin eq 'called') {
					$CallInExp++;
				} 
				if ($NotCallExp{$SNP} && $origin eq 'called') {
					$CallInNotExp++;
				} 
			}
		}
		my $score1 = $Correct/($Correct + $NotCorrect)*100;
		print OUT "ACC: $score1\n";
		my @CallExp;
		foreach $key (sort keys %CallExp) {
			push(@CallExp, $key);
		}	
		my $CallExp = $#CallExp + 1;
		my @NotCallExp;
		foreach $key (keys %NotCallExp) {
			push(@NotCallExp, $key);
		}
		my $NotCallExp = $#NotCallExp + 1;
		my $TP = $CallInExp;
		my $FP = $CallInNotExp;
		my $TN = $NotCallExp - $CallInNotExp;
		my $FN = $CallExp - $CallInExp;
		if (($TP+$FN) != 0) {
			my $score2 = $TP/($TP+$FN)*100;
			print OUT "SEN: $score2\n";
		} else {
			print OUT "SEN: NA\n";
		}
		if (($FP+$TN) != 0) {
			my $score3 = $TN/($FP+$TN)*100;
			print OUT "SPE: $score3\n";
		} else {
			print OUT "SPE: NA\n";
		}
		if (($TP+$FP) != 0) {
			$score5 = $TP/($TP+$FP)*100;
			print OUT "PRE: $score5\n";
		} else {
			print OUT "PRE: NA\n";
		}
		if (($TP+$FN) != 0) {
			$score6 = $TP/($TP+$FN)*100;
			print OUT "REC: $score6\n";
		} else {
			print OUT "REC: NA\n";
		}
		if (($TP+$FN) != 0 && ($TP+$FP) != 0) {
			my $score7 = ($score5 + $score6)/2;
			print OUT "Fsc: $score7\n";
		} else {
			print OUT "Fsc: NA\n";
		}
		if ((($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN)) != 0) {
			my $score4 = (($TP*$TN)-($FP*$FN))/sqrt(($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))*100;
			print OUT "MCC: $score4\n";
		} else {
			print OUT "MCC: NA\n";
		}
	} else {
		print OUT "When there is more than one haplogroup, SHELLY cannot determine the quality\n";
	}
	close OUT;
	close STATOUT;
}



#################################################################################
#				END OF SUBROUTINES				#
#################################################################################


##################### Create hash of quality control file. ######################

my %qc2snp;

open (QC, "< $qcFile") or die " Could not open QC file: $qcFile\n";
while (<QC>) {
	my $line = $_; 
	chomp($line);
	my ($name, @SNPs) = split(/\t/, $line); 
	my @defSNPs;
	foreach $snp (@SNPs) {
		if ($snp ne '') { push(@defSNPs, $snp);}
	}
	my $SNPs = join('_', @defSNPs);	
	$qc2snp{$name} = $SNPs;
}
close QC;


######################## Create hash of conversion file. ########################

my %name2anc;
my %name2mut;
my %name2ancestral;
my %name2mutant;
my %name2pos;

my %SNPsInConv;
my @SNPsInConv;
my %synonyms; 

my $flag = 0; 

open (CONV, "< $convFile") or die " Could not open conversion file: $convFile\n";
while (<CONV>) {
	my $line = $_; 
	$line =~ s/[\x0A\x0D]//g; 
	chomp($line);
	$line =~ s/\n//;
	if ($flag == 1) {
		my ($name, $pos36, $pos37, $mutation, $type, $ignore) = split(/\t/, $line); 
		if (trimSpace($ignore) eq 'no') {

			
			$mutation = trimSpace($mutation);
			my ($anc, $mut) = split(/->/, $mutation);
			
			$name = trimSpace($name);

			if ($version eq 'hg18') {
				push(@SNPsInConv, $name);
				$pos36 = trimSpace($pos36);
				$name2anc{$name} = $pos36."_".$anc;
				$name2mut{$name} = $pos36."_".$mut;
				my $nameMut = $pos36."_".$mut;
				$SNPsInConv{$nameMut} = 1;
				my $nameAnc = $pos36."_".$anc;
				$SNPsInConv{$nameAnc} = 1;
							
				my $synoName = $pos36."_".$mutation;
				my $prev = $synonyms{$synoName};
				my $next = $prev.$name."_";
				$synonyms{$synoName} = $next;

				$name2pos{$name} = $pos36;
				$name2ancestral{$name} = $anc;
				$name2mutant{$name} = $mut;
				
			} elsif ($version eq 'hg19') {
				push(@SNPsInConv, $name);
				$pos37 = trimSpace($pos37);
				$name2anc{$name} = $pos37."_".$anc;
				$name2mut{$name} = $pos37."_".$mut;			
				my $nameMut = $pos37."_".$mut;
				$SNPsInConv{$nameMut} = 1;
				my $nameAnc = $pos37."_".$anc;
				$SNPsInConv{$nameAnc} = 1;
												
				my $synoName = $pos37."_".$mutation;
				my $prev = $synonyms{$synoName};
				my $next = $prev.$name."_";
				$synonyms{$synoName} = $next;
				
				$name2pos{$name} = $pos37;
				$name2ancestral{$name} = $anc;
				$name2mutant{$name} = $mut;
			}
		}
	} else {
		$flag = 1;
	}
} 
close CONV;

my %reverseSynonyms;
foreach $synoName (keys %synonyms) {
	my $names =  $synonyms{$synoName};
	my @synos = split(/_/,$names);
	foreach $syno (@synos) {
		if ($syno ne '' && $syno ne '-') {
			$reverseSynonyms{$syno} = trimUnderscore($names);
		}
	}
}


############################ Create hash of tree file. ##########################
 
my %node2parent;
my %parent2nodes;
my %node2snp;
my @SNPsInTree;
my %SNPsInTree;
%altnames;

open (TREE, "< $treeFile") or die " Could not open tree file: $treeFile\n";
while (<TREE>) {
	my $line = $_; 
	chomp($line);
	my @splittedLine = split(/\t/, $line); 
	my $SNPs;
	foreach $elem (@splittedLine[3..$#splittedLine]) {
		if ($elem =~ m/\D*/ && $elem ne '-') {
			$SNPs = $SNPs.$elem."_";
			push(@SNPsInTree, $elem);
			my $posMut = $name2mut{$elem};
			my $posAnc = $name2anc{$elem};
			if ($posMut ne '') {
				$SNPsInTree{$posMut} = 1;
			} elsif ($posAnc ne '') {
				$SNPsInTree{$posAnc} = 1;
			} 
		}
	}
	my $node = $splittedLine[0];
	my $altName = $splittedLine[1];
	$altnames{$node} = $altName;
	my $parent = $splittedLine[2];
	$node2parent{$node} = $parent;
	my $children = $parent2nodes{$parent};
	$children = $children.$node."_";
	$parent2nodes{$parent} = $children;
	$node2snp{$node} = $SNPs;
} 
close TREE;

foreach $node (sort keys %node2snp) {
	my $snps = $node2snp{$node};
	$snps = trimUnderscore($snps);
	$node2snp{$node} = $snps;
}
 
foreach $parent (sort keys %parent2nodes) {
	my $nodes = $parent2nodes{$parent};
	$nodes = trimUnderscore($nodes);
	$parent2nodes{$parent} = $nodes;
}


############### Create array of SNPs in conversion File, not tree. ##############

@SNPsInConv = uniqueArray(@SNPsInConv);
@SNPsInTree = uniqueArray(@SNPsInTree);

my %tree;
foreach $elem (@SNPsInTree) {
	$elem = trimSpace($elem);
	$tree{$elem} = 1;
}
my @SNPsConvNotTree;
foreach $elem (@SNPsInConv) {
	$elem = trimSpace($elem);
	my @synonyms = split(/\_/, $reverseSynonyms{$elem});
	my $digit = 1;
	foreach $syno (@synonyms) {
		if (isInArray($syno, @SNPsInTree) == 1) {
			$digit = 0;
		}
	}
	if ($digit == 1) {
		push(@SNPsConvNotTree, $elem);
	}
}
# @SNPsConvNotTree contains elements from '@SNPsInConv' that are not in '@SNPsInTree'


######################## Get reference base for each SNP. #######################

# Get reference sequence.
my $refSeq;
my $in  = Bio::SeqIO->new(-file => "$hgFile" , -format => 'Fasta');
while(my $sobj = $in->next_seq) {
	my $id = $sobj->id;
	$refSeq = $sobj->seq();
}

# Get reference bases of SNP.
my %reference;
foreach $SNP (sort @SNPsInConv) {
	my $position = substr($name2mut{$SNP}, 0, -2) - 1;
	if ($position ne -1) {
		my $base = substr($refSeq, $position, 1);
		$base =~ tr/a-z/A-Z/;
		$position++;
		$reference{$SNP} = $position."_".$base;
	}
}

# Get reference bases of SNP (for quality determination after haplogroup determination).
my %referenceSNPs;
foreach $name (sort keys %name2pos) {
	my $position = $name2pos{$name} - 1;
	my $mutant = $name2mutant{$name};
	my $ancestral = $name2ancestral{$name};
	
	my $base = substr($refSeq, $position, 1);
	$base =~ tr/a-z/A-Z/;
	if ($base eq $mutant) {
		$referenceSNPs{$name} = 1;
	} elsif ($base eq $ancestral) {
		$referenceSNPs{$name} = 0;
	} else {
		$referenceSNPs{$name} = -1
	}
}

############################ Get leaves. ###########################

my %posLeaves;
my %negLeaves;
my @leaves; 

foreach $node (sort keys %node2parent) {
	my $children = $parent2nodes{$node};
	if ($children eq '') {
		my $SNPs = $node2snp{$node};
		$posLeaves{$node} = $SNPs;
		push(@leaves, $node);
	} else {
		my $negSNPs = '';
		foreach $child (split(/\_/, $children)) {
			my $SNPs = $node2snp{$child};
# 			$SNPs = substr($SNPs,0,-1);
			if ($SNPs ne '') {
				$negSNPs = $negSNPs.$SNPs.'*';
			}
		}

		my $negNode = $node."*";
		push(@leaves, $negNode);
		$negLeaves{$negNode} = $negSNPs;
		my $posSNPs = $node2snp{$node};
		if ($node ne 'Root') {
			while ($posSNPs eq '') {
				$node = $node2parent{$node};
				$posSNPs = $node2snp{$node};
			}
		}
		$posSNPs = $node2snp{$node};
		$posLeaves{$negNode} = $posSNPs;
	}
}


foreach $node (sort keys %posLeaves) {
	my $SNPs = $posLeaves{$node};
	$SNPs = trimUnderscore($SNPs);
	$posLeaves{$node} = $SNPs;
}

foreach $node (sort keys %negLeaves) {
	my $SNPs = $negLeaves{$node};
	$SNPs = trimStar($SNPs);
	$negLeaves{$node} = $SNPs;
}




################### Get called SNPs for input file. ##########################

# Crease hash with all SNPs in file. 
my %SNPsInFile;
my @newSNPs;

open (IN, "< $inFile") or die "Could not open input file: $inFile\n";
my $numberLines = 0;
while (<IN>) {
	my $line = $_; 
	chomp($line);
	my ($chr, $pos, $ref, $mut) = split(/\t/, $line);
 	if ($chr ne 'chrY' || $chr ne 'Y') {
		my $name = $pos.'_'.$mut;
		$SNPsInFile{$name} = 1;
		if ($SNPsInConv{$name} != 1) {
			push(@newSNPs, $name);
		}
 	}
	if ($chr eq 'chrY' || $chr eq 'Y') {$numberLines++;};
}
close IN;

######################## Get results for input file. #######################

if ($numberLines >= 0) {

	# Set path to results output file.
	my @inFile = split(/\//, $inFile);
	$outfile = $outDir.substr($inFile[-1], 0, -4).'_analysis.txt';
	open (OUT, "> $outfile");
	my @linesToPrint;

	my $refFile;
	$refFile = $stathgFile; 

	my $analysisFile = $outfile;
	my $statusFile = substr($analysisFile,0, -12)."status.txt";
	my $outFile = substr($analysisFile,0, -12)."quality.txt";
	my $statusOutFile = substr($analysisFile,0, -12)."statusSNPs.txt";

	if ($filePartialOption ne '') {
		my %partial;
		open (IN, "< $filePartialOption") ;
		while(<IN>) {
			my $line = $_;
			chomp($line);
			my (undef, $start,$end) = split(/\t/, $line);
			$partial{$start} = $end;
		}
		close IN;

		# Create hash for results. 
		my %results;
		my %called;
		my %origins;
		foreach $SNP (@SNPsInConv) {

			if ($SNPsInFile{$name2anc{$SNP}}) {
				$results{$SNP} = 0;
				$called{$SNP} = 0;
				$origins{$SNP} = 'called';
			} elsif ($SNPsInFile{$name2mut{$SNP}}) {
				$results{$SNP} = 1;
				$called{$SNP} = 1;
				$origins{$SNP} = 'called';
			} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
				$results{$SNP} = 0;
				$origins{$SNP} = 'reference';
			} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
				$results{$SNP} = 1;
				$origins{$SNP} = 'reference';
			} else {
				$results{$SNP} = -1;
				$origins{$SNP} = 'unknown';
			}
		}

		# Quality Control. 
		$cutoff = $cutoffLowQ;
		# Create hash of definite results. 
		my %definiteResults;
		my %inRegion;
		foreach $SNP (@SNPsInTree) {
			my ($position, undef) = split(/_/,$name2anc{$SNP});
			my $inRegion = 0;

			foreach $start (sort keys %partial) {
				my $end = $partial{$start};
				if ($position >= $start && $position <= $end) {
					$inRegion = 1 ;
				}
			}
			if ($inRegion == 1){ 
				$definiteResults{$SNP} = $results{$SNP};
				$inRegion{$SNP} = 1 ;
			} elsif ($inRegion == 0) {
				if ($SNPsInFile{$name2anc{$SNP}}) {
					$definiteResults {$SNP} = 0;
				} elsif ($SNPsInFile{$name2mut{$SNP}}) {
					$definiteResults {$SNP} = 1;
				} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
					$definiteResults {$SNP} = 0;
				} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
					$definiteResults {$SNP} = 0;
				} else {
					$definiteResults {$SNP} = -1;
				}
			}
		}

		my $statusfile = $outDir.substr($inFile[-1], 0, -4).'_status.txt';
		
		open (STATUS, "> $statusfile");
		foreach $key (keys %results) {
			print STATUS $key."\t".$results{$key}."\t".$origins{$key}."\n";
		}
		close STATUS;
		
		#########################################################################
		push(@linesToPrint,"QUALITY\n");
		push(@linesToPrint,"-------\n");
		push(@linesToPrint,"Option Partial (insufficient)\n\n");

		# Set path of SNPsConvNotTree output file. 
		my $CNTfile = $outDir.substr($inFile[-1], 0, -4).'_ConvNotTree.txt';
		open (CNT, "> $CNTfile") ;

		# Set path of SNPsNotConvNotTree output file. 
		my $NCNTfile = $outDir.substr($inFile[-1], 0, -4).'_newSNPs.txt';
		open (NCNT, "> $NCNTfile") ;
		
		# NotConvNotTree
		foreach $snp (sort @newSNPs) {
			print NCNT $snp."\n";
		}
		# ConvNotTree
		my @SNPsCNT;
		
		foreach $SNPname (sort @SNPsConvNotTree) {
			if ($SNPname ne '' && $SNPname ne '-') {
				if ($SNPsInFile{$name2anc{$SNPname}}) {
					# nothing
				} elsif ($SNPsInFile{$name2mut{$SNPname}}) {
					push(@SNPsCNT, $reverseSynonyms{$SNPname});
					$origins{$reverseSynonyms{$SNPname}} = 'called';
				} elsif ($reference{$SNP} eq $name2anc{$SNPname}) {
					# nothing
				} elsif ($reference{$SNP} eq $name2mut{$SNPname}) {
					push(@SNPsCNT, $reverseSynonyms{$SNPname});
					$origins{$reverseSynonyms{$SNPname}} = 'reference';
				} else {
					# nothing
				}
			}
		}

		foreach $SNPname (uniqueArray(@SNPsCNT)) {
			my $ori = $origins{$SNPname};
			my @SNPname = split(/_/, $SNPname);
			$SNPname =~ s/_/ - /g;
			
			my $id;
			if ($ori eq 'called') {
				$id = $name2pos{$SNPname[0]}.'_'.$name2mutant{$SNPname[0]};
			} else {
				$id = $name2pos{$SNPname[0]}.'_'.$name2ancestral{$SNPname[0]};
			}
			print CNT $id."\t".$SNPname."\t".$ori."\n";
		}
			
		close NCNT;
		close CNT;

		# Vertical method. 
		push(@linesToPrint,"VERTICAL\n");
		push(@linesToPrint,"--------\n");

		my @verticals;
		foreach $leaf (sort keys %posLeaves) {
			$leaf = trimSpace($leaf);
			if (substr($leaf, -1) ne '*') {
				my $SNPs = $posLeaves{$leaf};
				if ($SNPs ne '') { # don't include positive leaves that don't have SNPs
					my @SNPs = split(/\_/, $SNPs);
					my @resultsSNPs;
					foreach $SNP (@SNPs) {
						push(@resultsSNPs, $definiteResults{$SNP});
					}
					if ( getMetric(@resultsSNPs) > $cutoff ) {
						my $name = layout($leaf);
						$name = altName($name);
						push(@linesToPrint, "\$ $name\n");
						push(@verticals, $leaf);
					}
				}
			}
		}

		foreach $leaf (sort keys %negLeaves) {
			my $children = $negLeaves{$leaf};
			if ($children ne '') { # when some children have SNPs, only include leave when it is also positive && SNPs of children are too
				my @children = split(/\*/, $children);
				my @resChildren;
				foreach $child (@children) {
					my @SNPs = split(/\_/, $child);
					my @resultsSNPs;
					foreach $SNP (@SNPs) {
						push(@resultsSNPs, $definiteResults{"$SNP"});
					}
					push(@resChildren, getMetric(@resultsSNPs)) ;
				}

				my $res = 1;
				foreach $resChild (@resChildren) {
					if ($resChild > $cutoff) {
						$res = 0;
					}
				}

				my $currentSNPs = $posLeaves{$leaf};
				if ($currentSNPs ne '') {
					my @currentSNPs = split(/\_/, $currentSNPs);
					my @resCurrent;
					foreach $SNP (@currentSNPs) {
						push(@resCurrent, $definiteResults{$SNP});
					}
					if (getMetric(@resCurrent) < $cutoff) {
						$res = 0;
					}
				}

				if ($res == 1) {
					my $name = layout($leaf);
					$name = altName($name);
					push(@linesToPrint, "\$ $name\n");
					push(@verticals, $leaf);
				}
			} else {  # when all children don't have SNPs, only include leaf when it is also positive
				my $res = 1;
			
				my $currentSNPs = $posLeaves{$leaf};
				if ($currentSNPs ne '') {
					my @currentSNPs = split(/\_/, $currentSNPs);
					my @resCurrent;
					foreach $SNP (@currentSNPs) {
						push(@resCurrent, $definiteResults{$SNP});
					}
					if (getMetric(@resCurrent) < $cutoff) {
						$res = 0;
					}
				}

				if ($res == 1) {
					my $name = layout($leaf);
					$name = altName($name);
					push(@linesToPrint, "\$ $name\n");
					push(@verticals, $leaf);
				}	
			}			
		}

		# Horizontal method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "HORIZONTAL\n");
		push(@linesToPrint, "----------\n");
		my @toInvest = 'Root';
		my @horizontals;

		my %nodeScores;
		while(@toInvest) {
			my @temp;
			foreach $inv (@toInvest) {
				my @tempInv;
				my $children = $parent2nodes{$inv};
				foreach $child (split(/\_/, $children)) {
					my $SNPs = $node2snp{$child};
					if ($SNPs ne '') {
						my @SNPs = split(/\_/, $SNPs);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							if (substr($SNP, -1) eq "\r") { $SNP = substr($SNP, 0,-1) }
							push(@resultsSNPs, $definiteResults{$SNP});
						}
						if ( getMetric(@resultsSNPs) > $cutoff) {
							push(@tempInv, $child);
						}
						$nodeScores{$child} = getMetric(@resultsSNPs);
					} else {
						push(@tempInv, $child);
					}	
				}
				if (@tempInv) {
					foreach $elem (@tempInv) {
						push(@temp, $elem);
					}
				} else {
					push(@horizontals, $inv);
				}
			}
			@toInvest = @temp;
		}


		foreach $HR (@horizontals) {
			my $parents;
			my $parent = $node2parent{$HR};
			while ($parent ne '-') {
				$parents = $parents.$parent." ";
				$parent = $node2parent{$parent};
			}
			my $line = "! ".$HR."\t".$parents."\n";
			push(@linesToPrint, $line);
		}

		# Combi method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "COMBI\n");
		push(@linesToPrint, "-----\n");
		my @combis;
		if ($#verticals == -1) {
			@combis = @horizontals;
		} else {
			foreach $vert (@verticals) {
				my @path;
				if (substr($vert, -1) ne '*') {
					my $result = 0;
					my $parent = $vert;
					push(@path, $parent);
					while ($parent ne '-' && $parent ne 'Root') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
					foreach $hori (@horizontals) {
						if (isInArray($hori, @path)) {
							$result = 1;
						}
					}
					if ($result == 1) {
						push(@combis, $vert);
					}
				} else {
					my $result = 0;
					$vertWithoutStar = substr($vert, 0, -1);
					my $parent = $vertWithoutStar;
					push(@path, $parent);
					while ($parent ne '-' && $parent ne 'Root') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
					foreach $hori (@horizontals) {
						if (isInArray($hori, @path)) {
							$result = 1;
						}
					}
					if ($result == 1) {
						push(@combis, $vert);
					}
				}
			}
		}
		foreach $combi (@combis) {
			my $name = layout($combi);
			$name = altName($name);
			push(@linesToPrint, "% $name\n");
		}

		# Most specific method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "SPECIFIC\n");
		push(@linesToPrint, "--------\n");
		my %specific;
		foreach $combi (@combis) {
			my @path;
			my $result = 0;
			if (substr($combi, -1) ne '*') {
				my $parent = $combi;
				while ($parent ne '-') {
					$parent = $node2parent{$parent};
					push(@path, $parent);
				}

			} else {
				$combiWithoutStar = substr($combi, 0, -1);
				my $parent = $combiWithoutStar;
				while ($parent ne '-') {
					$parent = $node2parent{$parent};
					push(@path, $parent);
				}
			}
			my $hits = 0;
			foreach $comb (sort @combis) {

				if (substr($comb, -1) eq '*') {
					$combWithoutStar = substr($comb, 0,-1) ;
				} else {
					$combWithoutStar = $comb; 
				}
				if (isInArray($combWithoutStar, @path)) {
					$hits++;
				}
			}
			if ($hits > 0) {
				$specific{$combi} = $hits;
			}
		}
		my $prev;
		my @best;
		my $best = 1;
		foreach $spec (sort { $specific{$b} cmp $specific{$a} } keys %specific) {
			if ($best == 1) {
				push(@best, $spec);
				$best = $specific{$spec};
			} else {
				if ($best eq $specific{$spec}) {
					push(@best, $spec);
				}
			}
			$prev = $specific{$spec};
			my $name = layout($spec);
			$name = altName($name);
			my $line = "@ ".$name."\n"; 
			push(@linesToPrint, $line);
		}

		# Results. 
		push(@linesToPrint, "\n");
		push(@linesToPrint, "RESULTS\n");
		push(@linesToPrint, "-------\n");
		my @finalresults;
		if ($#best == -1) {
			foreach $res (@combis) {
				my $name = layout($res);
				$name = altName($name);
				$result = $name." ".$verticals{$res};
				push(@finalresults, $result);
			}

			if ($#combis == 0) {
				my $name = layout($combis[0]);
				$name = altName($name);

				push(@linesToPrint, "> $name\n");
			} else {
				my $bestGroup;
				my $bestGroupVerticalScore = 0;
							
				foreach $value (@combis) {
					my $name = layout($value);
					$name = altName($name);
					push(@linesToPrint, "> $name\n");
				}
			}
		} else {
			foreach $res (@best) {
				my $name = layout($res);
				$name = altName($name);
				$result = $name." ".$verticals{$res};
				push(@finalresults, $result);
			}
			if ($#best == 0) {
				my $name = layout($best[0]);
				$name = altName($name);
				push(@linesToPrint, "> $name\n");
			} else {
				my $bestGroup;
				my $bestGroupVerticalScore = 0;
				
				my @results;			
				foreach $res (@best) {
					push(@results, $res);
				}
				if ($#results == -1) {
					@results = @best;
				}
				
				foreach $res (@results) {
					my $name = layout($res);
					$name = altName($name);

					push(@linesToPrint, "> $name\n");
				}
			}
		}

		# Calculate quality.
		my $mcc = 100;
		my $runagain = 0;
		if ($#finalresults > 0) {
			$runagain = 1;
		} else {
			my $haplogroup;
			if ($#finalresults == 0) {
				$haplogroup = $finalresults[0];
			}
			# get SNPs on its path
			$haplogroup  =~ s/\[.*\]//;
			$haplogroup = substr($haplogroup, 0, -1);
			my $name = $haplogroup;
			my %SNPsInPath;
			my %SNPsToIgnore;

			if ( $name =~ m/\*/ ) {
				while($name =~ m/\*/) {
					$name =~ s/\*//;
				}
				if (substr($name, -1) eq ' ') {
					$name = substr($name, 0, -1);
				}
				my $parent = $node2parent{$name};
				my @children = split(/\_/, $parent2nodes{$name});

				# To ignore
				my @childrenToIgnore;
				while (@children) {
					my $child = $children[0];
					if ($child ne $name) {
						push(@childrenToIgnore, $child);
						my @grandchildren = split(/\_/, $parent2nodes{$child}); 
						foreach $grandchild (@grandchildren) {
							push(@children, $grandchild);
						}
					}
					shift(@children);
				}
				@children = split(/\_/, $parent2nodes{$parent});
				# To ignore
				while (@children) {
					my $child = $children[0];
					if ($child ne $name) {
						push(@childrenToIgnore, $child);
						my @grandchildren = split(/\_/, $parent2nodes{$child}); 
						foreach $grandchild (@grandchildren) {
							push(@children, $grandchild);
						}
					}
					shift(@children);
				}

				foreach $childToIgnore (@childrenToIgnore) {
					my @SNPs = split(/\_/, $node2snp{$childToIgnore}); 
					foreach $SNP (@SNPs) {
						$SNPsToIgnore{$SNP} = 1;
					}
				}

				# In Path
				my @SNPs = split(/\_/, $node2snp{$name});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}

				while($parent ne 'Root') {
					my @SNPs = split(/\_/, $node2snp{$parent});
					foreach $SNP (@SNPs) {
						$SNPsInPath{$SNP} = 1;
					}
					$parent = $node2parent{$parent};
				}
			} else {
				# In Path
				while(substr($name, -1) eq ' ') {
					$name = substr($name, 0, -1);
				}
				my $parent = $node2parent{$name};
				my @SNPs = split(/\_/, $node2snp{$name});
				foreach $SNP (@SNPs) {
					$SNPsInPath{$SNP} = 1;
				}
				while($parent ne 'Root') {
					my @SNPs = split(/\_/, $node2snp{$parent});
					foreach $SNP (@SNPs) {
						$SNPsInPath{$SNP} = 1;
					}
					$parent = $node2parent{$parent};
				}
			}

			# Create hash with expected called SNPs. 
			my %CallExp;
			my %NotCallExp;
			foreach $SNP (sort keys %inRegion) {
				if (!($SNPsToIgnore{$SNP})) {
					my $statRef = $referenceSNPs{$SNP};
					if ($SNPsInPath{$SNP}) {
						if ($statRef != 1) {
							$CallExp{$SNP} = 1;
						} else {
							$NotCallExp{$SNP} = 1;
						}
					} else {
						if ($statRef != 0) {
							$CallExp{$SNP} = 1;
						} else {
							$NotCallExp{$SNP} = 1;
						}
					}
				}
			}

			# Get qualities. 
			my $CallInExp = 0; 
			my $CallInNotExp = 0;
			my $Correct = 0;
			my $NotCorrect = 0;

			foreach $SNP (sort keys %inRegion) {
				if (!($SNPsToIgnore{$SNP})) {
					my $status = $results{$SNP};
					if ($SNPsInPath{$SNP}) {
						if ($status == 1) {
							$Correct++;
						} else {
							$NotCorrect++;
						} 
					} else {
						if ($status == 0 || $status == -1) {
							$Correct++;
						} else {
							$NotCorrect++;
						} 
					}
					my $origin = $origins{$SNP};
					if ($CallExp{$SNP} && $origin eq 'called') {
						$CallInExp++;
					} 

					if ($NotCallExp{$SNP} && $origin eq 'called') {
						$CallInNotExp++;
					} 
				}
			}

			my @CallExp;
			foreach $key (sort keys %CallExp) {
				push(@CallExp, $key);
			}	
			my $CallExp = $#CallExp + 1;
					
			my @NotCallExp;
			foreach $key (keys %NotCallExp) {
				push(@NotCallExp, $key);
			}
			my $NotCallExp = $#NotCallExp + 1;

			my $TP = $CallInExp;
			my $FP = $CallInNotExp;
			my $TN = $NotCallExp - $CallInNotExp;
			my $FN = $CallExp - $CallInExp;

			if ((($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN)) != 0) {
				$mcc = (($TP*$TN)-($FP*$FN))/sqrt(($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))*100;
			} else {
				$mcc = 'NA';
			}

			if ($mcc ne 'NA' && $mcc < 95) {
				$runagain = 1;
				push(@linesToPrint,"\nMCC\n");
				push(@linesToPrint,"---\n");
				push(@linesToPrint,"$mcc\n");
			}
		}

		if ($runagain == 1) {
			# Run again as insufficient (old method)
			foreach $lineToPrint (@linesToPrint) {
				print OUT "### $lineToPrint";
			}

			print OUT "\n\nRESULTS FOR INSUFFICIENT CALL QUALITY BASED ON MCC\n";
			print OUT "--------------------------------------------------\n\n";
			foreach $line (@linesToPrint[0..2]) {
				$line =~ s/^ //;
				$line =~ s/ $//;
				print OUT $line;
			}
			print OUT "0 \(insufficient\)\n\n";

			$QCtype = 'insufficient';
			$cutoff = $cutoffLowQ;
			foreach $SNP (@SNPsInConv) {
				if ($SNPsInFile{$name2anc{$SNP}}) {
					$definiteResults{$SNP} = 0;
				} elsif ($SNPsInFile{$name2mut{$SNP}}) {
					$definiteResults{$SNP} = 1;
				} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
					$definiteResults{$SNP} = 0;
				} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
					$definiteResults{$SNP} = 0;
				} else {
					$definiteResults{$SNP} = -1;
				}
			}

			# Vertical method. 
			print OUT "VERTICAL\n--------\n"; $| = 1;
			my @verticals;
			foreach $leaf (sort keys %posLeaves) {
				$leaf = trimSpace($leaf);
				if (substr($leaf, -1) ne '*') {
					my $SNPs = $posLeaves{$leaf};
					if ($SNPs ne '') { # don't include positive leaves that don't have SNPs
						my @SNPs = split(/\_/, $SNPs);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							push(@resultsSNPs, $definiteResults{$SNP});
						}

						if ( getMetric(@resultsSNPs) > $cutoff ) {
							my $name = layout($leaf);
							$name = altName($name);
							print OUT "\$ $name\n"; $| = 1;
							push(@verticals, $leaf);
						}
					}
				}
			}

			foreach $leaf (sort keys %negLeaves) {
				my $children = $negLeaves{$leaf};
				if ($children ne '') { # when some children have SNPs, only include leave when it is also positive && SNPs of children are too
					my @children = split(/\*/, $children);
					my @resChildren;
					foreach $child (@children) {
						my @SNPs = split(/\_/, $child);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							push(@resultsSNPs, $definiteResults{"$SNP"});
						}
						push(@resChildren, getMetric(@resultsSNPs)) ;
					}

					my $res = 1;
					foreach $resChild (@resChildren) {
						if ($resChild > $cutoff) {
							$res = 0;
						}
					}
				
					my $currentSNPs = $posLeaves{$leaf};
					if ($currentSNPs ne '') {
						my @currentSNPs = split(/\_/, $currentSNPs);
						my @resCurrent;
						foreach $SNP (@currentSNPs) {
							push(@resCurrent, $definiteResults{$SNP});
						}
						if (getMetric(@resCurrent) < $cutoff) {
							$res = 0;
						}
					}

					if ($res == 1) {
						my $name = layout($leaf);
						$name = altName($name);
						print OUT "\$ $name\n"; $| = 1;
						push(@verticals, $leaf);
					}
				} else {  # when all children don't have SNPs, only include leaf when it is also positive
					my $res = 1;
				
					my $currentSNPs = $posLeaves{$leaf};
					if ($currentSNPs ne '') {
						my @currentSNPs = split(/\_/, $currentSNPs);
						my @resCurrent;
						foreach $SNP (@currentSNPs) {
							push(@resCurrent, $definiteResults{$SNP});
						}
						if (getMetric(@resCurrent) < $cutoff) {
							$res = 0;
						}
					}

					if ($res == 1) {
						my $name = layout($leaf);
						$name = altName($name);
						print OUT "\$ $name\n"; $| = 1;
						push(@verticals, $leaf);
					}	
				}			
			}


			# Horizontal method.
			print OUT "\nHORIZONTAL\n----------\n"; $| = 1;
			my @toInvest = 'Root';
			my @horizontals;

			my %nodeScores;
			while(@toInvest) {
				my @temp;
				foreach $inv (@toInvest) {
					my @tempInv;
					my $children = $parent2nodes{$inv};
					foreach $child (split(/\_/, $children)) {
						my $SNPs = $node2snp{$child};
						if ($SNPs ne '') {
							my @SNPs = split(/\_/, $SNPs);
							my @resultsSNPs;
							foreach $SNP (@SNPs) {
								if (substr($SNP, -1) eq "\r") { $SNP = substr($SNP, 0,-1) }
								push(@resultsSNPs, $definiteResults{$SNP});
							}
							if ( getMetric(@resultsSNPs) > $cutoff) {
								push(@tempInv, $child);
							}
							$nodeScores{$child} = getMetric(@resultsSNPs);
						} else {
							push(@tempInv, $child);
						}	
					}
					if (@tempInv) {
						foreach $elem (@tempInv) {
							push(@temp, $elem);
						}
					} else {
						push(@horizontals, $inv);
					}
				}
				@toInvest = @temp;
			}


			foreach $HR (@horizontals) {
				my $parents;
				my $parent = $node2parent{$HR};
				while ($parent ne '-') {
					$parents = $parents.$parent." ";
					$parent = $node2parent{$parent};
				}
				print OUT "! ".$HR."\t".$parents."\n"; $| = 1;
			}

			# Combi method.
			print OUT "\nCOMBI\n-----\n"; $| = 1;
			my @combis;
			if ($#verticals == -1) {
				@combis = @horizontals;
			} else {
				foreach $vert (@verticals) {
					my @path;
					if (substr($vert, -1) ne '*') {
						my $result = 0;
						my $parent = $vert;
						push(@path, $parent);
						while ($parent ne '-' && $parent ne 'Root') {
							$parent = $node2parent{$parent};
							push(@path, $parent);
						}
						foreach $hori (@horizontals) {
							if (isInArray($hori, @path)) {
								$result = 1;
							}
						}
						if ($result == 1) {
							push(@combis, $vert);
						}
					} else {
						my $result = 0;
						$vertWithoutStar = substr($vert, 0, -1);
						my $parent = $vertWithoutStar;
						push(@path, $parent);
						while ($parent ne '-') {
							$parent = $node2parent{$parent};
							push(@path, $parent);
						}
						foreach $hori (@horizontals) {
							if (isInArray($hori, @path)) {
								$result = 1;
							}
						}
						if ($result == 1) {
							push(@combis, $vert);
						}
					}
				}
			}
			foreach $combi (@combis) {
				my $name = layout($combi);
				$name = altName($name);
				print OUT "% ".$name."\n"; $| = 1;
			}

			# Most specific method.
			print OUT "\nSPECIFIC\n--------\n"; $| = 1;
			my %specific;
			foreach $combi (@combis) {
				my @path;
				my $result = 0;
				if (substr($combi, -1) ne '*') {
					my $parent = $combi;
					while ($parent ne '-') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}

				} else {
					$combiWithoutStar = substr($combi, 0, -1);
					my $parent = $combiWithoutStar;
					while ($parent ne '-') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
				}
				my $hits = 0;
				foreach $comb (sort @combis) {

					if (substr($comb, -1) eq '*') {
						$combWithoutStar = substr($comb, 0,-1) ;
					} else {
						$combWithoutStar = $comb; 
					}
					if (isInArray($combWithoutStar, @path)) {
						$hits++;
					}
				}
				if ($hits > 0) {
					$specific{$combi} = $hits;
				}
			}
			my $prev;
			my @best;
			my $best = 1;
			foreach $spec (sort { $specific{$b} cmp $specific{$a} } keys %specific) {
				if ($best == 1) {
					push(@best, $spec);
					$best = $specific{$spec};
				} else {
					if ($best eq $specific{$spec}) {
						push(@best, $spec);
					}
				}
				$prev = $specific{$spec};
				my $name = layout($spec);
				$name = altName($name);
				print OUT "@ ".$name."\n"; $| = 1;
			}

			# Results. 
			print OUT "\nRESULTS\n-------\n"; $| = 1;
			if ($#best == -1) {
				if ($#combis == 0) {
					my $name = layout($combis[0]);
					$name = altName($name);
					print OUT "> ".$name."\n"; $| = 1;
				} else {
					my $bestGroup;
					my $bestGroupVerticalScore = 0;
								
					foreach $value (@combis) {
						my $name = layout($value);
						$name = altName($name);
						print OUT "> ".$name." ".$verticals{$res}."\n"; $| = 1;
					}
				}
			} else {
				if ($#best == 0) {
					my $name = layout($best[0]);
					$name = altName($name);
					print OUT "> ".$name."\n"; $| = 1;
				} else {
					my $bestGroup;
					my $bestGroupVerticalScore = 0;
					
					my @results;			
					foreach $res (@best) {
						push(@results, $res);
					}
					if ($#results == -1) {
						@results = @best;
					}
					
					foreach $res (@results) {
						my $name = layout($res);
						$name = altName($name);
						print OUT "> ".$name." ".$verticals{$res}."\n"; $| = 1;
					}
				}
			}
		} else {
			foreach $lineToPrint (@linesToPrint) {
				print OUT "$lineToPrint";
			}
		}
		
		close OUT;

		SHELLY2($analysisFile, $statusFile, $treeFile, $refFile, $filePartialOption, $convFile, $version, $outFile, $statusOutFile);
	 } else {
		# Create hash for results. 
		my %results;
		my %called;
		my %origins;
		foreach $SNP (@SNPsInConv) {

			if ($SNPsInFile{$name2anc{$SNP}}) {
				$results{$SNP} = 0;
				$called{$SNP} = 0;
				$origins{$SNP} = 'called';
			} elsif ($SNPsInFile{$name2mut{$SNP}}) {
				$results{$SNP} = 1;
				$called{$SNP} = 1;
				$origins{$SNP} = 'called';
			} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
				$results{$SNP} = 0;
				$origins{$SNP} = 'reference';
			} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
				$results{$SNP} = 1;
				$origins{$SNP} = 'reference';
			} else {
				$results{$SNP} = -1;
				$origins{$SNP} = 'unknown';
			}
		}

		my $statusfile = $outDir.substr($inFile[-1], 0, -4).'_status.txt';

		open (STATUS, "> $statusfile");
		foreach $key (keys %results) {
			print STATUS $key."\t".$results{$key}."\t".$origins{$key}."\n";
		}
		close STATUS;
		
		# Quality Control. 
		my @group1;
		my @allSNPs1 = split(/_/, $qc2snp{'group1'});
		my @SNPs1;
		foreach $SNP (@allSNPs1) {
			if (length($SNP) > 1) {push(@SNPs1, $SNP)}
		}			

		foreach $SNP1 (@SNPs1) {
			my $score = $results{$SNP1};
			push(@group1, $score);
		}
		
		my @group2; 
		my @allSNPs2 = split(/_/, $qc2snp{'group2'});			
		my @SNPs2;			
		foreach $SNP (@allSNPs2) {
			if (length($SNP) > 1) {push(@SNPs2, $SNP)}
		}	
		foreach $SNP2 (@SNPs2) {
			my $score = $results{$SNP2};
			push(@group2, $score);
		}

		my @R1a1a; 
		my @SNPsR1a1a = split(/_/, $qc2snp{'R1a1a'});			
		foreach $SNPR1a1a (@SNPsR1a1a) {
			my $score = $results{$SNPR1a1a};
			push(@R1a1a, $score);
		}

		my @extra1; 
		my @SNPSextra1 = split(/_/, $qc2snp{'extra1'});			
		foreach $SNPextra1 (@SNPSextra1) {
			my $score = $results{$SNPextra1};
			push(@extra1, $score);
		}

		my @extra2; 
		my @SNPSextra2 = split(/_/, $qc2snp{'extra2'});			
		foreach $SNPextra2 (@SNPSextra2) {
			my $score = $results{$SNPextra2};
			push(@extra2, $score);
		}

		my $QCscore;
		my $bestHaplo; 

		my $bestScore = 0;
		my @ones;
		my @f;

		foreach $qc (sort keys %qc2snp) {
			my @allSNPs = split(/_/, $qc2snp{$qc});
			my @SNPs;
			foreach $SNP (@allSNPs) {
				if (length($SNP) > 1) {push(@SNPs, $SNP)}
			}

			my @scores;
			foreach $SNP (@SNPs) {
				my $score = $results{$SNP};
				push(@scores, $score);
			}
			if (substr($qc, 0,5) ne 'group' && substr($qc, 0,5) ne 'extra' && $qc ne 'R1a1a') {
				my $newScore = getMetric(@scores);
				if ($newScore > $bestScore) {
					$bestScore = $newScore;
					$bestHaplo = $qc;
					@ones = @SNPs;
				}
				if ($qc eq 'F') {
					@f = @SNPs;
				}
			}
		}
		if ($bestHaplo ne 'F') {

			my @all; 
			foreach $SNP (@f) {
				push(@all, $SNP);
			}
			foreach $SNP (@ones) {
				push(@all, $SNP);
			}
			@all = uniqueArray(@all);
			my @notones = getDiffArrays(join('_', @all), join('_', @ones));
			my $total = $#all + 1;

			my $correct;
			foreach $snp (@ones) {
				my $score = $results{$snp};
				if ($score == 1) {
					$correct++;
				}
			}
			foreach $snp (@f) {
				my $score = $results{$snp};
				if ($score != 1) {
					$correct++;
				}
			}

			$QCscore = $correct / $total;
		} else { 

			my $bestScoreF = 0;
			my $bestHaploF; 

			# Tommy1
			my $correct1;
			my $all1;

			foreach $score (@group1) {
				$all1++;
				if ($score == 1) {
					$correct1++;
				}
			}
			foreach $score (@group2) {
				$all1++;
				if ($score == 0) {
					$correct1++;
				}
			}
			my $tommy1 = $correct1/$all1;

			if ($tommy1 > $bestScoreF) {
				$bestScoreF = $tommy1;
				$bestHaploF = 'G';
			}

			# Tommy2
			my $correct2;
			my $all2;
			foreach $score (@group1) {
				$all2++;
				if ($score == 0) {
					$correct2++;
				}
			}
			foreach $score (@group2) {
				$all2++;
				if ($score == 1) {
					$correct2++;
				}
			}
			my $tommy2 = $correct2/$all2;	

			if ($tommy2 > $bestScoreF) {
				$bestScoreF = $tommy2;
				$bestHaploF = 'R';
			}

			# Tommy3
			my $correct3;
			my $all3;
			foreach $score (@group1) {
				$all3++;
				if ($score == 0) {
					$correct3++;
				}
			}
			foreach $score (@group2) {
				$all3++;
				if ($score == 0) {
					$correct3++;
				}
			}
			my $tommy3 = $correct3/$all3;				
						
			if ($tommy3 > $bestScoreF) {
				$bestScoreF = $tommy3;
				$bestHaploF = 'not R & not G';
			}
			if ($bestHaploF eq 'tommy2' && $bestScoreF > 0.90) {
				foreach $snp (@group2) {
					if ($snp eq '0') {
						$bestScoreF = 0.0;
					}
				}
				foreach $snp (@R1a1a) {
					if ($snp eq '1') {
						$bestScoreF = 0.0;
					}
				}
			}
			$QCscore = $bestScoreF;
			$bestHaplo = $bestHaplo.": ".$bestHaploF;
		}
		my $QCtype = 'sufficient';

		#########################################################################
		my %definiteResults;
		my $cutoff;

		if ($QCscore < 0.90) {
			$QCtype = 'insufficient';
			$cutoff = $cutoffLowQ;
			foreach $SNP (@SNPsInConv) {
				if ($SNPsInFile{$name2anc{$SNP}}) {
					$definiteResults {$SNP} = 0;
				} elsif ($SNPsInFile{$name2mut{$SNP}}) {
					$definiteResults {$SNP} = 1;
				} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
					$definiteResults {$SNP} = 0;
				} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
					$definiteResults {$SNP} = 0;
				} else {
					$definiteResults {$SNP} = -1;
				}
			}
		} else {
			my $extra1 = 0;
			foreach $snp (@extra1) {
				if ($snp eq '1') {
					$extra1++;
				}
			}
			my $extra2 = 0;
			foreach $snp (@extra2) {
				if ($snp eq '1') {
					$extra2++;
				}
			}

			if ($extra1 == (length(@extra1)+1) or $extra2 == (length(@extra2)+1)) {
				$QCtype = 'insufficient';
				$cutoff = $cutoffLowQ;
				foreach $SNP (@SNPsInConv) {
					if ($SNPsInFile{$name2anc{$SNP}}) {
						$definiteResults {$SNP} = 0;
					} elsif ($SNPsInFile{$name2mut{$SNP}}) {
						$definiteResults {$SNP} = 1;
					} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
						$definiteResults {$SNP} = 0;
					} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
						$definiteResults {$SNP} = 0;
					} else {
						$definiteResults {$SNP} = -1;
					}
				}
			} else {
				$cutoff = $cutoffHighQ;
				%definiteResults = %results;
			}
		}
		#########################################################################
		push(@linesToPrint,"QUALITY\n");
		push(@linesToPrint,"-------\n");
		push(@linesToPrint,"$bestHaplo\n");
		push(@linesToPrint,"$QCscore \($QCtype\)\n");
		push(@linesToPrint,"\n");

		# Set path of SNPsConvNotTree output file. 
		my $CNTfile = $outDir.substr($inFile[-1], 0, -4).'_ConvNotTree.txt';
		open (CNT, "> $CNTfile") ;

		# Set path of SNPsNotConvNotTree output file. 
		my $NCNTfile = $outDir.substr($inFile[-1], 0, -4).'_newSNPs.txt';
		open (NCNT, "> $NCNTfile") ;
		
		# NotConvNotTree
		foreach $snp (sort @newSNPs) {
			print NCNT $snp."\n";
		}
		# ConvNotTree
		my @SNPsCNT;
		
		foreach $SNPname (sort @SNPsConvNotTree) {
			if ($SNPname ne '' && $SNPname ne '-') {
				if ($SNPsInFile{$name2anc{$SNPname}}) {
					# nothing
				} elsif ($SNPsInFile{$name2mut{$SNPname}}) {
					push(@SNPsCNT, $reverseSynonyms{$SNPname});
					$origins{$reverseSynonyms{$SNPname}} = 'called';
				} elsif ($reference{$SNP} eq $name2anc{$SNPname}) {
					# nothing
				} elsif ($reference{$SNP} eq $name2mut{$SNPname}) {
					push(@SNPsCNT, $reverseSynonyms{$SNPname});
					$origins{$reverseSynonyms{$SNPname}} = 'reference';
				} else {
					# nothing
				}
			}
		}

		foreach $SNPname (uniqueArray(@SNPsCNT)) {
			my $ori = $origins{$SNPname};
			my @SNPname = split(/_/, $SNPname);
			$SNPname =~ s/_/ - /g;
			
			my $id;
			if ($ori eq 'called') {
				$id = $name2pos{$SNPname[0]}.'_'.$name2mutant{$SNPname[0]};
			} else {
				$id = $name2pos{$SNPname[0]}.'_'.$name2ancestral{$SNPname[0]};
			}
			print CNT $id."\t".$SNPname."\t".$ori."\n";
		}
			
		close NCNT;
		close CNT;

		# Vertical method. 
		push(@linesToPrint,"VERTICAL\n");
		push(@linesToPrint,"--------\n");

		my @verticals;
		foreach $leaf (sort keys %posLeaves) {
			$leaf = trimSpace($leaf);
			if (substr($leaf, -1) ne '*') {
				my $SNPs = $posLeaves{$leaf};
				if ($SNPs ne '') { # don't include positive leaves that don't have SNPs
					my @SNPs = split(/\_/, $SNPs);
					my @resultsSNPs;
					foreach $SNP (@SNPs) {
						push(@resultsSNPs, $definiteResults{$SNP});
					}

					if ( getMetric(@resultsSNPs) > $cutoff ) {
						my $name = layout($leaf);
						$name = altName($name);
						push(@linesToPrint, "\$ $name\n");
						push(@verticals, $leaf);
					}
				}
			}
		}

		foreach $leaf (sort keys %negLeaves) {
			my $children = $negLeaves{$leaf};
			if ($children ne '') { # when some children have SNPs, only include leave when it is also positive && SNPs of children are too
				my @children = split(/\*/, $children);
				my @resChildren;
				foreach $child (@children) {
					my @SNPs = split(/\_/, $child);
					my @resultsSNPs;
					foreach $SNP (@SNPs) {
						push(@resultsSNPs, $definiteResults{"$SNP"});
					}
					push(@resChildren, getMetric(@resultsSNPs)) ;
				}

				my $res = 1;
				foreach $resChild (@resChildren) {
					if ($resChild > $cutoff) {
						$res = 0;
					}
				}
				my $currentSNPs = $posLeaves{$leaf};
				if ($currentSNPs ne '') {
					my @currentSNPs = split(/\_/, $currentSNPs);
					my @resCurrent;
					foreach $SNP (@currentSNPs) {
						push(@resCurrent, $definiteResults{$SNP});
					}
					if (getMetric(@resCurrent) < $cutoff) {
						$res = 0;
					}
				}

				if ($res == 1) {
					my $name = layout($leaf);
					$name = altName($name);
					push(@linesToPrint, "\$ $name\n");
					push(@verticals, $leaf);
				}
			} else {  # when all children don't have SNPs, only include leaf when it is also positive
				my $res = 1;
			
				my $currentSNPs = $posLeaves{$leaf};
				if ($currentSNPs ne '') {
					my @currentSNPs = split(/\_/, $currentSNPs);
					my @resCurrent;
					foreach $SNP (@currentSNPs) {
						push(@resCurrent, $definiteResults{$SNP});
					}
					if (getMetric(@resCurrent) < $cutoff) {
						$res = 0;
					}
				}

				if ($res == 1) {
					my $name = layout($leaf);
					$name = altName($name);
					push(@linesToPrint, "\$ $name\n");
					push(@verticals, $leaf);
				}	
			}			
		}

		# Horizontal method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "HORIZONTAL\n");
		push(@linesToPrint, "----------\n");
		my @toInvest = 'Root';
		my @horizontals;

		my %nodeScores;
		while(@toInvest) {
			my @temp;
			foreach $inv (@toInvest) {
				my @tempInv;
				my $children = $parent2nodes{$inv};
				foreach $child (split(/\_/, $children)) {
					my $SNPs = $node2snp{$child};
					if ($SNPs ne '') {
						my @SNPs = split(/\_/, $SNPs);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							if (substr($SNP, -1) eq "\r") { $SNP = substr($SNP, 0,-1) }
							push(@resultsSNPs, $definiteResults{$SNP});
						}
						if ( getMetric(@resultsSNPs) > $cutoff) {
							push(@tempInv, $child);
						}
						$nodeScores{$child} = getMetric(@resultsSNPs);
					} else {
						push(@tempInv, $child);
					}	
				}
				if (@tempInv) {
					foreach $elem (@tempInv) {
						push(@temp, $elem);
					}
				} else {
					push(@horizontals, $inv);
				}
			}
			@toInvest = @temp;
		}

		foreach $HR (@horizontals) {
			my $parents;
			my $parent = $node2parent{$HR};
			while ($parent ne '-') {
				$parents = $parents.$parent." ";
				$parent = $node2parent{$parent};
			}
			my $line = "! ".$HR."\t".$parents."\n";
			push(@linesToPrint, $line);
		}

		# Combi method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "COMBI\n");
		push(@linesToPrint, "-----\n");
		my @combis;
		if ($#verticals == -1) {
			@combis = @horizontals;
		} else {
			foreach $vert (@verticals) {
				my @path;
				if (substr($vert, -1) ne '*') {
					my $result = 0;
					my $parent = $vert;
					push(@path, $parent);
					while ($parent ne '-' && $parent ne 'Root') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
					foreach $hori (@horizontals) {
						if (isInArray($hori, @path)) {
							$result = 1;
						}
					}
					if ($result == 1) {
						push(@combis, $vert);
					}
				} else {
					my $result = 0;
					$vertWithoutStar = substr($vert, 0, -1);
					my $parent = $vertWithoutStar;
					push(@path, $parent);
					while ($parent ne '-' && $parent ne 'Root') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
					foreach $hori (@horizontals) {
						if (isInArray($hori, @path)) {
							$result = 1;
						}
					}
					if ($result == 1) {
						push(@combis, $vert);
					}
				}
			}
		}
		foreach $combi (@combis) {
			my $name = layout($combi);
			$name = altName($name);
			push(@linesToPrint, "% $name\n");
		}

		# Most specific method.
		push(@linesToPrint, "\n");
		push(@linesToPrint, "SPECIFIC\n");
		push(@linesToPrint, "--------\n");
		my %specific;
		foreach $combi (@combis) {
			my @path;
			my $result = 0;
			if (substr($combi, -1) ne '*') {
				my $parent = $combi;
				while ($parent ne '-') {
					$parent = $node2parent{$parent};
					push(@path, $parent);
				}

			} else {
				$combiWithoutStar = substr($combi, 0, -1);
				my $parent = $combiWithoutStar;
				while ($parent ne '-') {
					$parent = $node2parent{$parent};
					push(@path, $parent);
				}
			}
			my $hits = 0;
			foreach $comb (sort @combis) {

				if (substr($comb, -1) eq '*') {
					$combWithoutStar = substr($comb, 0,-1) ;
				} else {
					$combWithoutStar = $comb; 
				}
				if (isInArray($combWithoutStar, @path)) {
					$hits++;
				}
			}
			if ($hits > 0) {
				$specific{$combi} = $hits;
			}
		}
		my $prev;
		my @best;
		my $best = 1;
		foreach $spec (sort { $specific{$b} cmp $specific{$a} } keys %specific) {
			if ($best == 1) {
				push(@best, $spec);
				$best = $specific{$spec};
			} else {
				if ($best eq $specific{$spec}) {
					push(@best, $spec);
				}
			}
			$prev = $specific{$spec};
			my $name = layout($spec);
			$name = altName($name);
			my $line = "@ ".$name."\n"; 
			push(@linesToPrint, $line);
		}

		# Results. 
		push(@linesToPrint, "\n");
		push(@linesToPrint, "RESULTS\n");
		push(@linesToPrint, "-------\n");
		my @finalresults;
		if ($#best == -1) {
			foreach $res (@combis) {
				my $name = layout($res);
				$name = altName($name);
				$result = $name." ".$verticals{$res};
				push(@finalresults, $result);
			}

			if ($#combis == 0) {
				my $name = layout($combis[0]);
				$name = altName($name);

				push(@linesToPrint, "> $name\n");
			} else {
				my $bestGroup;
				my $bestGroupVerticalScore = 0;
							
				foreach $value (@combis) {
					my $name = layout($value);
					$name = altName($name);
					push(@linesToPrint, "> $name\n");
				}
			}
		} else {
			foreach $res (@best) {
				my $name = layout($res);
				$name = altName($name);
				$result = $name." ".$verticals{$res};
				push(@finalresults, $result);
			}
			if ($#best == 0) {
				my $name = layout($best[0]);
				$name = altName($name);
				push(@linesToPrint, "> $name\n");
			} else {
				my $bestGroup;
				my $bestGroupVerticalScore = 0;
				
				my @results;			
				foreach $res (@best) {
					push(@results, $res);
				}
				if ($#results == -1) {
					@results = @best;
				}
				
				foreach $res (@results) {
					my $name = layout($res);
					$name = altName($name);

					push(@linesToPrint, "> $name\n");
				}
			}
		}

		# Calculate quality.
		my $mcc = 100;
		my $runagain = 0;
		if ($#finalresults > 0  && $QCtype eq "sufficient") {
			$runagain = 1;
		} else {
			my $haplogroup;
			if ($#finalresults == 0) {
				$haplogroup = $finalresults[0];
			}

			if (substr($haplogroup, 0,3) eq 'R1b') {
				# get SNPs on its path
				$haplogroup  =~ s/\[.*\]//;
				$haplogroup = substr($haplogroup, 0, -1);
				my $name = $haplogroup;
				my %SNPsInPath;
				my %SNPsToIgnore;

				if ( $name =~ m/\*/ ) {
					while($name =~ m/\*/) {
						$name =~ s/\*//;
					}
					if (substr($name, -1) eq ' ') {
						$name = substr($name, 0, -1);
					}
					my $parent = $node2parent{$name};
					my @children = split(/\_/, $parent2nodes{$name});

					# To ignore
					my @childrenToIgnore;
					while (@children) {
						my $child = $children[0];
						if ($child ne $name) {
							push(@childrenToIgnore, $child);
							my @grandchildren = split(/\_/, $parent2nodes{$child}); 
							foreach $grandchild (@grandchildren) {
								push(@children, $grandchild);
							}
						}
						shift(@children);
					}
					@children = split(/\_/, $parent2nodes{$parent});
					# To ignore
					while (@children) {
						my $child = $children[0];
						if ($child ne $name) {
							push(@childrenToIgnore, $child);
							my @grandchildren = split(/\_/, $parent2nodes{$child}); 
							foreach $grandchild (@grandchildren) {
								push(@children, $grandchild);
							}
						}
						shift(@children);
					}

					foreach $childToIgnore (@childrenToIgnore) {
						my @SNPs = split(/\_/, $node2snp{$childToIgnore}); 
						foreach $SNP (@SNPs) {
							$SNPsToIgnore{$SNP} = 1;
						}
					}

					# In Path
					my @SNPs = split(/\_/, $node2snp{$name});
					foreach $SNP (@SNPs) {
						$SNPsInPath{$SNP} = 1;
					}

					while($parent ne 'Root') {
						my @SNPs = split(/\_/, $node2snp{$parent});
						foreach $SNP (@SNPs) {
							$SNPsInPath{$SNP} = 1;
						}
						$parent = $node2parent{$parent};
					}
				} else {
					# In Path
					while(substr($name, -1) eq ' ') {
						$name = substr($name, 0, -1);
					}
					my $parent = $node2parent{$name};
					my @SNPs = split(/\_/, $node2snp{$name});
					foreach $SNP (@SNPs) {
						$SNPsInPath{$SNP} = 1;
					}
					while($parent ne 'Root') {
						my @SNPs = split(/\_/, $node2snp{$parent});
						foreach $SNP (@SNPs) {
							$SNPsInPath{$SNP} = 1;
						}
						$parent = $node2parent{$parent};
					}
				}

				# Create hash with expected called SNPs. 
				my %CallExp;
				my %NotCallExp;
				foreach $SNP (@SNPsInTree) {
					if (!($SNPsToIgnore{$SNP})) {
						my $statRef = $referenceSNPs{$SNP};
						if ($SNPsInPath{$SNP}) {
							if ($statRef != 1) {
								$CallExp{$SNP} = 1;
							} else {
								$NotCallExp{$SNP} = 1;
							}
						} else {
							if ($statRef != 0) {
								$CallExp{$SNP} = 1;
							} else {
								$NotCallExp{$SNP} = 1;
							}
						}
					}
				}

				# Get qualities. 
				my $CallInExp = 0; 
				my $CallInNotExp = 0;
				my $Correct = 0;
				my $NotCorrect = 0;

				foreach $SNP (@SNPsInTree) {
					if (!($SNPsToIgnore{$SNP})) {
						my $status = $results{$SNP};
						if ($SNPsInPath{$SNP}) {
							if ($status == 1) {
								$Correct++;
							} else {
								$NotCorrect++;
							} 
						} else {
							if ($status == 0 || $status == -1) {
								$Correct++;
							} else {
								$NotCorrect++;
							} 
						}
						my $origin = $origins{$SNP};
						if ($CallExp{$SNP} && $origin eq 'called') {
							$CallInExp++;
						} 

						if ($NotCallExp{$SNP} && $origin eq 'called') {
							$CallInNotExp++;
						} 
					}
				}

				my @CallExp;
				foreach $key (sort keys %CallExp) {
					push(@CallExp, $key);
				}	
				my $CallExp = $#CallExp + 1;
						
				my @NotCallExp;
				foreach $key (keys %NotCallExp) {
					push(@NotCallExp, $key);
				}
				my $NotCallExp = $#NotCallExp + 1;

				my $TP = $CallInExp;
				my $FP = $CallInNotExp;
				my $TN = $NotCallExp - $CallInNotExp;
				my $FN = $CallExp - $CallInExp;

				$mcc = (($TP*$TN)-($FP*$FN))/sqrt(($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN))*100;
				push(@linesToPrint,"\n");
			}
			if ($mcc < 95 && $QCtype eq "sufficient") {
				$runagain = 1;
				push(@linesToPrint,"MCC\n");
				push(@linesToPrint,"---\n");
				push(@linesToPrint,"$mcc\n");
			}
		}

		if ($runagain == 1) {
		# Run again as insufficient.
			foreach $lineToPrint (@linesToPrint) {
				print OUT "### $lineToPrint";
			}

			print OUT "\n\nRESULTS FOR INSUFFICIENT CALL QUALITY BASED ON MCC\n";
			print OUT "--------------------------------------------------\n\n";
			foreach $line (@linesToPrint[0..2]) {
				$line =~ s/^ //;
				$line =~ s/ $//;
				print OUT $line;
			}
			print OUT "0 \(insufficient\)\n\n";

			$QCtype = 'insufficient';
			$cutoff = $cutoffLowQ;
			foreach $SNP (@SNPsInConv) {
				if ($SNPsInFile{$name2anc{$SNP}}) {
					$definiteResults{$SNP} = 0;
				} elsif ($SNPsInFile{$name2mut{$SNP}}) {
					$definiteResults{$SNP} = 1;
				} elsif ($reference{$SNP} eq $name2anc{$SNP}) {
					$definiteResults{$SNP} = 0;
				} elsif ($reference{$SNP} eq $name2mut{$SNP}) {
					$definiteResults{$SNP} = 0;
				} else {
					$definiteResults{$SNP} = -1;
				}
			}

			# Vertical method. 
			print OUT "VERTICAL\n--------\n"; $| = 1;
			my @verticals;
			foreach $leaf (sort keys %posLeaves) {
				$leaf = trimSpace($leaf);
				if (substr($leaf, -1) ne '*') {
					my $SNPs = $posLeaves{$leaf};
					if ($SNPs ne '') { # don't include positive leaves that don't have SNPs
						my @SNPs = split(/\_/, $SNPs);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							push(@resultsSNPs, $definiteResults{$SNP});
						}

						if ( getMetric(@resultsSNPs) > $cutoff ) {
							my $name = layout($leaf);
							$name = altName($name);
							print OUT "\$ $name\n"; $| = 1;
							push(@verticals, $leaf);
						}
					}
				}
			}

			foreach $leaf (sort keys %negLeaves) {
				my $children = $negLeaves{$leaf};
				if ($children ne '') { # when some children have SNPs, only include leave when it is also positive && SNPs of children are too
					my @children = split(/\*/, $children);
					my @resChildren;
					foreach $child (@children) {
						my @SNPs = split(/\_/, $child);
						my @resultsSNPs;
						foreach $SNP (@SNPs) {
							push(@resultsSNPs, $definiteResults{"$SNP"});
						}
						push(@resChildren, getMetric(@resultsSNPs)) ;
					}

					my $res = 1;
					foreach $resChild (@resChildren) {
						if ($resChild > $cutoff) {
							$res = 0;
						}
					}
				
					my $currentSNPs = $posLeaves{$leaf};
					if ($currentSNPs ne '') {
						my @currentSNPs = split(/\_/, $currentSNPs);
						my @resCurrent;
						foreach $SNP (@currentSNPs) {
							push(@resCurrent, $definiteResults{$SNP});
						}
						if (getMetric(@resCurrent) < $cutoff) {
							$res = 0;
						}
					}

					if ($res == 1) {
						my $name = layout($leaf);
						$name = altName($name);
						print OUT "\$ $name\n"; $| = 1;
						push(@verticals, $leaf);
					}
				} else {  # when all children don't have SNPs, only include leaf when it is also positive
					my $res = 1;
				
					my $currentSNPs = $posLeaves{$leaf};
					if ($currentSNPs ne '') {
						my @currentSNPs = split(/\_/, $currentSNPs);
						my @resCurrent;
						foreach $SNP (@currentSNPs) {
							push(@resCurrent, $definiteResults{$SNP});
						}
						if (getMetric(@resCurrent) < $cutoff) {
							$res = 0;
						}
					}

					if ($res == 1) {
						my $name = layout($leaf);
						$name = altName($name);
						print OUT "\$ $name\n"; $| = 1;
						push(@verticals, $leaf);
					}	
				}			
			}


			# Horizontal method.
			print OUT "\nHORIZONTAL\n----------\n"; $| = 1;
			my @toInvest = 'Root';
			my @horizontals;

			my %nodeScores;
			while(@toInvest) {
				my @temp;
				foreach $inv (@toInvest) {
					my @tempInv;
					my $children = $parent2nodes{$inv};
					foreach $child (split(/\_/, $children)) {
						my $SNPs = $node2snp{$child};
						if ($SNPs ne '') {
							my @SNPs = split(/\_/, $SNPs);
							my @resultsSNPs;
							foreach $SNP (@SNPs) {
								if (substr($SNP, -1) eq "\r") { $SNP = substr($SNP, 0,-1) }
								push(@resultsSNPs, $definiteResults{$SNP});
							}
							if ( getMetric(@resultsSNPs) > $cutoff) {
								push(@tempInv, $child);
							}
							$nodeScores{$child} = getMetric(@resultsSNPs);
						} else {
							push(@tempInv, $child);
						}	
					}
					if (@tempInv) {
						foreach $elem (@tempInv) {
							push(@temp, $elem);
						}
					} else {
						push(@horizontals, $inv);
					}
				}
				@toInvest = @temp;
			}


			foreach $HR (@horizontals) {
				my $parents;
				my $parent = $node2parent{$HR};
				while ($parent ne '-') {
					$parents = $parents.$parent." ";
					$parent = $node2parent{$parent};
				}
				print OUT "! ".$HR."\t".$parents."\n"; $| = 1;
			}

			# Combi method.
			print OUT "\nCOMBI\n-----\n"; $| = 1;
			my @combis;
			if ($#verticals == -1) {
				@combis = @horizontals;
			} else {
				foreach $vert (@verticals) {
					my @path;
					if (substr($vert, -1) ne '*') {
						my $result = 0;
						my $parent = $vert;
						push(@path, $parent);
						while ($parent ne '-' && $parent ne 'Root') {
							$parent = $node2parent{$parent};
							push(@path, $parent);
						}
						foreach $hori (@horizontals) {
							if (isInArray($hori, @path)) {
								$result = 1;
							}
						}
						if ($result == 1) {
							push(@combis, $vert);
						}
					} else {
						my $result = 0;
						$vertWithoutStar = substr($vert, 0, -1);
						my $parent = $vertWithoutStar;
						push(@path, $parent);
						while ($parent ne '-') {
							$parent = $node2parent{$parent};
							push(@path, $parent);
						}
						foreach $hori (@horizontals) {
							if (isInArray($hori, @path)) {
								$result = 1;
							}
						}
						if ($result == 1) {
							push(@combis, $vert);
						}
					}
				}
			}
			foreach $combi (@combis) {
				my $name = layout($combi);
				$name = altName($name);
				print OUT "% ".$name."\n"; $| = 1;
			}

			# Most specific method.
			print OUT "\nSPECIFIC\n--------\n"; $| = 1;
			my %specific;
			foreach $combi (@combis) {
				my @path;
				my $result = 0;
				if (substr($combi, -1) ne '*') {
					my $parent = $combi;
					while ($parent ne '-') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}

				} else {
					$combiWithoutStar = substr($combi, 0, -1);
					my $parent = $combiWithoutStar;
					while ($parent ne '-') {
						$parent = $node2parent{$parent};
						push(@path, $parent);
					}
				}
				my $hits = 0;
				foreach $comb (sort @combis) {

					if (substr($comb, -1) eq '*') {
						$combWithoutStar = substr($comb, 0,-1) ;
					} else {
						$combWithoutStar = $comb; 
					}
					if (isInArray($combWithoutStar, @path)) {
						$hits++;
					}
				}
				if ($hits > 0) {
					$specific{$combi} = $hits;
				}
			}
			my $prev;
			my @best;
			my $best = 1;
			foreach $spec (sort { $specific{$b} cmp $specific{$a} } keys %specific) {
				if ($best == 1) {
					push(@best, $spec);
					$best = $specific{$spec};
				} else {
					if ($best eq $specific{$spec}) {
						push(@best, $spec);
					}
				}
				$prev = $specific{$spec};
				my $name = layout($spec);
				$name = altName($name);
				print OUT "@ ".$name."\n"; $| = 1;
			}

			# Results. 
			print OUT "\nRESULTS\n-------\n"; $| = 1;
			if ($#best == -1) {
				if ($#combis == 0) {
					my $name = layout($combis[0]);
					$name = altName($name);
					print OUT "> ".$name."\n"; $| = 1;
				} else {
					my $bestGroup;
					my $bestGroupVerticalScore = 0;
								
					foreach $value (@combis) {
						my $name = layout($value);
						$name = altName($name);
						print OUT "> ".$name." ".$verticals{$res}."\n"; $| = 1;
					}
				}
			} else {
				if ($#best == 0) {
					my $name = layout($best[0]);
					$name = altName($name);
					print OUT "> ".$name."\n"; $| = 1;
				} else {
					my $bestGroup;
					my $bestGroupVerticalScore = 0;
					
					my @results;			
					foreach $res (@best) {
						push(@results, $res);
					}
					if ($#results == -1) {
						@results = @best;
					}
					
					foreach $res (@results) {
						my $name = layout($res);
						$name = altName($name);
						print OUT "> ".$name." ".$verticals{$res}."\n"; $| = 1;
					}
				}
			}
		} else {
			foreach $lineToPrint (@linesToPrint) {
				print OUT "$lineToPrint";
			}
		}
		close OUT;
		
		SHELLY($analysisFile,$statusFile,$treeFile,$refFile,$outFile,$statusOutFile);
	}
} else {
	print "This file does not contain any SNPs.\n";
}



