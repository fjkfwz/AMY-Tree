#!/Perl/bin/perl 
use Bio::SeqIO;

#################################################################################
#										#
#	Script name:	getStatusReference.pl					#
#	Author: 	Anneleen Van Geystelen					#
#	Version:	v1.0							#
#	Created:	23/08/2011						#
#	Last Modified:	23/08/2011						#
#										#
#	Description:	This Perl script will create a status file of a givern	#
#			reference fasta file. 					#
#										#
#										#
#	INPUT:									#
#	    variable:								#
# 	    fixed: 	path of reference fasta file				#
#			path of conversion file					#
#			version (i.e. hg18 or hg19)				#
#										#
#	OUTPUT: 	status file						#
#										#
# 	PERL MODULES:								#
#										#
#################################################################################

#################################################################################
#				VARIABLES					#
#################################################################################

if ($#ARGV != 3) {
	print "Please, give valid arguments\n"; die;
}

# Set path of reference fasta file.
my $refFile = $ARGV[0];
if (!(-f $refFile)) { die "$refFile does not exist\n"; }

# Set path of the conversion file. 
my $convFile = $ARGV[1];
if (!(-f $convFile)) { die "$convFile does not exist\n"; }

# Set version variable. 
my $version = $ARGV[2];

# Set path of output status file. 
my $statFile = $ARGV[3];

#################################################################################
#				END OF VARIABLES				#
#################################################################################

#################################################################################
#				SUBROUTINES					#
#################################################################################

# Function to trim leading and trailing spaces.
sub trimSpace {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

#################################################################################
#				END OF SUBROUTINES				#
#################################################################################


######################## Create hash of conversion file. ########################

my %name2anc;
my %name2mut;
my %name2pos;

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
				$pos36 = trimSpace($pos36);
				$name2anc{$name} = $anc;
				$name2mut{$name} = $mut;
				$name2pos{$name} = $pos36;
			} elsif ($version eq 'hg19') {
				$pos37 = trimSpace($pos37);
				$name2anc{$name} = $anc;
				$name2mut{$name} = $mut;
				$name2pos{$name} = $pos37;
			}
		}
	} else {
		$flag = 1;
	}
} 
close CONV;


######################## Create hash of reference. #######################

# Get reference sequence.
my $refSeq;
my $in  = Bio::SeqIO->new(-file => "$refFile" , -format => 'Fasta');
while(my $sobj = $in->next_seq) {
	my $id = $sobj->id;
	$refSeq = $sobj->seq();
}
 
################### Create hash with status of each SNP. ##################

# Get reference bases of SNP.
my %reference;
foreach $name (sort keys %name2pos) {
	my $position = $name2pos{$name} - 1;
	my $mutant = $name2mut{$name};
	my $ancestral = $name2anc{$name};
	
	my $base = substr($refSeq, $position, 1);
	$base =~ tr/a-z/A-Z/;
	if ($base eq $mutant) {
		$reference{$name} = 1;
	} elsif ($base eq $ancestral) {
		$reference{$name} = 0;
	} else {
		$reference{$name} = -1
	}
}

######################### Write output. ####################################

open (OUT, "> $statFile");

foreach $name (sort keys %reference) {
	my $status = $reference{$name};
	print OUT $name."\t".$status."\n";
}

close OUT;
