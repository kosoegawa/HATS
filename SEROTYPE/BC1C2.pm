#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: BC1C2.pm 
# This module was developed to capture Bw4 and Bw6 alleles
# last modified and documented on January 27 2022

package BC1C2;
use strict;

#my @c1c2 = (43,52,76,79,82,83);		#used 79 instead of 80
#my @c1c2 = (2,76,79,82,83);		#used 79 instead of 80
my @c1c2 = (2,76,77,80,82);		#used 79 instead of 80
my %ref;
my %c1c2;
#$ref{"B-4601"} = "HLA00331";	#B*46:01:01:01, Negative
#$c1c2{"B-4601"} = "C1";
$ref{"C-0102"} = "HLA00401";	# C*01:02:01:01
$ref{"C-0202"} = "HLA00404";	# C*02:02:01
$ref{"C-0310"} = "HLA:HLA01076";	# C*03:10
$ref{"C-1404"} = "HLA00465";	# C*14:04
$ref{"B-7301"} = "HLA00392";	# B*73:01:01:01
$c1c2{"C-0102"} = "C1"; $c1c2{"C-0202"} = "C2";
$c1c2{"C-0310"} = "Uncertain";
$c1c2{"C-1404"} = "Uncertain";
$c1c2{"B-7301"} = "C1";

sub HLAB_LEADER {
	my $leader = 23;		# A specific
	return $leader;
}

sub RESIDUES {
	my @residues = @c1c2;
	my $residues_ref = \@residues;
	return $residues_ref;
}

sub REF {
	my $ref_ref = \%ref;
	return $ref_ref;
}

sub C1C2 {
	my $c1c2_ref = \%c1c2;
	return $c1c2_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
		
	return $partial_ref;
}


1;
