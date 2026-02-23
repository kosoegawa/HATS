#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: CC1C2.pm 
# This module was developed to capture HLA-C Bw4 and Bw6 alleles
# last reviewed on February 19 2026

package CC1C2;
use strict;

my @c1c2 = (2,76,77,80);
#my @c1c2 = (2,76,80,82);
my %c1c2;
my %ref;
$ref{"Cw0102"} = "HLA00401";	# C*01:02:01:01
$ref{"Cw0202"} = "HLA00404";	# C*02:02:01
$ref{"Cw0310"} = "HLA01076";	# C*03:10
$ref{"Cw1404"} = "HLA00465";	# C*14:04
$ref{"B7301"} = "HLA00392";	# B*73:01:01:01
#$ref{"Cw1212"} = "HLA01878";	# C*12:12, C1 Bw6
#$ref{"Cw0669"} = "HLA07131";	# C*06:69, C2 Bw6, this is required, do not delete
$c1c2{"Cw0102"} = "C1"; $c1c2{"Cw0202"} = "C2";
$c1c2{"Cw0310"} = "Uncertain";
$c1c2{"Cw1404"} = "Uncertain";
$c1c2{"B7301"} = "Uncertain";
#$c1c2{"C-1212"} = "Bw6"; $c1c2{"C-0669"} = "Bw6";

sub HLAC {
	my $gene = "C";
	return $gene;
}

sub HLAC_LEADER {
	my $leader = 23;		# C specific
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
	my $seq = "N" x 25;	#change the number of missing nucleotide
		
	return $partial_ref;
}


1;
