#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: BC1C2.pm 
# This module was developed to capture Bw4 and Bw6 alleles
# last reviewed on February 20 2026

package BC1C2;
use strict;

my @c1c2 = (76,77,80,81,82,83);		#used 79 instead of 80
my %ref;
my %c1c2;
$ref{"Cw0102"} = "HLA00401";	# C*01:02:01:01
$ref{"Cw0202"} = "HLA00404";	# C*02:02:01
$ref{"Cw0310"} = "HLA01076";	# C*03:10
$ref{"Cw1404"} = "HLA00465";	# C*14:04
$ref{"B7301"} = "HLA00392";	# B*73:01:01:01
$c1c2{"Cw0102"} = "C1"; $c1c2{"Cw0202"} = "C2";
$c1c2{"Cw0310"} = "Uncertain";
$c1c2{"Cw1404"} = "Uncertain";
$c1c2{"B7301"} = "C1";

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
