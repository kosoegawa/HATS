#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: CBw.pm 
# This module was developed to capture HLA-C Bw4 and Bw6 alleles
# last reviewed on November 14 2023

package CBw;
use strict;

#my @bw = (76,77,80,82,83);
my @bw = (2,76,80,82,83);
my %bw;
my %ref;
$ref{"C-0102"} = "HLA00401";	# C*01:02:01:01
$ref{"C-0202"} = "HLA00404";	# C*02:02:01
#$ref{"C-0310"} = "HLA:HLA01076";	# C*03:10
$ref{"C-1212"} = "HLA01878";	# C*12:12, C1 Bw6
$ref{"C-0669"} = "HLA07131";	# C*06:69, C2 Bw6, this is required, do not delete
$bw{"C-0102"} = "Negative"; $bw{"C-0202"} = "Negative";
#$bw{"C-0310"} = "C2";
$bw{"C-1212"} = "Bw6"; $bw{"C-0669"} = "Bw6";

sub HLAC {
	my $gene = "C";
	return $gene;
}

sub HLAC_LEADER {
	my $leader = 23;		# C specific
	return $leader;
}

sub RESIDUES {
	my @residues = @bw;
	my $residues_ref = \@residues;
	return $residues_ref;
}

sub REF {
	my $ref_ref = \%ref;
	return $ref_ref;
}

sub BW {
	my $bw_ref = \%bw;
	return $bw_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 25;	#change the number of missing nucleotide
	$partial{ "C-0669" } = $seq;
		
	return $partial_ref;
}


1;
