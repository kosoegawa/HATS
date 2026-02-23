#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: BBw.pm 
# This module was developed to capture Bw4 and Bw6 alleles
# last reviewed on February 20 2026

package BBw;
use strict;

my @bw = (2,76,79,82,83);		#used 79 instead of 80, added 77
my %ref;
my %bw;
$ref{"B0702"} = "HLA00132";		#B*07:02:01:01, Bw6, removed dash (-)
$ref{"B1302"} = "HLA00153";	#B*13:02:01:01, Bw4
$ref{"B4601"} = "HLA00331";	#B*46:01:01:01, Negative
$bw{"B1302"} = "Bw4"; $bw{"B0702"} = "Bw6"; $bw{"B4601"} = "Negative";
#$bw{"B54"} = "Bw6";

sub HLAB_LEADER {
	my $leader = 23;		# A specific
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
#	my $seq = "N" x 25;	#change the number of missing nucleotide
#	$partial{ "Cw669" } = $seq;
		
	return $partial_ref;
}


1;
