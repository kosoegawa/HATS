#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: ABw.pm 
# This module was developed to capture HLA-A Bw4 and Bw6 alleles
# last modified and documented on April 20 2020

package ABw;
use strict;

#my @bw = (76,80,82,83);
my @bw = (2,76,82,83);		#used 79 instead of 80
#my @bw = (43,76,80,82,83);
my %ref;
my %bw;
$ref{"A-3201"} = "HLA00101";		# ref for Bw4: A*32:01:01:01
$ref{"A-3002"} = "HLA00090";		# ref for Bw6: A*30:02:01:01
$bw{"A-3201"} = "Bw4"; $bw{"A-3002"} = "Bw6";

sub HLAA {
	my $gene = "A";
	return $gene;
}

sub HLAA_LEADER {
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
