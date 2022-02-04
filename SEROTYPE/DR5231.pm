#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DR5231.pm 
# This module was developed to capture DRB1 alleles with DR52, DR53 and DR51 serotypes
# last modified and documented on May 15 2020

package DR5231;
use strict;

my @dr52 = (9, 10, 12, 13);
my @dr53 = (9, 10, 11, 12, 13);
my @dr51 = (9, 10, 11, 12, 13);


my %ref;
my %dr5231;
my %group;
my %base;
$ref{"DR52"} = "HLA00887";		# DRB3*01:01:02:01
$ref{"DR53"} = "HLA00905";		# DRB4*01:01:01:01
$ref{"DR51"} = "HLA00915";		# DRB5*01:01:01:01
$dr5231{"DR52"} = "DR52";
$dr5231{"DR53"} = "DR53";
$dr5231{"DR51"} = "DR51";

sub DRB1 {
	my $gene = "DRB1";
	return $gene;
}

sub DRB1_LEADER {
	my $leader = 28;		# DRB1 specific
	return $leader;
}

sub RESIDUES {
	my ( $serotype ) = @_;
	my @residues = ();
	my $residues_ref = \@residues;
	if ( $serotype eq "DR52" ) {
		@residues = @dr52; 
	}
	elsif ( $serotype eq "DR53" ) {
		@residues = @dr53; 
	}
	elsif ( $serotype eq "DR51" ) {
		@residues = @dr51; 
	}
	return $residues_ref;
}

sub REF {
	my $ref_ref = \%ref;
	return $ref_ref;
}

sub DR5231 {
	my $dr5231_ref = \%dr5231;
	return $dr5231_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 34;
#	$partial{ "DR1209" } = $seq;
		
	return $partial_ref;
}

1;
