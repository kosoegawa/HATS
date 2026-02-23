#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DRB3_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on February 1 2026

package DRB3_INFO;
use strict;

my @dr52 = (9, 10, 12, 13);
my @extra = (11,16,47,58,60,67,70,71,74);	# FULL only
#added residue 60

my %dr52;
my %group;
my %base;
$dr52{"DR5201"} = "HLA00887";		# DRB3*01:01:02:01
$dr52{"DR5202"} = "HLA00895";		# DRB3*02:02:01:01
$dr52{"DR5203"} = "HLA00902";		# DRB3*03:01:01:01
$group{"DR5201"} = "DR52";
$group{"DR5202"} = "DR52";
$group{"DR5203"} = "DR52";
$base{"DR5201"} = "DR52";
$base{"DR5202"} = "DR52";
$base{"DR5203"} = "DR52";


my @subtype = ("DR5202","DR5203");	# modify here if serotype modified

sub DRB3 {
	my $gene = "DRB3";
	return $gene;
}

sub DRB3_LEADER {
	my $leader = 28;		# DRB3 specific
	return $leader;
}

sub GROUP {
	my $group_ref = \%group;
	return $group_ref;
}

sub BASE {
	my $base_ref = \%base;
	return $base_ref;
}

sub BASETYPE {
	my @basetype = ("DR52");

	my $basetype_ref = \@basetype;
}

sub PARENT {
	my %parent;
	my $parent_ref = \%parent;
	foreach my $key ( keys %base ) {	# $key = A0101
		$parent{ $key } = $base{ $key };	# $base( $key } = "A1";
	}
	return $parent_ref;
}

sub BROAD {
	my %broad;
	foreach my $base ( keys %base ) {
		$broad{ $base } = $base{ $base };
	}
	my $broad_ref = \%broad;
	return $broad_ref;
}

sub RESIDUES {
	my ( $serotype ) = @_;
	my @combined = ();
	push @combined, @dr52;
	push @combined, @extra;
	my %seen;
	my @unique;
	foreach my $value ( sort { $a <=> $b } @combined ) {
		unless ( exists $seen{ $value } ) {
			push @unique, $value;
			$seen{ $value } = 0;
		}
	}
	
	my @residues = ();
	my $residues_ref = \@residues;
	if ( $serotype =~ /DR52/ ) {
		@residues = @dr52;
	}
	else {
		@residues = @unique;
	}
	return $residues_ref;
}

sub REF {
	my ( $serotype ) = @_;
	my %ref;
	my $ref_ref = \%ref;

	if ( $serotype =~ /DR52/ ) {
		%ref = %dr52;
	}
	else {
		%ref = (%dr52);
	}
	
	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dr52);
	my @tmp = sort keys %ref;
	for ( my $index = 0; $index < scalar @tmp; $index++ ) {
		$sero[0][$index] = $tmp[$index];
	}
	for ( my $index = 0; $index < scalar @subtype; $index++ ) {
		$sero[1][$index] = $subtype[$index];
	}
	my $sero_ref = \@sero;
	return $sero_ref;
}

sub KEY {
	my %tmp = (%dr52);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /DR52/ ) {
			$ref{$key} = "DRB3\\*";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "X" x 34;
	$partial{ "general" } = $seq;
	$partial{ "DRB3*01:02" } = "X" x 29;
	$partial{ "DRB3*02:11" } = "X" x 29;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "DR5202" } = 0;
	$known_cross{ "DR5203" } = 0;
	return $known_cross_ref;
}

1;
