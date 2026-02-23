#!/usr/bin/perl -w
#
# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DRB5_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last reviewed, modified and documented on February 1 2026

package DRB5_INFO;
use strict;

my @dr51 = (9, 10, 11, 12, 13);
my @extra = (16,47,58,60,67,70,71,74);	# FULL only

my %dr51;
my %group;
my %base;
$dr51{"DR5101"} = "HLA00915";		# DRB5*01:01:01:01
$dr51{"DR5102"} = "HLA00926";		# DRB5*02:02:01
$dr51{"DR5103"} = "HLA00918";		# DRB5*01:03
$group{"DR5101"} = "DR51";
$group{"DR5102"} = "DR51";
$group{"DR5103"} = "DR51";
$base{"DR5101"} = "DR51";
$base{"DR5102"} = "DR51";
$base{"DR5103"} = "DR51";


my @subtype = ("DR5102","DR5103");	# modify here if serotype modified

sub DRB5 {
	my $gene = "DRB5";
	return $gene;
}

sub DRB5_LEADER {
	my $leader = 28;		# DRB5 specific
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
	my @basetype = ("DR51");

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
	push @combined, @dr51;
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
	if ( $serotype =~ /DR51/ ) {
		@residues = @dr51;
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

	if ( $serotype =~ /DR51/ ) {
		%ref = %dr51;
	}
	else {
		%ref = (%dr51);
	}
	
	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dr51);
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
	my %tmp = (%dr51);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /DR51/ ) {
			$ref{$key} = "DRB5\\*";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 34;
	$partial{ "general" } = $seq;
	$partial{ "DRB5*01:04" } = "N" x 29;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "DR5102" } = 0;
	$known_cross{ "DR5103" } = 0;
	return $known_cross_ref;
}

1;
