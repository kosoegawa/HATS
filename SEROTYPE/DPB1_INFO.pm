#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DPB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last reviewed, modified and documented on March 7 2026

package DPB1_INFO;
use strict;

my @dpb01 = (56,57,69,84,85,86);	# required for the alleles missing exon 3
my @dpb15 = (56,57,69,84,85,86,96);

my %group;
my %base;
my %dpb01;
my %dpb15;
$dpb01{"DPB01"} = "HLA00514";	# DP-B1*01:01:01:01 
$dpb01{"DPB0201"} = "HLA00517";	# DP-B1*02:01:02:01
$dpb01{"DPB0202"} = "HLA00519";	# DP-B1*02:02:01:01
$dpb01{"DPB03"} = "HLA00520";	# DP-B1*03:01:01:01
$dpb01{"DPB0401"} = "HLA00521";	# DP-B1*04:01:01:01
$dpb01{"DPB0402"} = "HLA00522";	# DP-B1*04:02:01:01
$dpb01{"DPB06"} = "HLA00524";	# DP-B1*06:01:01:01
$dpb01{"DPB10"} = "HLA00527";	# DP-B1*10:01:01:01
$dpb01{"DPB13"} = "HLA00530";	# DP-B1*13:01:01:01
$dpb01{"DPB45"} = "HLA00564";	# DP-B1*45:01
$dpb01{"DPB46"} = "HLA00565";	# DP-B1*46:01:01
$dpb01{"DPB80"} = "HLA00599";	# DP-B1*80:01

$dpb15{"DPB15"} = "HLA00532";	# DP-B1*15:01:01:01
$dpb15{"DPB17"} = "HLA00534";	# DPB1*17:01:01:01 258 bp
$dpb15{"DPB18"} = "HLA00535";	# DP-B1*18:01:01:01
$dpb15{"DPB30"} = "HLA00549";	# DPB1*30:01:01:01 258 bp
$dpb15{"DPB31"} = "HLA00550";	# DPB1*31:01:01:01 258 bp
 
$group{"DPB01"} = "DPB01"; $group{"DPB0201"} = "DPB01"; $group{"DPB0202"} = "DPB01"; $group{"DPB03"} = "DPB01"; $group{"DPB0401"} = "DPB01"; $group{"DPB0402"} = "DPB01";
$group{"DPB06"} = "DPB01"; $group{"DPB10"} = "DPB01"; $group{"DPB13"} = "DPB01"; $group{"DPB45"} = "DPB01"; $group{"DPB46"} = "DPB01"; $group{"DPB80"} = "DPB01";

$group{"DPB15"} = "DPB15"; $group{"DPB17"} = "DPB15"; $group{"DPB18"} = "DPB15"; $group{"DPB30"} = "DPB15"; $group{"DPB31"} = "DPB15";

$base{"DPB01"} = "DPB01"; $base{"DPB0201"} = "DPB0201"; $base{"DPB0202"} = "DPB0202"; $base{"DPB03"} = "DPB03"; $base{"DPB0401"} = "DPB0401"; $base{"DPB0402"} = "DPB0402";
$base{"DPB06"} = "DPB06"; $base{"DPB10"} = "DPB10"; $base{"DPB13"} = "DPB13"; $base{"DPB45"} = "DPB45"; $base{"DPB46"} = "DPB46"; $base{"DPB80"} = "DPB80";
$base{"DPB15"} = "DPB15"; $base{"DPB17"} = "DPB17"; $base{"DPB18"} = "DPB18"; $base{"DPB30"} = "DPB30"; $base{"DPB31"} = "DPB31";

my @subtype = ();	# modify here if serotype modified

sub DPB1 {
	my $gene = "DPB1";
	return $gene;
}

sub DPB1_LEADER {
	my $leader = 28;		# DPB1 specific
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
	my @basetype = ("DPB01", "DPB0201", "DPB0202", "DPB03", "DPB0401", "DPB0402", "DPB06",
		"DPB10", "DPB13", "DPB15", "DPB17", "DPB18", "DPB30", "DPB31", "DPB45", "DPB46", "DPB80");	#
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
	push @combined, @dpb01;
	push @combined, @dpb15;
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
	@residues = @dpb01;
		if ( $serotype eq "DPB01" ) {
		@residues = @dpb01;
	}
	elsif ( $serotype eq "DPB15" ) {
		@residues = @dpb15;
	}
	else {
		@residues = @unique;
	}
	return $residues_ref;
}

sub REF {
	my ( $serotype ) = @_;	# unnecessary for DPB1
	my %ref;
	my $ref_ref = \%ref;
	
	if ( $serotype eq "DPB01" ) {
		%ref = %dpb01;
	}
	elsif ( $serotype eq "DPB15" ) {
		%ref = %dpb15;
	}
	else {	# all together
		%ref = (%dpb01,%dpb15);
	}

	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dpb01,%dpb15);
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
	my %tmp = (%dpb01,%dpb15);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key eq "DPB01" ) {
			$ref{$key} = "DPB1\\*01";
		}
		elsif (( $key eq "DPB0201" ) || ( $key eq "DPB0202" )) {
			$ref{$key} = "DPB1\\*02";
		}
		elsif ( $key eq "DPB03" ) {
			$ref{$key} = "DPB1\\*03";
		}
		elsif (( $key eq "DPB0401" ) || ( $key eq "DPB0402")) {
			$ref{$key} = "DPB1\\*04";
		}
		elsif ( $key eq "DPB06" ) {
			$ref{$key} = "DPB1\\*06";
		}
		elsif ( $key eq "DPB10" ) {
			$ref{$key} = "DPB1\\*10:";
		}
		elsif ( $key eq "DPB13" ) {
			$ref{$key} = "DPB1\\*13:";
		}
		elsif ( $key eq "DPB15" ) {
			$ref{$key} = "DPB1\\*15:";
		}
		elsif ( $key eq "DPB17" ) {
			$ref{$key} = "DPB1\\*17:";
		}
		elsif ( $key eq "DPB18" ) {
			$ref{$key} = "DPB1\\*18:";
		}
		elsif ( $key eq "DPB30" ) {
			$ref{$key} = "DPB1\\*30:";
		}
		elsif ( $key eq "DPB31" ) {
			$ref{$key} = "DPB1\\*31:";
		}
		elsif ( $key eq "DPB45" ) {
			$ref{$key} = "DPB1\\*45:";
		}
		elsif ( $key eq "DPB46" ) {
			$ref{$key} = "DPB1\\*46:";
		}
		elsif ( $key eq "DPB80" ) {
			$ref{$key} = "DPB1\\*80:";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 34;	#change the number of missing nucleotide
	$partial{ "general" } = $seq;
	my $seq2 = "N" x 36;	#change the number of missing nucleotide
	$partial{ "DPB1*37:01" } = $seq2;
	$partial{ "DPB1*56:01" } = $seq2;
	$partial{ "DPB1*74:01" } = $seq2;
	$partial{ "DPB1*75:01" } = $seq2;
	$partial{ "tail" } = "X" x 100;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "DPB01" } = 0;
	$known_cross{ "DPB0201" } = 0;
	$known_cross{ "DPB0202" } = 0;
	$known_cross{ "DPB03" } = 0;
	$known_cross{ "DPB0401" } = 0;
	$known_cross{ "DPB0402" } = 0;
	$known_cross{ "DPB06" } = 0;
	$known_cross{ "DPB10" } = 0;
	$known_cross{ "DPB13" } = 0;
	$known_cross{ "DPB15" } = 0;
	$known_cross{ "DPB17" } = 0;
	$known_cross{ "DPB18" } = 0;
	$known_cross{ "DPB30" } = 0;
	$known_cross{ "DPB31" } = 0;
	$known_cross{ "DPB45" } = 0;
	$known_cross{ "DPB46" } = 0;
	$known_cross{ "DPB80" } = 0;
	return $known_cross_ref;
}

1;
