#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DPB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on March 25 2024

package DPB1_INFO;
use strict;

#my @dpb1 = (56,57,69,82,84,85);	# 82 is included to be DPB1 specific
my @dpb1 = (56,57,69,84,85,86);	# 82 is included to be DPB1 specific

my %group;
my %base;
my %dpb1;
$dpb1{"DP-0101"} = "HLA00514";	# DP-B1*01:01:01:01 
$dpb1{"DP-0201"} = "HLA00517";	# DP-B1*02:01:02:01
$dpb1{"DP-0301"} = "HLA00520";	# DP-B1*03:01:01:01
$dpb1{"DP-0401"} = "HLA00521";	# DP-B1*04:01:01:01
$dpb1{"DP-1001"} = "HLA00527";	# DP-B1*10:01:01:01
$dpb1{"DP-1501"} = "HLA00532";	# DP-B1*15:01:01:01
$dpb1{"DP-1801"} = "HLA00535";	# DP-B1*18:01:01:01
$dpb1{"DP-4601"} = "HLA00565";	# DP-B1*46:01:01
$dpb1{"DP-0402"} = "HLA00522";	# DP-B1*04:02:01:01
$dpb1{"DP-0202"} = "HLA00519";	# DP-B1*02:02:01:01
$dpb1{"DP-0601"} = "HLA00524";	# DP-B1*06:01:01:01
$dpb1{"DP-1301"} = "HLA00530";	# DP-B1*13:01:01:01
$dpb1{"DP-4501"} = "HLA00564";	# DP-B1*45:01
$dpb1{"DP-8001"} = "HLA00599";	# DP-B1*80:01

$dpb1{"DP-1101"} = "HLA00528";	# DP-B1*11:01:01:01, 69R is equivalent to 69K, this is DP-1
$dpb1{"DP-3401"} = "HLA00553";	# DP-B1*34:01:01:01, 69K is equivalent to 69R, this is DP-15
$dpb1{"DP-13601"} = "HLA07366";	# DP-B1*136:01, 69R is equivalent to 69K, this is equivalent to DP-45
$dpb1{"DP-6901"} = "HLA00588";	# DP-B1*69:01:01:01

$group{"DP-0101"} = "DP-01"; $group{"DP-0201"} = "DP-01"; $group{"DP-0301"} = "DP-01"; $group{"DP-0401"} = "DP-01"; $group{"DP-1001"} = "DP-01";
$group{"DP-1501"} = "DP-01"; $group{"DP-1801"} = "DP-01"; $group{"DP-4601"} = "DP-01"; $group{"DP-0402"} = "DP-01"; $group{"DP-0202"} = "DP-01";
$group{"DP-0601"} = "DP-01"; $group{"DP-1301"} = "DP-01"; $group{"DP-1101"} = "DP-01"; $group{"DP-3401"} = "DP-01"; $group{"DP-4501"} = "DP-01";
$group{"DP-13601"} = "DP-01"; $group{"DP-8001"} = "DP-01"; $group{"DP-6901"} = "DP-01";
$base{"DP-0101"} = "DP-0101"; $base{"DP-0201"} = "DP-0201"; $base{"DP-0301"} = "DP-0301"; $base{"DP-0401"} = "DP-0401";
$base{"DP-1001"} = "DP-1001"; $base{"DP-1501"} = "DP-0401"; $base{"DP-1801"} = "DP-0402"; $base{"DP-4601"} = "DP-4601";	
$base{"DP-0402"} = "DP-0402"; $base{"DP-0202"} = "DP-0202"; $base{"DP-0601"} = "DP-0601"; $base{"DP-1301"} = "DP-1301";
$base{"DP-4501"} = "DP-4501"; $base{"DP-8001"} = "DP-8001";
$base{"DP-1101"} = "DP-0101"; $base{"DP-3401"} = "DP-0401"; $base{"DP-13601"} = "DP-4501";  $base{"DP-6901"} = "DP-0301";
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
	my @basetype = ("DP-0101", "DP-0201", "DP-0301", "DP-0401", "DP-1001", "DP-4601",
		"DP-0402", "DP-0202", "DP-0601", "DP-1301", "DP-4501", "DP-8001");	#"DP-1501","DP-1801",
	my $basetype_ref = \@basetype;
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
	push @combined, @dpb1;
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
	if ( $serotype =~ /DP/ ) {
		@residues = @dpb1;
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

	%ref = %dpb1; 

	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = %dpb1;
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
	my %tmp = %dpb1;
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key eq "DP-0101" ) {
			$ref{$key} = "DPB1\\*01";
		}
		elsif (( $key eq "DP-0201" ) || ( $key eq "DP-0202" )) {
			$ref{$key} = "DPB1\\*02";
		}
		elsif ( $key eq "DP-0301" ) {
			$ref{$key} = "DPB1\\*03";
		}
		elsif (( $key eq "DP-0401" ) || ( $key eq "DP-0402")) {
			$ref{$key} = "DPB1\\*04";
		}
		elsif ( $key eq "DP-0601" ) {
			$ref{$key} = "DPB1\\*06";
		}
		elsif ( $key eq "DP-1001" ) {
			$ref{$key} = "DPB1\\*10:";
		}
		elsif ( $key eq "DP-1101" ) {
			$ref{$key} = "DPB1\\*11:";
		}
		elsif ( $key eq "DP-1301" ) {
			$ref{$key} = "DPB1\\*13:";
		}
		elsif ( $key eq "DP-1501" ) {
			$ref{$key} = "DPB1\\*15:";
		}
		elsif ( $key eq "DP-1801" ) {
			$ref{$key} = "DPB1\\*18:";
		}
		elsif ( $key eq "DP-3401" ) {
			$ref{$key} = "DPB1\\*34:";
		}
		elsif ( $key eq "DP-4501" ) {
			$ref{$key} = "DPB1\\*45:";
		}
		elsif ( $key eq "DP-4601" ) {
			$ref{$key} = "DPB1\\*46:";
		}
		elsif ( $key eq "DP-6901" ) {
			$ref{$key} = "DPB1\\*69:";
		}
		elsif ( $key eq "DP-8001" ) {
			$ref{$key} = "DPB1\\*80:";
		}
		elsif ( $key eq "DP-13601" ) {
			$ref{$key} = "DPB1\\*136:";
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
	#$partial{ "DPB1*26:01" } = $seq2;
	$partial{ "DPB1*37:01" } = $seq2;
	$partial{ "DPB1*56:01" } = $seq2;
	$partial{ "DPB1*74:01" } = $seq2;
	$partial{ "DPB1*75:01" } = $seq2;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	return $known_cross_ref;
}

1;
