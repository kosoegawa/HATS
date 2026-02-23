#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DPB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last reviewed, modified and documented on February 1 2026

package DPB1_INFO;
use strict;

my @dp1 = (56,57,69,84,85,86,96);
my @dp3 = (56,57,69,84,85,86);	# required for the alleles missing exon 3

my %group;
my %base;
my %dp1;
my %dp3;
$dp3{"DPB01"} = "HLA00514";	# DP-B1*01:01:01:01 
$dp3{"DPB0201"} = "HLA00517";	# DP-B1*02:01:02:01
$dp3{"DPB03"} = "HLA00520";	# DP-B1*03:01:01:01
$dp1{"DPB0401"} = "HLA00521";	# DP-B1*04:01:01:01
$dp3{"DPB10"} = "HLA00527";	# DP-B1*10:01:01:01
$dp3{"DPB15"} = "HLA00532";	# DP-B1*15:01:01:01
$dp1{"DPB18"} = "HLA00535";	# DP-B1*18:01:01:01
$dp3{"DPB46"} = "HLA00565";	# DP-B1*46:01:01
$dp3{"DPB0402"} = "HLA00522";	# DP-B1*04:02:01:01
$dp3{"DPB0202"} = "HLA00519";	# DP-B1*02:02:01:01
$dp3{"DPB06"} = "HLA00524";	# DP-B1*06:01:01:01
$dp3{"DPB13"} = "HLA00530";	# DP-B1*13:01:01:01
$dp3{"DPB45"} = "HLA00564";	# DP-B1*45:01
$dp3{"DPB80"} = "HLA00599";	# DP-B1*80:01
#$dp3{"DP-1101"} = "HLA00528";	# DP-B1*11:01:01:01, 69R is equivalent to 69K, this is DPB1
#$dp1{"DP-3401"} = "HLA00553";	# DP-B1*34:01:01:01, 69K is equivalent to 69R, this is DPB15
#$dp3{"DP-13601"} = "HLA07366";	# DP-B1*136:01, 69R is equivalent to 69K, this is equivalent to DPB45
#$dp3{"DP-6901"} = "HLA00588";	# DP-B1*69:01:01:01, 69R is equivalent to 69K, this is DPB3

$dp1{"DPB17"} = "HLA00534";	# DPB1*17:01:01:01 258 bp
$dp1{"DPB30"} = "HLA00549";	# DPB1*30:01:01:01 258 bp
$dp1{"DPB31"} = "HLA00550";	# DPB1*31:01:01:01 258 bp

$group{"DPB01"} = "DPB03"; $group{"DPB0201"} = "DPB03"; $group{"DPB03"} = "DPB03"; $group{"DPB0401"} = "DPB01"; $group{"DPB10"} = "DPB03";
$group{"DPB15"} = "DPB03"; $group{"DPB18"} = "DPB01"; $group{"DPB46"} = "DPB03"; $group{"DPB0402"} = "DPB03"; $group{"DPB0202"} = "DPB03";
$group{"DPB06"} = "DPB03"; $group{"DPB13"} = "DPB03"; $group{"DPB45"} = "DPB03"; $group{"DPB80"} = "DPB03";# $group{"DPB40"} = "DPB01";
# $group{"DPB11"} = "DPB03"; $group{"DPB34"} = "DPB01"; $group{"DPB13601"} = "DPB03"; $group{"DPB6901"} = "DPB03";
$group{"DPB17"} = "DPB01"; $group{"DPB30"} = "DPB01"; $group{"DPB31"} = "DPB01";
# $group{"DPB28"} = "DPB01"; $group{"DPB35001"} = "DPB01"; $group{"DPB46301"} = "DPB01";
$base{"DPB01"} = "DPB01"; $base{"DPB0201"} = "DPB0201"; $base{"DPB03"} = "DPB03"; $base{"DPB0401"} = "DPB0401"; $base{"DPB10"} = "DPB10";
$base{"DPB15"} = "DPB15"; $base{"DPB18"} = "DPB18"; $base{"DPB46"} = "DPB46"; $base{"DPB0402"} = "DPB0402"; $base{"DPB0202"} = "DPB0202";
$base{"DPB06"} = "DPB06"; $base{"DPB13"} = "DPB13"; $base{"DPB45"} = "DPB45"; $base{"DPB80"} = "DPB80";
#$base{"DPB1101"} = "DPB0101"; $base{"DPB3401"} = "DPB1501"; $base{"DPB13601"} = "DPB4501";  $base{"DPB6901"} = "DPB0301";
$base{"DPB17"} = "DPB17"; $base{"DPB30"} = "DPB30"; $base{"DPB31"} = "DPB31";
#$base{"DPB2801"} = "DPB0402"; $base{"DPB4001"} = "DPB0401"; $base{"DPB46301"} = "DPB1801"; $base{"DPB35001"} = "DPB1501";

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
	push @combined, @dp1;
	push @combined, @dp3;
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
	@residues = @dp1;
		if ( $serotype eq "DPB01" ) {
		@residues = @dp1;
	}
	elsif ( $serotype eq "DPB03" ) {
		@residues = @dp3;
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
		%ref = %dp1;
	}
	elsif ( $serotype eq "DPB03" ) {
		%ref = %dp3;
	}
	else {	# all together
		%ref = (%dp1,%dp3);
	}

	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dp1,%dp3);
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
	my %tmp = (%dp1,%dp3);
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
		#elsif ( $key eq "DPB1101" ) {
		#	$ref{$key} = "DPB1\\*11:";
		#}
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
		#elsif ( $key eq "DPB2801" ) {
		#	$ref{$key} = "DPB1\\*28:";
		#}
		elsif ( $key eq "DPB30" ) {
			$ref{$key} = "DPB1\\*30:";
		}
		elsif ( $key eq "DPB31" ) {
			$ref{$key} = "DPB1\\*31:";
		}
		#elsif ( $key eq "DPB3401" ) {
		#	$ref{$key} = "DPB1\\*34:";
		#}
		#elsif ( $key eq "DPB4001" ) {
		#	$ref{$key} = "DPB1\\*40:";
		#}
		elsif ( $key eq "DPB45" ) {
			$ref{$key} = "DPB1\\*45:";
		}
		elsif ( $key eq "DPB46" ) {
			$ref{$key} = "DPB1\\*46:";
		}
		#elsif ( $key eq "DPB6901" ) {
		#	$ref{$key} = "DPB1\\*69:";
		#}
		elsif ( $key eq "DPB80" ) {
			$ref{$key} = "DPB1\\*80:";
		}
		#elsif ( $key eq "DPB13601" ) {
		#	$ref{$key} = "DPB1\\*136:";
		#}
		#elsif ( $key eq "DPB35001" ) {
		#	$ref{$key} = "DPB1\\*350:";
		#}
		#		elsif ( $key eq "DPB46301" ) {
		#	$ref{$key} = "DPB1\\*463:";
		#}
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
	#$known_cross{ "DPB1101" } = 0;
	$known_cross{ "DPB13" } = 0;
	$known_cross{ "DPB15" } = 0;
	$known_cross{ "DPB17" } = 0;
	$known_cross{ "DPB18" } = 0;
	$known_cross{ "DPB30" } = 0;
	$known_cross{ "DPB31" } = 0;
	#$known_cross{ "DPB3401" } = 0;	# did not change to $dp3, SEROTYPE is same as DPB4001
	#$known_cross{ "DPB4001" } = 0;
	$known_cross{ "DPB45" } = 0;
	$known_cross{ "DPB46" } = 0;
	#$known_cross{ "DPB6901" } = 0;
	$known_cross{ "DPB80" } = 0;
	#$known_cross{ "DPB13601" } = 0;
	#$known_cross{ "DPB35001" } = 0;
	#$known_cross{ "DPB46301" } = 0;
	return $known_cross_ref;
}

1;
