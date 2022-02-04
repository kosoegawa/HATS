#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DPB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on July 31 2021

package DPB1_INFO;
use strict;

my @dpb1 = (56,57,69,82,84,85);	# 82 is included to be DPB1 specific

my %group;
my %base;
my %dpb1;
$dpb1{"DP0101"} = "HLA00514";	# DPB1*01:01:01:01 
$dpb1{"DP0201"} = "HLA00517";	# DPB1*02:01:02:01
$dpb1{"DP0301"} = "HLA00520";	# DPB1*03:01:01:01
$dpb1{"DP0401"} = "HLA00521";	# DPB1*04:01:01:01
$dpb1{"DP1001"} = "HLA00527";	# DPB1*10:01:01:01
$dpb1{"DP1501"} = "HLA00532";	# DPB1*15:01:01:01
$dpb1{"DP1801"} = "HLA00535";	# DPB1*18:01:01:01
$dpb1{"DP4601"} = "HLA00565";	# DPB1*46:01:01
$dpb1{"DP0402"} = "HLA00522";	# DPB1*04:02:01:01
$dpb1{"DP0202"} = "HLA00519";	# DPB1*02:02:01:01
$dpb1{"DP0601"} = "HLA00524";	# DPB1*06:01:01:01
$dpb1{"DP1301"} = "HLA00530";	# DPB1*13:01:01:01
$dpb1{"DP4501"} = "HLA00564";	# DPB1*45:01
$dpb1{"DP8001"} = "HLA00599";	# DPB1*80:01

$dpb1{"DP1101"} = "HLA00528";	# DPB1*11:01:01:01, 69R is equivalent to 69K, this is DP1
$dpb1{"DP3401"} = "HLA00553";	# DPB1*34:01:01:01, 69K is equivalent to 69R, this is DP15
$dpb1{"DP13601"} = "HLA07366";	# DPB1*136:01, 69R is equivalent to 69K, this is equivalent to DP45
$dpb1{"DP6901"} = "HLA00588";	# DPB1*69:01:01:01

$group{"DP0101"} = "DP1"; $group{"DP0201"} = "DP1"; $group{"DP0301"} = "DP1"; $group{"DP0401"} = "DP1"; $group{"DP1001"} = "DP1";
$group{"DP1501"} = "DP1"; $group{"DP1801"} = "DP1"; $group{"DP4601"} = "DP1"; $group{"DP0402"} = "DP1"; $group{"DP0202"} = "DP1";
$group{"DP0601"} = "DP1"; $group{"DP1301"} = "DP1"; $group{"DP1101"} = "DP1"; $group{"DP3401"} = "DP1"; $group{"DP4501"} = "DP1";
$group{"DP13601"} = "DP1"; $group{"DP8001"} = "DP1"; $group{"DP6901"} = "DP1";
$base{"DP0101"} = "DP01"; $base{"DP0201"} = "DP0201"; $base{"DP0301"} = "DP03"; $base{"DP0401"} = "DP0401";
$base{"DP1001"} = "DP10"; $base{"DP1501"} = "DP15"; $base{"DP1801"} = "DP18"; $base{"DP4601"} = "DP46";	
$base{"DP0402"} = "DP0402"; $base{"DP0202"} = "DP0202"; $base{"DP0601"} = "DP06"; $base{"DP1301"} = "DP13";
$base{"DP4501"} = "DP45";$base{"DP8001"} = "DP80";
$base{"DP1101"} = "DP01"; $base{"DP3401"} = "DP15"; $base{"DP13601"} = "DP45";  $base{"DP6901"} = "DP03";

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
	my @basetype = ("DP0101", "DP0201", "DP0301", "DP0401", "DP1001","DP1501", "DP1801","DP4601",
		"DP0402", "DP0202", "DP0601", "DP1301", "DP4501", "DP8001");
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
		if ( $key eq "DP0101" ) {
			$ref{$key} = "DPB1\\*01";
		}
		elsif (( $key eq "DP0201" ) || ( $key eq "DP0202" )) {
			$ref{$key} = "DPB1\\*02";
		}
		elsif ( $key eq "DP0301" ) {
			$ref{$key} = "DPB1\\*03";
		}
		elsif (( $key eq "DP0401" ) || ( $key eq "DP0402")) {
			$ref{$key} = "DPB1\\*04";
		}
		elsif ( $key eq "DP0601" ) {
			$ref{$key} = "DPB1\\*06";
		}
		elsif ( $key eq "DP1001" ) {
			$ref{$key} = "DPB1\\*10:";
		}
		elsif ( $key eq "DP1101" ) {
			$ref{$key} = "DPB1\\*11:";
		}
		elsif ( $key eq "DP1301" ) {
			$ref{$key} = "DPB1\\*13:";
		}
		elsif ( $key eq "DP1501" ) {
			$ref{$key} = "DPB1\\*15:";
		}
		elsif ( $key eq "DP1801" ) {
			$ref{$key} = "DPB1\\*18:";
		}
		elsif ( $key eq "DP3401" ) {
			$ref{$key} = "DPB1\\*34:";
		}
		elsif ( $key eq "DP4501" ) {
			$ref{$key} = "DPB1\\*45:";
		}
		elsif ( $key eq "DP4601" ) {
			$ref{$key} = "DPB1\\*46:";
		}
		elsif ( $key eq "DP6901" ) {
			$ref{$key} = "DPB1\\*69:";
		}
		elsif ( $key eq "DP8001" ) {
			$ref{$key} = "DPB1\\*80:";
		}
		elsif ( $key eq "DP13601" ) {
			$ref{$key} = "DPB1\\*136:";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 34;	#change the number of missing nucleotide
		
	return $partial_ref;
}

1;
