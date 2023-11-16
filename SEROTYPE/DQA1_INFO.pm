#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DQA1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on January 7 2021

package DQA1_INFO;
use strict;


my @dqa13 = (51,53,56,74,76);	#DQA1*01 & DQA1*03
my @dqa2456 = (25,51,52,53,74,75);	#DQA1*01 & DQA1*03

my %group;
my %base;
my %dqa13;
$dqa13{"DQA-01"} = "HLA00601";	# DQA1*01:01:01:01
$dqa13{"DQA-03"} = "HLA00608";	# DQA1*03:01:01
my %dqa2456;
$dqa2456{"DQA-02"} = "HLA00607";	# DQA1*02:01:01:01
$dqa2456{"DQA-04"} = "HLA00612";	# DQA1*04:01:01:01
$dqa2456{"DQA-05"} = "HLA00613";	# DQA1*05:01:01:01
$dqa2456{"DQA-06"} = "HLA00620";	# DQA1*06:01:01:01

$group{"DQA-01"} = "DQA01"; $group{"DQA-03"} = "DQA01";
$group{"DQA-02"} = "DQA02"; $group{"DQA-04"} = "DQA02"; $group{"DQA-05"} = "DQA02"; $group{"DQA-06"} = "DQA02";
$base{"DQA-01"} = "DQA01"; $base{"DQA-03"} = "DQA03";
$base{"DQA-02"} = "DQA02"; $base{"DQA-04"} = "DQA04"; $base{"DQA-05"} = "DQA05"; $base{"DQA-06"} = "DQA06";

my @subtype = ();	# modify here if serotype modified

sub DQA1 {
	my $gene = "DQA1";
	return $gene;
}

sub DQA1_LEADER {
	my $leader = 22;		# DQA1 specific
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
	my @basetype = ("DQA01","DQA02","DQA03","DQA04","DQA05","DQA06");
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
	push @combined, @dqa13;
	push @combined, @dqa2456;
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
	if (( $serotype eq "DQA01" ) || ( $serotype eq "DQA03" )) {
		@residues = @dqa13; 
	}
	elsif (( $serotype eq "DQA02" ) || ( $serotype eq "DQA04" ) || ( $serotype eq "DQA05" ) || ( $serotype eq "DQA06" )) {
		@residues = @dqa2456;
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

	if (( $serotype eq "DQA01" ) || ( $serotype eq "DQA03" )) {
		%ref = %dqa13; 
	}
	elsif (( $serotype eq "DQA02" ) || ( $serotype eq "DQA04") || ( $serotype eq "DQA05" ) || ( $serotype eq "DQA06" )) {
		%ref = %dqa2456;
	}
	else {
		%ref = (%dqa13, %dqa2456);
	}
	
	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dqa13, %dqa2456);
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
	my %tmp = (%dqa13, %dqa2456);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /DQA-01/ ) {
			$ref{$key} = "DQA1\\*01";
		}
		elsif ( $key =~ /DQA-02/ ) {
			$ref{$key} = "DQA1\\*02";
		}
		elsif ( $key =~ /DQA-03/ ) {
			$ref{$key} = "DQA1\\*03";
		}
		elsif ( $key =~ /DQA-04/ ) {
			$ref{$key} = "DQA1\\*04";
		}
		elsif ( $key =~ /DQA-05/ ) {
			$ref{$key} = "DQA1\\*05";
		}
		else {
			$ref{$key} = "DQA1\\*06";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 28;	#change the number of missing nucleotide
	$partial{ "general" } = $seq;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	return $known_cross_ref;
}

1;
