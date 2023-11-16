#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DPA1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on January 8 2021

package DPA1_INFO;
use strict;


my @dpa1 = (50,55,58 .. 60);	#DPA1*01 & DPA1*03

my %group;
my %base;
my %dpa1;
$dpa1{"DPA-01"} = "HLA00499";	# DPA1*01:03:01:01
$dpa1{"DPA-02"} = "HLA00504";	# DPA1*02:01:01:01

$group{"DPA-01"} = "DPA01"; $group{"DPA-02"} = "DPA01";
$base{"DPA-01"} = "DPA01"; $base{"DPA-02"} = "DPA02";

my @subtype = ();	# modify here if serotype modified

sub DPA1 {
	my $gene = "DPA1";
	return $gene;
}

sub DPA1_LEADER {
	my $leader = 30;		# DPA1 specific
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
	my @basetype = ("DPA01","DPA02");
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
	push @combined, @dpa1;
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
	if ( $serotype =~ /DPA/ ) {
		@residues = @dpa1;
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

	%ref = %dpa1; 

	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = %dpa1;
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
	my %tmp = %dpa1;
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /DPA-01/ ) {
			$ref{$key} = "DPA1\\*01";
		}
		elsif ( $key =~ /DPA-02/ ) {
			$ref{$key} = "DPA1\\*02";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 34;	#change the number of missing nucleotide
	$partial{ "general" } = $seq;
	my $seq2 = "N" x 41;	#change the number of missing nucleotide
	$partial{ "DPA1*01:03:02" } = $seq2;
		
	return $partial_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	return $known_cross_ref;
}

1;
