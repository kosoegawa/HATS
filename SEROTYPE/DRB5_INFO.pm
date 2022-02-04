#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DRB5_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on December 16 2019

package DRB5_INFO;
use strict;

my @dr51 = (9, 10, 11, 12, 13,71);

my %dr51;
my %group;
my %base;
$dr51{"DR5101"} = "HLA00915";		# DRB5*01:01:01:01
$dr51{"DR5102"} = "HLA00926";		# DRB5*02:02:01
$group{"DR5101"} = "DR51";
$group{"DR5102"} = "DR51";
$base{"DR5101"} = "DR51";
$base{"DR5102"} = "DR51";


my @subtype = ("DR5102");	# modify here if serotype modified

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
		
	return $partial_ref;
}

1;
