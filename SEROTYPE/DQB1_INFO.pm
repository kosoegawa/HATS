#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DQB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on January 31 2022

package DQB1_INFO;
use strict;


my @dq1 = (84,85,86,87,89,90);	#89 or 90
my @dq2= (46,47,71,74,85,86,87);
my @dq3 = (45,57,74,84,85);
my @dq4 = (57,70,71,84,85,87);	# 57 was included to remove noise DRP1*01:34


my %group;
my %base;
my %dq1;
$dq1{"DQ-0501"} = "HLA00638";	#DQB1*05:01:01:01
$dq1{"DQ-0602"} = "HLA00646";	#DQB1*06:02:01:01
$dq1{"DQ-0604"} = "HLA00648";	#DQB1*06:04:01:01
$group{"DQ-0501"} = "DQ1"; $group{"DQ-0602"} = "DQ1"; $group{"DQ-0604"} = "DQ1";
$base{"DQ-0501"} = "DQ5"; $base{"DQ-0602"} = "DQ6"; $base{"DQ-0604"} = "DQ6";
my %dq2;
$dq2{"DQ-0201"} = "HLA00622";	#DQB1*02:01:01
$group{"DQ-0201"} = "DQ2";
$base{"DQ-0201"} = "DQ2";
my %dq3;
$dq3{"DQ-0301"} = "HLA00625";	#DQB1*03:01:01:01
$dq3{"DQ-0304"} = "HLA00630";	#DQB1*03:04:01
$dq3{"DQ-0302"} = "HLA00627";	#DQB1*03:02:01:01
$dq3{"DQ-0303"} = "HLA00629";	#DQB1*03:03:02:01
$group{"DQ-0301"} = "DQ3"; $group{"DQ-0304"} = "DQ3"; $group{"DQ-0302"} = "DQ3"; $group{"DQ-0303"} = "DQ3";
$base{"DQ-0301"} = "DQ7"; $base{"DQ-0304"} = "DQ7"; $base{"DQ-0302"} = "DQ8"; $base{"DQ-0303"} = "DQ9";
my %dq4;
$dq4{"DQ-0401"} = "HLA00636";	#DQB1*04:01:01:01
$group{"DQ-0401"} = "DQ4";
$base{"DQ-0401"} = "DQ4";

my @subtype = ("DQ-0604","DQ-0304","DQ-0302","DQ-0303");	# modify here if serotype modified

sub DQB1 {
	my $gene = "DQB1";
	return $gene;
}

sub DQB1_LEADER {
	my $leader = 31;		# DQB1 specific
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
	my @basetype = ("DQ5","DQ6","DQ2","DQ7","DQ8","DQ9","DQ4");
	my $basetype_ref = \@basetype;
}

sub BROAD {
	my %broad;
	foreach my $base ( keys %base ) {
		if (( $base{ $base } eq "DQ5" ) || ( $base{ $base } eq "DQ6" )) {
			$broad{ $base } = "DQ1";
		}
		elsif (( $base{ $base } eq "DQ7" ) || ( $base{ $base } eq "DQ8" ) || ( $base{ $base } eq "DQ9" )) {
			$broad{ $base } = "DQ3";
		}
		else {
			$broad{ $base } = $base{ $base };
		}
	}
	my $broad_ref = \%broad;
	return $broad_ref;
}

sub RESIDUES {
	my ( $serotype ) = @_;
	my @combined = ();
	push @combined, @dq1;
	push @combined, @dq2;
	push @combined, @dq3;
	push @combined, @dq4;
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
	if ( $serotype eq "DQ1" ) {
		@residues = @dq1; 
	}
	elsif ( $serotype eq "DQ2" ) {
		@residues = @dq2;
	}
	elsif ( $serotype eq "DQ3" ) {
		@residues = @dq3;
	}
	elsif ( $serotype eq "DQ4" ) {
		@residues = @dq4;
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

	if ( $serotype eq "DQ1" ) {
		%ref = %dq1; 
	}
	elsif ( $serotype eq "DQ2" ) {
		%ref = %dq2;
	}
	elsif ( $serotype eq "DQ3" ) {
		%ref = %dq3;
	}
	elsif ( $serotype eq "DQ4" ) {
		%ref = %dq4;
	}
	else {
		%ref = (%dq1, %dq2, %dq3, %dq4);
	}
	
	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dq1, %dq2, %dq3, %dq4);
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
	my %tmp = (%dq1, %dq2, %dq3, %dq4);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /DQ-02/ ) {
			$ref{$key} = "DQB1\\*02";
		}
		elsif ( $key =~ /DQ-03/ ) {
			$ref{$key} = "DQB1\\*03";
		}
		elsif ( $key =~ /DQ-04/ ) {
			$ref{$key} = "DQB1\\*04";
		}
		elsif ( $key =~ /DQ-05/ ) {
			$ref{$key} = "DQB1\\*05";
		}
		else {
			$ref{$key} = "DQB1\\*06";
		}
	}
	return $key_ref;

}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 25;	#change the number of missing nucleotide
		
	return $partial_ref;
}

1;
