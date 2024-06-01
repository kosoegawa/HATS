#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: BBw.pm 
# This module was developed to capture Bw4 and Bw6 alleles
# last reviewed on November 14 2023

package COPYRESULT;
use strict;
use File::Copy;


sub COPYRESULT {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("RESULTS/" . $gene . "_Serotype_Table_IMGT_HLA_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "RESULTS/") or die "Copy failed: $!";
}

sub COPYTWORESULT {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("TWORESULTS/" . $gene . "_TwoField_Serotype_Table_IMGT_HLA_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "TWORESULTS/") or die "Copy failed: $!";
}

sub COPYRESIDUE {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("RESIDUES/" . $gene . "_Allele_Residues_*");
	my @target = glob("RESIDUES/target_" . $gene . "_*");
	push(@csv,@target);
	my @target_lax = glob("RESIDUES/target_LAX_" . $gene . "_*");
	push(@csv,@target_lax);
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Allele_Residues_" . $database . "_" . $date . ".csv", "RESIDUES/") or die "Copy failed: $!";
	copy("output/target_" . $gene . "_" . $date . ".csv", "RESIDUES/") or die "Copy gailed: $!";
	copy("output/target_LAX_" . $gene . "_" . $date . ".csv", "RESIDUES/") or die "Copy gailed: $!";
}

1;
