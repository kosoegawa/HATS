#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: BBw.pm 
# This module was developed to capture Bw4 and Bw6 alleles
# last reviewed on February 9 2026

package COPYRESULT;
use strict;
use File::Copy;


sub COPYRESULT {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("RESULTS/" . $gene . "*_Table_IMGT_HLA_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Allele_Antigen_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "RESULTS/") or die "Copy failed: $!";

	@csv = glob("COUNT/" . $gene . "_Allele_Antigen_Count_*");
	$csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Allele_Antigen_Count_" . $date . ".csv", "COUNT/") or die "Copy failed: $!";
}

sub COPYTWORESULT {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("TWORESULTS/" . $gene . "*_Table_IMGT_HLA_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Protein_Antigen_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "TWORESULTS/") or die "Copy failed: $!";

	@csv = glob("COUNT/" . $gene . "_Protein_Antigen_Count_*");
	$csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Protein_Antigen_Count_" . $date . ".csv", "COUNT/") or die "Copy failed: $!";
}

sub COPYPRACTICAL {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("PRACTICAL/" . $gene . "*_Table_IMGT_HLA_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Allele_Antigen_Practical_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "PRACTICAL/") or die "Copy failed: $!";

	@csv = glob("PRACTICAL_PRO/" . $gene . "*_Table_IMGT_HLA_*");
	$csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Protein_Antigen_Practical_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "PRACTICAL_PRO/") or die "Copy failed: $!";
}

sub COPYRESIDUE {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("RESIDUES/" . $gene . "_DEP_*");
	my @target = glob("RESIDUES/target_" . $gene . "_*");
	push(@csv,@target);
	my @target_lax = glob("RESIDUES/target_LAX_" . $gene . "_*");
	push(@csv,@target_lax);
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_DEP_" . $database . "_" . $date . ".csv", "RESIDUES/") or die "Copy failed: $!";
	copy("output/target_" . $gene . "_" . $date . ".csv", "RESIDUES/") or die "Copy gailed: $!";
	copy("output/target_LAX_" . $gene . "_" . $date . ".csv", "RESIDUES/") or die "Copy gailed: $!";
}

sub COPYLEGACY {
	my ( $gene, $database, $date ) = @_;
	my @csv = glob("LEGACY/" . $gene . "_Legacy_*");
	my $csvs = scalar @csv;
	if ( $csvs > 0 ) {
		unlink @csv;
	}
	copy("output/" . $gene . "_Legacy_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "LEGACY/") or die "Copy failed: $!";
	copy("output/" . $gene . "_Legacy_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "LEGACY/") or die "Copy failed: $!";


}

1;
