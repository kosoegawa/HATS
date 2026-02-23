#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: ASSIGNED_SHORT.pm 
# This module was developed to print table
# last reviewed, modified and documented on February 19 2026

package ASSIGNED_SHORT;
use strict;
use GROUP_SORT;
use COMBINE;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

my @broad = ("A9","A10","A19","A28","B5","B12","B14","B15","B16","B17","B21","B22","B40","B70","Cw3","DR2","DR3","DR5","DR6","DQ1","DQ3");

# print SHORT
sub PRINT {
	my ( $unassigned_ref, $short_ref ) = @_;
	print "SHORT\n";
	open(FILE, ">output/Assigned_SHORT_" . $date . ".csv");
	foreach my $allele (sort @$unassigned_ref) {
		if ( exists $short_ref->{ $allele } ) {
			print FILE $allele . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			for (my $index = 0; $index < $num; $index++) {
				print FILE $short_ref->{ $allele }->[$index];
				if ( $index < $num - 1) {
					print FILE ",";
				}
				else {
					print FILE "\n";
				}
			}
		}
		else {
			print FILE $allele . "\n";
		}
	}
	close FILE;
}

sub PRINT_RESIDUES {
	my ( $elements_ref,$gene,$residues_all_ref,$database ) = @_;

	my $antigen_ref = COMBINE::ANTIGEN();

	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$elements_ref;
	my $elements_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically
	my %residue; 
	open(FILE, ">output/" . $gene . "_DEP_" . $database . "_" . $date . ".csv");
	print FILE "Protein,AA/AN/BR,Qualifier,";
	# print residues
	for (my $index = 0; $index < scalar @$residues_all_ref; $index++) {
		print FILE $residues_all_ref->[ $index ];
		my $limit = scalar @$residues_all_ref - 1;
		if ( $index < $limit ) {
			print FILE ",";
		}
		else {
			print FILE "\n";
		}
	}
	
	foreach my $twoField ( @$elements_sorted_ref ) {		#go through all alleles
		unless ( exists $residue{ $twoField } ) {
			if ( exists $antigen_ref->{ $elements_ref->{ $twoField }->[0] } ) {
				print FILE $twoField . ",". $antigen_ref->{ $elements_ref->{ $twoField }->[0] } . "," . $elements_ref->{ $twoField }->[1];	# paragrarh mark is included in elements
			}
			else {
				print FILE $twoField . ",". $elements_ref->{ $twoField }->[0] . "," . $elements_ref->{ $twoField }->[1];	# paragrarh mark is included in elements
			}
			$residue{ $twoField } = 0;
		}
	}
	close FILE;
}


1;
