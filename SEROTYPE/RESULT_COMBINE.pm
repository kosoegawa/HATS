#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: RESULT_COMBINE.pm 
# This module was developed to generate Summary table
# last modified and documented on February 20 2026

package RESULT_COMBINE;
use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
use Openfile;
use GROUP_SORT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

sub COMBINE {		# deal with remaining serotypes with strict mode
	my ($results_dir) = @_;
	my @csv = glob($results_dir . "*.csv");

	my $header = "";
	my %combined;
	my $combined_ref = \%combined;
	my @alleles;
	my $alleles_ref = \@alleles;
	my $gene ="";
	my $db = "";
	my $count = 0;
	foreach my $csv ( sort @csv ) {
		print $csv . "\n";
		if ( $csv =~ /(\w+)_(Allele_Antigen_Table_IMGT_HLA_\d+\.\d+\.\d+)_/ ) {
			$gene = $gene . "_" . $1;
			$db = $2;
		}
		elsif ( $csv =~ /(\w+)_(Protein_Antigen_Table_IMGT_HLA_\d+\.\d+\.\d+)_/ ) {
			$gene = $gene . "_" . $1;
			$db = $2;
		}
	
		my @list = Openfile::open_file_from_list( $csv );
		$header = shift @list;
		foreach my $list ( @list ) {
			my @elements = split(",", $list);
			my $allele = shift @elements;

			if ( exists $combined_ref->{ $allele } ) {
				if ( $combined_ref->{ $allele }->[0] eq "U" ) {	
					delete $combined_ref->{ $allele };	# delete existing record
					for ( my $index = 0; $index < scalar @elements; $index++ ) {
						$combined_ref->{ $allele }->[$index] = $elements[ $index ];
					}
				}
			}
			else {
				push @alleles, $allele;
				for ( my $index = 0; $index < scalar @elements; $index++ ) {
					$combined_ref->{ $allele }->[$index] = $elements[ $index ];
				}
			}
		}
	}

	print $gene . "\n";
	open( FILE, ">COMBINED/combined" . $gene . "_" . $db . "_" . $date . ".csv" );
	print FILE $header . "\n";
	my $alleles_sorted_ref = GROUP_SORT::SORT( $alleles_ref );
	foreach my $allele ( @$alleles_sorted_ref ) {
		print FILE $allele . ",";
		my $limit = scalar @{$combined_ref->{ $allele }};
		for ( my $index = 0; $index < $limit; $index++ ) {
			if ( $index < $limit - 1 ) {
				print FILE $combined_ref->{ $allele }->[ $index ] . ",";
			}
			else {
				print FILE $combined_ref->{ $allele }->[ $index ] . "\n";
			}
		}
	}
	close FILE;

}


1;
