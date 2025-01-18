#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: combine.pl
# This is to merge all serotype tables in RESULTS directory
# last modified and documented on January 18 2025

use strict;
use lib 'SEROTYPE';
use Openfile;
use GROUP_SORT;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

my $results_dir = "RESULTS/";
my @csv = glob($results_dir . "*.csv");

# remove combined
my @tmp;
foreach my $csv ( @csv ) {
	
	if ( $csv =~ /combined/ ) {
		unlink $csv;
	}
	else {
		push @tmp, $csv;
	}
}

my $header = "";
my %combined;
my $combined_ref = \%combined;
my @alleles;
my $alleles_ref = \@alleles;
my $gene ="";
my $db = "";
my $count = 0;
foreach my $csv ( sort @tmp ) {
	print $csv . "\n";
	
	if ( $csv =~ /(\w+)_(Serotype_Table_IMGT_HLA_\d+\.\d+\.\d+)_/ ) {
		$gene = $gene . "_" . $1;
		$db = $2;
	}
	my @list = Openfile::open_file_from_list( $csv );
	$header = shift @list;
	foreach my $list ( @list ) {
		my @elements = split(",", $list);
		my $allele = shift @elements;

		if ( exists $combined_ref->{ $allele } ) {
			if ( $combined_ref->{ $allele }->[0] eq "UNA" ) {	
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

open( FILE, ">RESULTS/combined" . $gene . "_" . $db . "_" . $date . ".csv" );
print FILE $header . "\n";
	
	my $alleles_sorted_ref = GROUP_SORT::SORT( $alleles_ref );
	foreach my $allele ( @$alleles_sorted_ref ) {
		print FILE $allele . ",";
		my $limit = scalar @{$combined_ref->{ $allele }};
#		print $limit . "\n";
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

