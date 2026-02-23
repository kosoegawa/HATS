#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: SUMCOUNT.pm 
# This module was developed to generate Summary table
# Separated from COUNT.pm
# last reviewed on February 20 2026

package SUMCOUNT;
use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
use Openfile;
use GROUP_SORT;
use ASSIGNED_SHORT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character


sub SUMMARY {
	my ( $csv_ref, $gene, $null_ref, $qallele_ref, $out_name ) = @_;
	my @list = Openfile::open_file_from_list($csv_ref);

	my $header = shift @list;
	my %full;
	my $full_ref = \%full;

	my $common_full = 0;
	my $inter_full = 0;
	my $wd_full = 0;
		
	my $common_sero = 0;
	my $inter_sero = 0;
	my $wd_sero = 0;
		
	my $common_short = 0;
	my $inter_short = 0;
	my $wd_short = 0;
		
	my $common_ins = 0;	# in silico
	my $inter_ins = 0;
	my $wd_ins = 0;

	foreach my $line ( @list ) {
		my @elements = split( ",", $line );
		my $type = $elements[1];
#		print $type . "\n";
		if (( $type eq "NULL" ) || ( $type eq "Questionable" )) {
#		if (( $type eq "N" ) || ( $type eq "Q" )) {
			next;
		}
		unless ( exists $full_ref->{ $type } ) {
			$full_ref->{ $type }->[0] = 0;		# FULL
			$full_ref->{ $type }->[1] = 0;		# Serotype
			$full_ref->{ $type }->[2] = 0;		# S
			$full_ref->{ $type }->[3] = 0;		# InSilico
		}

		if ( $elements[2] eq "U" ) {
			$full_ref->{ $type }->[0] = $full_ref->{ $type }->[0] + 1;
			next;
		}
		elsif ( $elements[2] eq "F" ) {
			$full_ref->{ $type }->[0] = $full_ref->{ $type }->[0] + 1;
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_full++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_full++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_full++;
			}
			next;
		}
		elsif ( $elements[2] eq "S" ) {
			$full_ref->{ $type }->[1] = $full_ref->{ $type }->[1] + 1;
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_sero++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_sero++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_sero++;
			}
			next;
		}
		elsif ( $elements[2] eq "I" ) {
			$full_ref->{ $type }->[2] = $full_ref->{ $type }->[2] + 1;
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_short++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_short++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_short++;
			}
			next;
		}
		elsif ( $elements[2] eq "InSilico" ) {
			$full_ref->{ $type }->[3] = $full_ref->{ $type }->[3] + 1;
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_ins++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_ins++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_ins++;
			}
			next;
		}
	}

	open(FILE, ">output/" . $gene . $out_name . $date . ".csv"); 
	print FILE "ANTIGEN,FULL,SEROTYPE,INCOMPLETE,InSilico\n";
	my $full = 0;
	my $sero = 0;
	my $short = 0;
	my $sc = 0;
	my $ins = 0;
	foreach my $type ( sort keys %full ) {
		print FILE $type . ",";
		print FILE $full_ref->{ $type }->[0] . ",";		# FULL
		print FILE $full_ref->{ $type }->[1] . ",";		# Serotype
		print FILE $full_ref->{ $type }->[2] . ",";		# S
		print FILE $full_ref->{ $type }->[3] . "\n";		# InSilico

		$full = $full + $full_ref->{ $type }->[0];
		$sero = $sero + $full_ref->{ $type }->[1];
		$short = $short + $full_ref->{ $type }->[2];
		$ins = $ins + $full_ref->{ $type }->[3];

	}
	print FILE "\n";

	print FILE "AssignedTotal," . $full . "," . $sero . "," . $short . "," . $ins . "\n";
	print FILE "CommonTotal," . $common_full . "," . $common_sero . "," . $common_short . "," . $common_ins . "\n";
	print FILE "IntermediateTotal," . $inter_full . "," . $inter_sero . "," . $inter_short . "," . $inter_ins . "\n";
	print FILE "WellDocumentedTotal," . $wd_full . "," . $wd_sero . "," . $wd_short . "," . $wd_ins . "\n";

	my $rare_full = $full - $common_full - $inter_full - $wd_full;
	my $rare_sero = $sero - $common_sero - $inter_sero - $wd_sero;
	my $rare_short = $short - $common_short - $inter_short - $wd_short;
	my $rare_ins = $ins - $common_ins - $inter_ins - $wd_ins; 
	print FILE "RareTotal," . $rare_full . "," . $rare_sero . "," . $rare_short . "," . $rare_ins . "\n\n";
	
	my $null_count = scalar keys %$null_ref;
	my $qallele_count = scalar keys %$qallele_ref;
	my $subtotal = $full + $sero + $short + $ins;
	print FILE "TOTAL EXPRESSED ALLELES," . $subtotal . "\n";
	print FILE "Null alleles," . $null_count . "\n";
	print FILE "Q alleles," . $qallele_count . "\n";
	my $total = $full + $sero + $short + $ins + $null_count + $qallele_count;
	print FILE "TOTAL ALLELES," . $total . "\n";

	close FILE;
}

1;
