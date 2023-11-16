#!/usr/bin/perl -w
#
# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: NullAllele.pm 
# This module was developed to capture Null alleles
# last reviewed on November 14 2023

package NullAllele;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub all {		# deal with remaining serotypes with strict mode
	my ($fasta_ref, $gene) = @_;

	my %nullAllele;
	my $nullAllele_ref = \%nullAllele;
	foreach my $head ( keys %$fasta_ref ) {
		if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
			my $allele = $1;
			if ( $allele =~ /N/ ) {	# exclude Null allele
				$nullAllele { $allele } = "NULL";
			}
		}
	}

	print "Number on NUll Alleles of " . $gene . ": " . scalar ( keys %nullAllele ) . "\n";
	open(FILE, ">output/NullAllele_" . $gene . "_" . $date . ".csv");
	foreach my $allele ( sort ( keys %nullAllele ) ) {
		print FILE $allele . "\n";
	}
	close FILE;

	return $nullAllele_ref;
}

1;
