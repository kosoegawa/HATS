#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: READ_MSF.pm
# Read msn file and generate stuffer sequence
# last reviewed on April 5 2026

package READ_MSF;
use strict;
use lib 'SEROTYPE';
use Openfile;
use GROUP_SORT;

sub MISSING_SEQ {
	my ( $gene, $null_ref, $qallele_ref, $ref_allele, $msf_ref ) = @_;

	my @list = Openfile::open_file_from_list( $msf_ref );

	my %allele_seq;
	my $allele_seq_ref = \%allele_seq;
	my @combined;
	my $combined_ref = \@combined;
	foreach my $line ( @list ) {
		if ( $line =~ /Name:/ ) {
			next;
		}
		elsif (( exists $allele_seq{ $ref_allele } ) && $line =~ /$ref_allele/) { #exist from the loop in the second incidence of the reference
			last;
		}

		if ( $line =~ /($gene\*\d+:\d+:*\d*:*\d*\w*)\s+(.+)/ ) {	# allele name, this should filter out extra gene seq
			my $allele = $1;
			my $seq = $2;
			$seq =~ s/ //g;	# remove space
			$seq =~ s/\.//g;	#remove .

			unless (( exists $null_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
				if ( exists $allele_seq{ $allele } ) {
					next;
				}
				else {	# capture only the first incidence
					push @combined, $allele;
					$allele_seq{ $allele } = $seq;
				}
			}
		}
	}

	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically
	my $ref_len = length( $allele_seq{ $ref_allele } );

	my %missing_seq;
	my $missing_seq_ref = \%missing_seq;
	foreach my $allele ( @$alleles_sorted_ref ) {
		my $difference = $ref_len - length ( $allele_seq{ $allele} );
		$missing_seq{ $allele } = "X" x $difference; # added X for the length of the missing sequence
	}

	return $missing_seq_ref;
}



1;
