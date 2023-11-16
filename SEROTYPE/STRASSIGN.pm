#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: STRASSIGN.pm
# stringent assign
# This module was developed to convert HLA allele to HLA serotype using strict mode
# last modified and documented on October 25 2023
#

package STRASSIGN;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub all {		# deal with remaining serotypes with strict mode
	my ($fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref ) = @_;
	# the larst argument $partial_ref can be ignored, if no partial reference seqence is present

	my %assigned;	# this is retured at the end of this module
	my $assigned_ref = \%assigned;
	my %all;
	my %nullAllele;
	my %qAllele;
	foreach my $key ( keys %$ref_ref ) {	# go through each serotype
		my $target = "";
		foreach my $head ( keys %$fasta_ref ) {	# go through fasta
			if ( $head =~ /$ref_ref->{ $key }/ ) {	# matching HLA ID, not HLA allele name
				my $elements = scalar @$residues_ref;	# all residues
				for ( my $index = 0; $index < $elements; $index++ ) {
					my $position = $residues_ref->[ $index ] + $leader;	# residue position including leader peptide
					
					my $seq = "";
					$seq =  $fasta_ref->{ $head };
					unless ( $seq =~ /^M[A-Z]+/ ) {		#protein sequence does not start with M
						$seq = $partial_ref->{ $key } . $seq;	# add sequence
					}
					$target = $target . substr($seq, $position, 1);		# AA residue at target position
					if ( $index != $elements - 1 ) {
						my $num = $residues_ref->[ $index + 1 ] - $residues_ref->[ $index ] - 1;	# add random AA between key residues
						$target = $target . "[A-Z]{$num}";
					}
				}
			}
		}

		my @alleles;
		foreach my $head ( keys %$fasta_ref ) {
			if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
				my $allele = $1;
				$all{ $allele } = "";	# add HLA all gene specific allele
				if ( $allele =~ /[0-9]+N/ ) {	# exclude Null allele
					$nullAllele { $allele } = "";
				}
				elsif ( $allele =~ /[0-9]+Q/ ) {	# exclude Null allele
					$qAllele { $allele } = "";
				}
			}

			if ( $fasta_ref->{ $head } =~ /$target/ ) {
				if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
					my $allele = $1;
					unless (( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {	# exclude Null and Q allele
						push @alleles, $allele;
					}
				}
			}
		}
		open(FILE, ">output/" . $key . "_" . $date . ".csv");
		foreach my $allele ( sort @alleles ) {
			print FILE  $allele . "\n";
			$assigned{ $allele } = $key;
		}
		close FILE;
	}

	my @unassigned = ();
	my $alleleCount = scalar ( keys %all );
	print "Number of " . $gene . " alleles: " . $alleleCount . "\n";
	foreach my $all ( keys %all ) {
		unless (( exists $assigned{ $all }) || ( exists $nullAllele { $all } ) || ( exists $qAllele { $all } )) {
			push(@unassigned, $all);
		}
	}
	print "Assigned FULL: " . scalar ( keys %assigned ) . "\n";
	print "Unassigned FULL: " . scalar ( @unassigned ) . "\n";
	
	open (FILE, ">output/Stringent_Unclassified_" . $gene . "_" . $date . ".csv");
	foreach my $unassigned ( sort @unassigned ) {
		print FILE $unassigned . "\n";
	}
	close FILE;

	return $assigned_ref;
}


1;
