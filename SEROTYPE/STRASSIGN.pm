#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: STRASSIGN.pm
# stringent assign
# This module was developed to convert HLA allele to HLA serotype using strict mode
# last modified and documented on February 22 2026
#

package STRASSIGN;
use strict;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
#my $date = `date +%F`;          # invoke bash date command
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

					my $aa = substr($seq, $position, 1);		# AA residue at target position
					if (( $gene eq "B" ) && ( $residues_ref->[ $index ] == 67 )) {
						if (( $aa eq "S" ) || ( $aa eq "C" )) {
							$target = $target . "[S,C]";
						}
						elsif (( $aa eq "Y" ) || ( $aa eq "F" )) {
							$target = $target . "[Y,F]";
						}
						else {
							$target = $target . $aa;		# AA residue at target position
						}
					}
					elsif (( $gene eq "B" ) && ( $residues_ref->[ $index ] == 167) && (( $aa eq "S" ) || ( $aa eq "G" ))) {
						$target = $target . "[S,G]";
					}
					elsif (( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 67 ) && (( $aa eq "I" ) || ( $aa eq "L" ))) {
						$target = $target . "[I,L]";
					}
					elsif (((( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 71 )) || (( $gene =~ /DPB/ ) && ( $residues_ref->[ $index ] == 69 ))) && (( $aa eq "K" ) || ( $aa eq "R" ))) {
						$target = $target . "[K,R]";
					}
					elsif (( $gene =~ /DRB/ ) && (!( $key eq "DR0301")) && (!( $key eq "DR0302" )) && ( $residues_ref->[ $index ] == 47 ) && (( $aa eq "F" ) || ( $aa eq "Y" ))) {
						$target = $target . "[F,Y]";
					}
					elsif (($gene eq "DQB1") && ($residues_ref->[ $index ] == 57) && (!( $key eq "DQ0302")) && (!( $key eq "DQ0303" ))) {
						if ((( $aa eq "D" ) || ( $aa eq "V" ) || ( $aa eq "S" )) || ( $aa eq "A" )) {
							$target = $target . "[D,V,S,A]";
						}
					}
					elsif (( $gene eq "DPB1" ) && ( $residues_ref->[ $index ] == 84 ) && (( $aa eq "G" ) || ( $aa eq "V" ))) {
						$target = $target . "[G,V]";
					}
					else {
						$target = $target . $aa;		# AA residue at target position
					}

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
				unless ( $gene =~ /DRB/ ) {
					if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
						my $allele = $1;
						unless (( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {	# exclude Null and Q allele
							push @alleles, $allele;
						}
					}
				}
				else {
					if ( $head =~ /HLA:\S+ (DRB[1,3,4,5]\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
						my $allele = $1;
						unless (( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {	# exclude Null and Q allele
							push @alleles, $allele;
						}
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
