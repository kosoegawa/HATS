#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: ASSIGN.pm 
# This module was developed to convert HLA allele to HLA serotype using relaxed mode
# last reviewed, modified and documented on October 6 2023

package ASSIGN;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub ASSIGN {		# deal with remaining serotypes with strict mode
	my ( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $known_cross_ref ) = @_;
	#, $known_cross_ref_in
	# $fasta_ref: original fasta file, hla_prot.fasta
	# $assigned_ref: assigned alleles
	# $gene: HLA gene
	# $leader specifies the length of leader peptide - 1
	# $ref_ref: reference alleles, hash reference
	# $residues_ref: key residue positions in array reference
	# partial_ref is an optional if no partial sequence is used as reference
	my %known_cross = %$known_cross_ref;

	my %assigned;	#assigned alleles
	my $assigned2_ref = \%assigned;
	my %all;	# all alleles in a specific gene
	my %allHLA;	# all alleles in the hla_prot.fasta file
	foreach my $key ( keys %$ref_ref ) {	# $key indicates serotype, e.g. A-0201
		print $key . "\n";
		my $target = "";	# target specifies AA sequence

		# This loop captures target AA sequence
		foreach my $head ( keys %$fasta_ref ) {		# go through original hla_prot.fasta
			if ( $head =~ /$ref_ref->{ $key }/ ) {
				my $elements = scalar @$residues_ref;	# determine the number of key residues
				for ( my $index = 0; $index < $elements; $index++ ) {	# go through residues
					my $position = $residues_ref->[ $index ] + $leader;	# residue position in AA sequence

					my $seq = "";
					$seq =  $fasta_ref->{ $head };
					unless ( $seq =~ /^M[A-Z]+/ ) {
						$seq = $partial_ref->{ $key } . $seq;
					}
					$target = $target . substr($seq, $position, 1);
					if ( $index != $elements - 1 ) {
						# calculate the number of AA between key residues => any AA sequence, but length is important
						my $num = $residues_ref->[ $index + 1 ] - $residues_ref->[ $index ] - 1;
						$target = $target . "[A-Z]{$num}";	# add any sequence in a specific length
					}
				}
			}
		}

		my @alleles;
		foreach my $head ( keys %$fasta_ref ) {		# go through original hla_prot.fasta
			# capture all alleles
			if ( $head =~ /HLA:\S+ (\S+\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
				my $allele = $1;
				# exclude already assigned, Null and Q alleles
				unless (( exists $assigned_ref->{ $allele } ) || ( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {
					$allHLA{ $allele } = "";	# add HLA allele
				}
			}
			# capture gene specific alleles
			if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
				my $allele = $1;
				# exclude already assigned, Null and Q alleles
				unless (( exists $assigned_ref->{ $allele } ) || ( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {
					$all{ $allele } = "";	# add gene specific HLA allele
				}
			}

			if ( $fasta_ref->{ $head } =~ /$target/ ) {	# include target
				if ( $head =~ /HLA:\S+ (\S+\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
					my $allele = $1;
					# exclude already assigned, Null and Q alleles
					unless (( exists $assigned_ref->{ $allele } ) || ( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {
						push @alleles, $allele;		# add allele containing target
					}
					if ( exists $assigned_ref->{ $allele } ) {
						if ( exists($known_cross{ $assigned_ref->{ $allele }} )) {
							next;
						}
						unless ( ( $assigned_ref->{ $allele } =~ /$key/ )) {
							print "Cross Reactive:" . $allele . "\t" . $assigned_ref->{ $allele } . "\n";
						}
					}
				}
			}
		}
		open(FILE, ">output/" . $key .  "_LAX_" . $date . ".csv");
		foreach my $allele ( sort @alleles ) {
			print FILE  $allele . "\n";
			$assigned{ $allele } = $key .  "_LAX";
		}
		close FILE;
	}

	my @unassigned = ();
	my $alleleCount = scalar ( keys %all);
	foreach my $all ( keys %all ) {
		unless ( exists $assigned{ $all }) {
			push(@unassigned, $all);
		}
	}
	print "Unassigned " . $gene . ": " . scalar ( @unassigned ) . "\n";
	
	#combine both assigned alleles
	foreach my $allele ( keys %$assigned_ref ) {
		$assigned{ $allele } = $assigned_ref->{ $allele };
	}
	
	return $assigned2_ref;
}

sub CROSS {	# populate cross-reactive alleles	
	my ( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref ) = @_;
	#, $known_cross_ref_in
	# $fasta_ref: original fasta file, hla_prot.fasta
	# $assigned_ref: assigned alleles
	# $gene: HLA gene
	# $leader specifies the length of leader peptide - 1
	# $ref_ref: reference alleles, hash reference
	# $residues_ref: key residue positions in array reference
	# partial_ref is an optional if no partial sequence is used as reference
	# corss-reactive alleles

	foreach my $key ( keys %$ref_ref ) {	# $key indicates serotype, e.g. A2
#		print $key . " CROSS\n";
		my $target = "";	# target specifies AA sequence

		# This loop captures target AA sequence
		foreach my $head ( keys %$fasta_ref ) {		# go through original hla_prot.fasta
			if ( $head =~ /$ref_ref->{ $key }/ ) {
				my $elements = scalar @$residues_ref;	# determine the number of key residues
				for ( my $index = 0; $index < $elements; $index++ ) {	# go through residues
					my $position = $residues_ref->[ $index ] + $leader;	# residue position in AA sequence

					my $seq = "";
					$seq =  $fasta_ref->{ $head };
					unless ( $seq =~ /^M[A-Z]+/ ) {
						$seq = $partial_ref->{ $key } . $seq;
					}
					$target = $target . substr($seq, $position, 1);
					if ( $index != $elements - 1 ) {
						# calculate the number of AA between key residues => any AA sequence, but length is important
						my $num = $residues_ref->[ $index + 1 ] - $residues_ref->[ $index ] - 1;
						$target = $target . "[A-Z]{$num}";	# add any sequence in a specific length
					}
				}
			}
		}

		foreach my $head ( keys %$fasta_ref ) {		# go through original hla_prot.fasta
			if ( $fasta_ref->{ $head } =~ /$target/ ) {	# include target
				if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
					my $allele = $1;
					if ( exists $assigned_ref->{ $allele } ) {
						unless ( $assigned_ref->{ $allele } =~ /$key/ ) {
							if ( exists $cross_ref->{ $allele } ) {
								my $count = scalar @{$cross_ref->{ $allele }};
								$cross_ref->{ $allele }->[ $count ] = $key;
							}
							else {
								$cross_ref-> { $allele }->[0] = $key
							}
						}
					}
				}
			}
		}
	}

	return $cross_ref;
}


sub UNASSIGNED {
	my ( $fasta_ref, $assigned_ref, $gene ) = @_;
	my @unassigned;
	my $unassigned_ref = \@unassigned;
	foreach my $head ( keys %$fasta_ref ) {
		if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
			my $allele = $1;
			unless (( exists $assigned_ref->{ $allele }) || ( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {
				push @unassigned, $allele;	# add HLA DQB1 allele
			}
		}
	}

	open (FILE, ">output/unassigned_" . $gene . "_" . $date . ".csv");
	foreach my $unassigned (sort @unassigned) {
		print FILE $unassigned . "\n";
	}
	close FILE;

	return $unassigned_ref;
}


# missing one key resiidue
sub SHORT {	
	my ( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref ) = @_;
	# $short_ref: keep adding short serotype

	my $res_elements = scalar @$residues_ref;
	my $last_pos = $res_elements - 1;
	my @positions = (0..$last_pos);
	foreach my $pos ( @positions ) {	# go through array position
		my @tmp = @$residues_ref;
		my $removed = 0;
		if ( $res_elements > 4 ) {	# do short for the sero types that has 5 or more elements
			$removed = splice @tmp, $pos, 1;	# remove an element
#			print $removed . "\n";
		}
		my $elements = scalar @tmp;

		foreach my $key ( keys %$ref_ref ) {
			my $target = "";	# define target
			foreach my $head ( keys %$fasta_ref ) {	# go through fasta
				if ( $head =~ /$ref_ref->{ $key }/ ) {
					for ( my $index = 0; $index < $elements; $index++ ) {
						my $position = $tmp[ $index ] + $leader;

						my $seq = "";
						$seq =  $fasta_ref->{ $head };
						unless ( $seq =~ /^M[A-Z]+/ ) {
							$seq = $partial_ref->{ $key } . $seq;
						}
						$target = $target . substr($seq, $position, 1);
						if ( $index != $elements - 1 ) {
							my $num = $tmp[ $index + 1 ] - $tmp[ $index ] - 1;
							$target = $target . "[A-Z]{$num}";
						}
					}
				}
			}
			unless ( $target eq "" ) {
			#	print $target . "\n";
				my @alleles;
				foreach my $head ( keys %$fasta_ref ) {
					if ( $fasta_ref->{ $head } =~ /$target/ ) {
						if ( $head =~ /HLA:\S+ (\S+\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
							my $allele = $1;
							# exclude already assigned, Null and Q alleles
							unless (( exists $assigned_ref->{ $allele } ) || ( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {
								push @alleles, $allele;
							}
						}
					}
				}
				foreach my $allele ( sort @alleles ) {
					unless ( exists $short_ref->{ $allele } ) {	# does not exist
						my $value = $key . "_". $removed;
						# when 58E does not exist, it is not considered as DR11 group
						unless (( $value =~ /DR11\d*_58/ ) || ( $value eq "A265_144")) {	
							$short_ref->{ $allele }->[0] = $value;
						}
					}
					else {		# exists
						my $num = scalar @{$short_ref->{ $allele }};
						my $value = $key . "_". $removed;
						unless (( $value =~ /DR11\d*_58/ ) || ( $value eq "A265_144")) {	
							$short_ref->{ $allele }->[ $num ] = $value;
						}
					}
				}
			}
		}
	}
	return $short_ref;
}

1;
