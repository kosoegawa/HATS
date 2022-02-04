#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: DR5231ASSIGN.pm
# stringent assign
# This module was developed to convert DRB1 allele to DR52, DR53 and DR51 serotypes using strict mode
# last modified and documented on May 15 2020
#

package DR5231ASSIGN;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub all {		# deal with remaining serotypes with strict mode
	my ($fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_seq, $dr5231_ref ) = @_;
	# the larst argument $partial_seq can be ignored, if no partial reference seqence is present

	my %assigned;	# this is retured at the end of this module
	my $assigned_ref = \%assigned;
	foreach my $key ( keys %$ref_ref ) {
		my $target = "";
		foreach my $head ( keys %$fasta_ref ) {	# go through fastq
			if ( $head =~ /$ref_ref->{ $key }/ ) {	# matching HLA accessing number
				my $elements = scalar @$residues_ref;	# all residues
				for ( my $index = 0; $index < $elements; $index++ ) {
					my $position = $residues_ref->[ $index ] + $leader;	# residue position including leader peptide
					
					my $seq = "";
					$seq =  $fasta_ref->{ $head };
					unless ( $seq =~ /^M[A-Z]+/ ) {
						$seq = $partial_seq . $seq;
					}
					$target = $target . substr($seq, $position, 1);
					if ( $index != $elements - 1 ) {
						my $num = $residues_ref->[ $index + 1 ] - $residues_ref->[ $index ] - 1;	# add random AA between key residues
						$target = $target . "[A-Z]{$num}";
					}
				}
			}
		}

		my @alleles;
		foreach my $head ( keys %$fasta_ref ) {
			if ( $fasta_ref->{ $head } =~ /$target/ ) {
				if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
					my $allele = $1;
					unless (( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {	# exclude Null and Q allele
						push @alleles, $allele;
					}
				}
			}
		}
		open(FILE, ">output/" . $key . "_" . $dr5231_ref->{$key} . "_" . $date . ".csv");
		foreach my $allele ( sort @alleles ) {
			print FILE  $allele . "\n";
			$assigned{ $allele } = $dr5231_ref->{$key};
		}
		close FILE;
	}
	return $assigned_ref;
}


1;
