#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: Bw46ASSIGN.pm
# stringent assign
# This module was developed to convert HLA allele to HLA serotype using strict mode
# last modified and documented on March 12 2020
#

package Bw46ASSIGN;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub all {		# deal with remaining serotypes with strict mode
	#my ($fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $bw_ref ) = @_;
	my ($fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $bw_ref, $partial_ref ) = @_;
	# the larst argument $partial_ref can be ignored, if no partial reference seqence is present

	my %assigned;	# this is retured at the end of this module
	my $assigned_ref = \%assigned;
	foreach my $key ( keys %$ref_ref ) {
		my $target = "";
		foreach my $head ( keys %$fasta_ref ) {	# go through fasta
			if ( $head =~ /$ref_ref->{ $key }/ ) {	# matching HLA accessing number
				my $elements = scalar @$residues_ref;	# all residues
				for ( my $index = 0; $index < $elements; $index++ ) {
					my $position = $residues_ref->[ $index ] + $leader;	# residue position including leader peptide
					
					my $seq = "";
					if ( exists $partial_ref->{ $key } ) {	# add missing sequence if only partial sequence is available
						$seq = $partial_ref->{ $key } . $fasta_ref->{ $head };
					}
					else {
						$seq =  $fasta_ref->{ $head };
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
		open(FILE, ">output/" . $key . "_" . $bw_ref->{$key} . "_" . $date . ".csv");
		foreach my $allele ( sort @alleles ) {
			print FILE  $allele . "\n";
			$assigned{ $allele } = $bw_ref->{$key};
		}
		close FILE;
	}
	return $assigned_ref;
}


sub Bw {		# deal with remaining serotypes with strict mode
	#my ($fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $bw_ref ) = @_;
	my ($fasta_ref, $leader, $ref_ref, $residues_ref, $bw_ref, $partial_ref ) = @_;
	# the larst argument $partial_ref can be ignored, if no partial reference seqence is present

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
					if ( exists $partial_ref->{ $key } ) {	# add missing sequence if only partial sequence is available
						$seq = $partial_ref->{ $key } . $fasta_ref->{ $head };
					}
					else {
						$seq =  $fasta_ref->{ $head };
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
				if ( $head =~ /HLA:\S+ (\w+\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
					my $allele = $1;
					unless (( $allele =~ /[0-9]+N/ ) || ( $allele =~ /[0-9]+Q/ )) {	# exclude Null and Q allele
						push @alleles, $allele;
					}
				}
			}
		}
		open(FILE, ">output/" . $key . "_" . $bw_ref->{$key} . "_" . $date . ".csv");
		foreach my $allele ( sort @alleles ) {
			print FILE  $allele . "\n";
			$assigned{ $allele } = $bw_ref->{$key};
		}
		close FILE;
	}
	return $assigned_ref;
}

1;
