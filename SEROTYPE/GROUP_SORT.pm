#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: GROUP_SORT.pm 
# This module was developed to sort HLA alleles
# last reviewed on November 14 2023

package GROUP_SORT;
use strict;

sub SORT {
	my ( $combined_ref ) = @_;	# original allele list
	my %gene;
	my %group;	# first field group
	foreach my $allele ( @$combined_ref ) {
		my $gene = "";
		my $group = "";
		if ( $allele =~ /(\w+)\*(\d+):\d+/ ) {
			$gene = $1;
			$group = $2;
			$gene{ $gene } = $gene;
			if ( $group =~ /^0(\d+)/ ) {
				$group{ $1 } = $group;	
			}
			else {
				$group{ $group } = $group;
			}
		}
	}

	my @combined_alleles;
	my $combined_alleles_ref = \@combined_alleles;
	foreach my $gene ( sort keys %gene ) {
		foreach my $key ( sort { $a <=> $b } keys %group ) {	# go through group, numeric sort
			my @alleles;
			my $alleles_ref = \@alleles;
			foreach my $allele ( @$combined_ref ) {
				my $locus = "";
				my $group = "";
				if ( $allele =~ /(\w+)\*(\d+):\d+/ ) {
					$locus = $1;
					$group = $2;
				}

				if (( $gene eq $locus ) && ( $group{ $key } eq $group )) {	# first field match
					push @alleles, $allele;
				}
			}
			my $tmp = SECOND_SORT( $alleles_ref);
			push @combined_alleles, @$tmp;
		}

	}
	return $combined_alleles_ref;
}

# second field sort
sub SECOND_SORT {
	my ( $combined_ref ) = @_;
	my %group;
	foreach my $allele ( @$combined_ref ) {
		my $group = "";
		if ( $allele =~ /\w+\*\d+:(\d+)/ ) {
			$group = $1;
			if ( $group =~ /^0(\d+)/ ) {
				$group{ $1 } = $group;	
			}
			else {
				$group{ $group } = $group;
			}
		}
	}

	my @combined_alleles;
	my $combined_alleles_ref = \@combined_alleles;
	foreach my $key ( sort { $a <=> $b } keys %group ) {	# go through group
		my @alleles;
		my $alleles_ref = \@alleles;
		foreach my $allele ( @$combined_ref ) {
			my $group = "";
			if ( $allele =~ /\w+\*\d+:(\d+)/ ) {
				$group = $1;
			}

			if ( $group{ $key } eq $group ) {	# second field match
				push @alleles, $allele;
			}
		}
		my $tmp = THIRD_SORT( $alleles_ref);
		push @combined_alleles, @$tmp;
	}
	my $size = scalar @combined_alleles;
	if ( $size == 0 ) {
		$combined_alleles_ref = $combined_ref;
	}
	return $combined_alleles_ref;
}

# third field sort
sub THIRD_SORT {
	my ( $combined_ref ) = @_;	# alleles (first, second, field allele names are identical)
	my %group;
	my $two_field = "";
	my $two_count = 0;
	foreach my $allele ( @$combined_ref ) {
		my $group = "";
		if ( $allele =~ /\w+\*\d+:\d+:(\d+)/ ) {
			$group = $1;
			if ( $group =~ /^0(\d+)/ ) {
				$group{ $1 } = $group;	
			}
			else {
				$group{ $group } = $group;
			}
		}
		elsif ( $allele =~ /\w+\*\d+:\d+/ ) {
			$two_field = $allele;
			$two_count++;
		}
	}

	my @combined_alleles;
	if ( $two_count > 0 ) {
		push @combined_alleles, $two_field;
	}
	my $combined_alleles_ref = \@combined_alleles;
	foreach my $key ( sort { $a <=> $b } keys %group ) {	# go through group
		my @alleles;
		my $alleles_ref = \@alleles;
		foreach my $allele ( @$combined_ref ) {
			my $group = "";
			if ( $allele =~ /\w+\*\d+:\d+:(\d+)/ ) {
				$group = $1;
			}

			if ( $group{ $key } eq $group ) {	# first field match
				push @alleles, $allele;
			}
		}
		my $tmp = FOURTH_SORT( $alleles_ref);
		push @combined_alleles, @$tmp;
	}
	my $size = scalar @combined_alleles;
	if ( $size == 0 ) {
		$combined_alleles_ref = $combined_ref;
	}
	return $combined_alleles_ref;
}

# fourth-field sort
sub FOURTH_SORT {
	my ( $combined_ref ) = @_;	# alleles (first, second, thrid field allele names are identical)
	my %group;
	my $third_field = "";
	my $third_count = 0;
	foreach my $allele ( @$combined_ref ) {
		my $group = "";
		if ( $allele =~ /\w+\*\d+:\d+:\d+:(\d+)/ ) {
			$group = $1;
			if ( $group =~ /^0(\d+)/ ) {
				$group{ $1 } = $group;	
			}
			else {
				$group{ $group } = $group;
			}
		}
		elsif ( $allele =~ /\w+\*\d+:\d+:\d+/ ) {
			$third_field = $allele;
			$third_count++;
		}
	}

	my @alleles;	# sorted allele
	if ( $third_count > 0 ) {
		push @alleles, $third_field;
	}
	my $alleles_ref = \@alleles;
	foreach my $key ( sort { $a <=> $b } keys %group ) {	# go through group
		foreach my $allele ( @$combined_ref ) {
			my $group = "";
			if ( $allele =~ /\w+\*\d+:\d+:\d+:(\d+)/ ) {
				$group = $1;
			}

			if ( $group{ $key } eq $group ) {	# first field match
				push @alleles, $allele;
			}
		}
	}
	my $size = scalar @alleles;
	if ( $size == 0 ) {
		$alleles_ref = $combined_ref;
	}
	return $alleles_ref;	# return sorted allele list
}

1;
