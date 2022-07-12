#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: ORGANIZE.pm 
# This module was developed to organize fasta file
# CIWD was developed to add CIWD information in the final output
# last modified and documented on February 24 2022

package ORGANIZE;
use strict;
use Openfile;

sub fasta {		
	my ( $arg ) = @_;
	my @list = Openfile::open_file_from_list( $arg );

	my %fasta;
	my $fasta_ref = \%fasta;
	my $key = "";
	my $seq = "";
	foreach my $line (@list) {
		my $_ = $line;
		if (/>/) {
			$key = $line;
			$seq = "";
		}
		else {
			$seq = $seq . $line;
		}
		$fasta{ $key } = $seq;
	}
	return $fasta_ref;
}

sub CIWD {
	my ( $gene ) = @_;
	my @list = ();
	my $file = "CIWD/" . $gene . "_CIWD.csv";
	
	my %ciwd;
	my $ciwd_ref = \%ciwd;

	if ( -e $file ) {
		@list = Openfile::open_file_from_list( $file );
	}

	foreach my $line (@list) {
		my @order = ( "C", "I", "W", "D" );

		unless (( $line =~ /N/ ) || ( $line =~ /\d+Q/ )) {
			if ( $line =~ /($gene\*\d+:\d+).*,([CIWD]+)/ ) {
				my $allele = $1;
				my $score = $2;
				unless ( exists $ciwd{ $allele } ) {
					$ciwd{ $allele } = $score;
				}
				else {
					unless ( $ciwd{ $allele } =~ /$score/ ) {	# avoid duplicate
						my $string = $ciwd{ $allele } . $score;
						my @chars = split //, $string;	# split by characters
						my $ordered = "";
						foreach my $order ( @order ) {
							foreach my $chars ( @chars ) {
								if ( $order eq $chars ) {
									$ordered = $ordered . $order;
									last;
								}
							}
						}
						if ( $ordered =~ /C/ ) {
							$ordered = "C";
						}
						elsif ( $ordered =~ /I/ ) {
							$ordered = "I";
						}
						$ciwd{ $allele } = $ordered;
					}
				}
			}
		}
	}
	return $ciwd_ref;
}

sub CWD {
	my ( $gene ) = @_;
	my @list;
	my $file = "CWD2/cwd200_alleles_noheader.txt";
	@list = Openfile::open_file_from_list( $file );

	my %cwd;
	my $cwd_ref = \%cwd;

	foreach my $line (@list) {
		unless ( $line =~ /(^$gene)\s+HLA\d+\s+/ ) {
			next;
		}
		my @order = ( "C", "WD" );

		unless (( $line =~ /\d+N/ ) || ( $line =~ /\d+Q/ )) {
			if ( $line =~ /$gene\s+HLA\d+\s+($gene\*\d+:\d+):*\d*:*\d*\s+\w+\s+\w+\s+([CWD]+)\s+/ ) {
				my $allele = $1;
				my $score = $2;
				unless ( exists $cwd{ $allele } ) {
					$cwd{ $allele } = $score;
				}
				else {
					unless ( $cwd{ $allele } =~ /$score/ ) {	# avoid duplicate
						my $string = $cwd{ $allele } . $score;
						my @chars = split //, $string;	# split by characters
						my $ordered = "";
						foreach my $order ( @order ) {
							foreach my $chars ( @chars ) {
								if ( $order eq $chars ) {
									$ordered = $ordered . $order;
									last;
								}
							}
						}
						if ( $ordered =~ /C/ ) {
							$ordered = "C";
						}
						$cwd{ $allele } = $ordered;
					}
				}
			}
		}
	}
	return $cwd_ref;
}

sub EURCWD {
	my ( $gene ) = @_;
	my @list;
	my $file = "EURCWD/EURCWD.csv";
	@list = Openfile::open_file_from_list( $file );

	my %cwd;
	my $cwd_ref = \%cwd;

	foreach my $line (@list) {
		if ($line =~ /($gene\*\d+:\d+),([CWD]+)/ ) {
			$cwd{$1} = $2;
		}
	}
	return $cwd_ref;
}

sub ALLELE {
	my ( $arg, $gene ) = @_;
	my @list = Openfile::open_file_from_list( $arg );
	
	my %alleles;
	my $alleles_ref = \%alleles;
	foreach my $line (@list) {
		unless (( $line =~ /N/ ) || ( $line =~ /Q/ )) {
			if ( $line =~ /($gene\*\d+:\d+:\d*:\d*\w*),/ ) {
				$alleles{ $1 } = 0;
			}
		}
	}
	return $alleles_ref;
}

sub PGROUP {
	my ( $file, $gene )= @_;
	my @list = Openfile::open_file_from_list( $file );

	my %group;
	my $group_ref = \%group;
	foreach my $line ( @list ) {
		if ( $line =~ /#/ ) {
			next;
		}
		elsif ( $line =~ /^$gene\*/ ) {
			my @element = split(";", $line);
			my @group = split("/", $element[1]);
			foreach my $allele ( @group ) {
				unless (( $allele =~ /N/ ) || ( $allele =~ /Q/ )) {
					my $name = $element[0] . $allele;
					if ($name =~ /($gene\*\d+:\d+)/ ) {
						my $two = $1;
						if ( scalar @element > 2 ) {
							$group{ $two } = $gene . "*" .  $element[2];
						}
						else {
							$group{ $two } = "";
						}
					}
				}
			}
		}
	}
	return $group_ref;
}


1;
