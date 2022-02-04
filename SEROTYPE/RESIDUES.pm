#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: RESIDUES.pm 
# This module was developed to print key residues
# last modified and documented on February 2 2022

package RESIDUES;
use strict;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character
my $header = "Serotype,WHOAccepted,Allele,CIWD3.0,CWD2.0,EURCWD";

sub pattern {
	my ( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $basetype_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref ) = @_;
	open(FILE, ">output/target_" . $gene . "_" . $date . ".csv");
	print FILE $header . ",";
	# print residues
	for (my $index = 0; $index < scalar @$residues_ref; $index++) {
#		print $residues_ref->[ $index ] . "\n";
		print FILE $residues_ref->[ $index ];
		my $limit = scalar @$residues_ref - 1;
		if ( $index < $limit ) {
			print FILE ",";
		}
		else {
			print FILE "\n";
		}
	}
	
	my %org_fasta;
	foreach my $base ( @$basetype_ref ) {
		foreach my $type ( sort keys %$base_ref ) {
			if ( $base_ref->{ $type } eq $base ) {
#				print $type . "\n";
				print FILE $type . "," . $base . ",";

#				my $target = "";	# define target
				foreach my $head ( keys %$fasta_ref ) {	# go through fasta
					if ( $head =~ /$ref_ref->{ $type }/ ) {		# check accession number
						if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {	# allele name
							$org_fasta{ $type } =  $fasta_ref->{ $head };
							my $allele = $1;
							
							my $twoField = "";
							if ( $allele =~ /($gene\*\d+:\d+)/ ) {
								$twoField = $1;
							}
							if ( exists $ciwd_ref->{ $twoField } ) {
								print FILE $allele . "," . $ciwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE $allele . ",,";
							}

							if ( exists $cwd_ref->{ $twoField } ) {
								print FILE $cwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE ",";
							}
							if ( exists $ecwd_ref->{ $twoField } ) {
								print FILE $ecwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE ",";
							}
						}
						my $elements = scalar @$residues_ref;
						for ( my $index = 0; $index < $elements; $index++ ) {
							my $position = $residues_ref->[ $index ] + $leader;
					
							my $seq = "";
							$seq =  $fasta_ref->{ $head };
							unless ( $seq =~ /^M[A-Z]+/ ) {
								$seq = $partial_ref->{ $type } . $seq;
							}

							print FILE substr($seq, $position, 1);
							if ( $index != $elements - 1 ) {
								print FILE ",";
							}
						}
					}
				}
				print FILE "\n";
			}
		}
	}
	close FILE;

}

sub LAX {
	my ( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $basetype_ref, $base_ref, $group_ref, $ciwd_ref, $cwd_ref, $ecwd_ref ) = @_;
	open(FILE, ">output/target_LAX_" . $gene . "_" . $date . ".csv");
	print FILE $header . ",";
	# print residues
	for (my $index = 0; $index < scalar @$residues_ref; $index++) {
		print FILE $residues_ref->[ $index ];
		my $limit = scalar @$residues_ref - 1;
		if ( $index < $limit ) {
			print FILE ",";
		}
		else {
			print FILE "\n";
		}
	}
	

	my %ref;
	foreach my $base ( @$basetype_ref ) {
		foreach my $type ( sort keys %$base_ref ) {
			if ( $base_ref->{ $type } eq $base ) {
				print FILE $type . "," . $base . ",";
				my $target = "";	# define target
				my $lax_res_ref;
				# LAX residues
				if ( $gene eq "A" ) {
					$lax_res_ref = HLAA_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "B" ) {
					$lax_res_ref = HLAB_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "C" ) {
					$lax_res_ref = HLAC_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB1" ) {
					$lax_res_ref = DRB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DQB1" ) {
					$lax_res_ref = DQB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DQA1" ) {
					$lax_res_ref = DQA1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DPB1" ) {
					$lax_res_ref = DPB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DPA1" ) {
					$lax_res_ref = DPA1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				foreach my $head ( keys %$fasta_ref ) {	# go through fasta
					if ( $head =~ /$ref_ref->{ $type }/ ) {
						if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {
							my $allele = $1;

							my $twoField = "";
							if ( $allele =~ /($gene\*\d+:\d+)/ ) {
								$twoField = $1;
							}
							if ( exists $ciwd_ref->{ $twoField } ) {
								print FILE $allele . "," . $ciwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE $allele . ",,";
							}
							if ( exists $cwd_ref->{ $twoField } ) {
								print FILE $cwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE ",";
							}
							if ( exists $ecwd_ref->{ $twoField } ) {
								print FILE $ecwd_ref->{ $twoField } . ",";
							}
							else {
								print FILE ",";
							}
							$ref{ $type } = $allele;
						}
						my $elements = scalar @$lax_res_ref;

						for (my $index = 0; $index < scalar @$residues_ref; $index++) {
							my $position = $residues_ref->[ $index ] + $leader;

							my $test = 0;
							foreach my $lax_num ( @$lax_res_ref ) {
								if ( $lax_num ==  $residues_ref->[ $index ] ) {
									$test = 1;
								}
							}

							my $str_num = scalar @$residues_ref;
							if ( $test == 1 ) {

								my $seq = "";
								$seq =  $fasta_ref->{ $head };
								unless ( $seq =~ /^M[A-Z]+/ ) {
									$seq = $partial_ref->{ $type } . $seq;
								}

								print FILE substr($seq, $position, 1);
								if ( $index != $str_num - 1 ) {
									print FILE ",";
								}
							}
							else {
								if ( $index != $str_num - 1 ) {
									print FILE ",";
								}
							}
						}
					}
				}
				print FILE "\n";
			}
		}
	}
	close FILE;

	open(FILE, ">output/target_position_LAX_" . $gene . "_" . $date . ".csv");
	print FILE $header . "\n";
	foreach my $base ( @$basetype_ref ) {
		foreach my $type ( sort keys %$base_ref ) {
			if ( $base_ref->{ $type } eq $base ) {
				
				my $twoField = "";
				if ( $ref{ $type } =~ /($gene\*\d+:\d+)/ ) {
					$twoField = $1;
				}
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE $type . "," . $base . "," . $ref{ $type } . "," . $ciwd_ref->{ $twoField } . ",";
				}
				else {
					print FILE $type . "," . $base . "," . $ref{ $type } . ",,";
				}
				if ( exists $cwd_ref->{ $twoField } ) {
					print FILE $cwd_ref->{ $twoField } . ",";
				}
				else {
					print FILE ",";
				}
				if ( exists $ecwd_ref->{ $twoField } ) {
					print FILE $ecwd_ref->{ $twoField } . ",";
				}
				else {
					print FILE ",";
				}

				my $lax_res_ref;
				if ( $gene eq "A" ) {
					$lax_res_ref = HLAA_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "B" ) {
					$lax_res_ref = HLAB_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "C" ) {
					$lax_res_ref = HLAC_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB1" ) {
					$lax_res_ref = DRB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DQB1" ) {
					$lax_res_ref = DQB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DQA1" ) {
					$lax_res_ref = DQA1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DPB1" ) {
					$lax_res_ref = DPB1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DPA1" ) {
					$lax_res_ref = DPA1_INFO::RESIDUES( $group_ref->{ $type } );
				}
				my $limit = scalar @$lax_res_ref;
				for ( my $index = 0; $index < $limit; $index++ ) {
					print FILE $lax_res_ref->[ $index ];
					if ( $index < $limit -1 ) {
						print FILE ",";
					}
				}
				print FILE "\n";
			}
		}
	}
	close FILE;
}


1;
