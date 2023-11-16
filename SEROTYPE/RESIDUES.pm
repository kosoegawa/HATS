#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: RESIDUES.pm 
# This module was developed to print key residues
# last modified and documented on November 7 2023

package RESIDUES;
use strict;
use ASSIGNED_SHORT;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character
my $header = "Serotype,WHOAccepted,Allele,CIWD3.0,CWD2.0,EURCWD";

sub pattern {
	my ( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $basetype_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref ) = @_;
	open(FILE, ">output/target_" . $gene . "_" . $date . ".csv");
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
	my %elements_full;
	my $elements_full_ref = \%elements_full;
	my %org_fasta;
	foreach my $base ( @$basetype_ref ) {
		foreach my $type ( sort keys %$base_ref ) {
			my $line = "FULL,";
			if ( $base_ref->{ $type } eq $base ) {
				print FILE $type . "," . $base . ",";

				foreach my $head ( keys %$fasta_ref ) {	# go through fasta
					my $allele = "";
					my $twoField = "";
					if ( $head =~ /$ref_ref->{ $type }/ ) {		# check accession number
						if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {	# allele name
							$org_fasta{ $type } =  $fasta_ref->{ $head };
							$allele = $1;
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
							$line = $line . substr($seq, $position, 1);
							if ( $index != $elements - 1 ) {
								print FILE ",";
								$line = $line . ",";
							}
						}
						$line = $line . "\n";
					}
				}
				print FILE "\n";
				$elements_full{ $type } = $line;
			}
		}
	}
	close FILE;
	return $elements_full_ref;
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
				elsif ( $gene eq "DRB3" ) {
					$lax_res_ref = DRB3_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB4" ) {
					$lax_res_ref = DRB4_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB5" ) {
					$lax_res_ref = DRB5_INFO::RESIDUES( $group_ref->{ $type } );
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
				elsif ( $gene eq "DRB3" ) {
					$lax_res_ref = DRB3_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB4" ) {
					$lax_res_ref = DRB4_INFO::RESIDUES( $group_ref->{ $type } );
				}
				elsif ( $gene eq "DRB5" ) {
					$lax_res_ref = DRB5_INFO::RESIDUES( $group_ref->{ $type } );
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

sub ELEMENTS {
	my ( $elements_full_ref,$fasta_ref,$gene,$null_ref,$qallele_ref,$residues_ref,$leader,$partial_ref,$assigned_ref,$short_ref) = @_;
	
	my %antigen;
	my $antigen_ref = \%antigen;
	$antigen_ref = ASSIGNED_SHORT::ANTIGEN();
	my %elements2;
	my $elements2_ref = \%elements2;
	foreach my $head ( sort keys %$fasta_ref ) {	# go through fasta
		my $allele = "";
		my $twoField = "";
		my $line = "";
		my $type = "";
		if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {	# allele name
			$allele = $1;
			if (( $allele eq "DRB1*11:13:01" ) || ( $allele eq "DQB1*06:18:01" ) || ( $allele eq "DPB1*26:01:01" )) {	# DRB1*11:13:01 is partial, but DRB1*11:13:02 is full length
				next;
			}
			unless (( exists $null_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
				if ( $allele =~ /($gene\*\d+:\d+)/ ) {
					$twoField = $1;
				}
				if ( exists $elements2_ref->{ $twoField } ) {
					next;
				}
				else {
					my $elements = scalar @$residues_ref;
					for ( my $index = 0; $index < $elements; $index++ ) {
						my $position = $residues_ref->[ $index ] + $leader;
					
						my $seq = "";
						$seq =  $fasta_ref->{ $head };
						unless ( $seq =~ /^M[A-Z]+/ ) {
							if ( exists $partial_ref->{ $allele } ) {
								$seq = $partial_ref->{ $allele } . $seq;
							}
							else {
								$seq = $partial_ref->{ "general" } . $seq;
							}
							my $end_position = $residues_ref->[ $elements - 1 ] + $leader;
							if ( length( $seq ) < $end_position ) {
								$seq = $seq . $partial_ref->{ "tail" }; 
								#print $allele . "," . $seq . "\n";
							} 
						}
		
						$line = $line . substr($seq, $position, 1);
						if ( $index != $elements - 1 ) {
							$line = $line . ",";
						}
					}

					if ( exists $assigned_ref->{ $allele } ) {
						if ( exists $elements_full_ref->{ $assigned_ref->{ $allele } } ) {	# elements_ref contains FULL serotype
							$elements2_ref->{ $twoField }->[0] = $assigned_ref->{ $allele };	# $type
							$elements2_ref->{ $twoField }->[1] = $elements_full_ref->{ $assigned_ref->{ $allele } };	#residues
						}
						elsif ( $assigned_ref->{ $allele } =~ /(\S+-*\d+)_(LAX)/ ) {	# LAX
							$type = $1;
							my $lax = $2;
							$elements2_ref->{ $twoField }->[0] = $type;	# $type
							$line = "SEROTYPE," . $line;
							$elements2_ref->{ $twoField }->[1] = $line . "\n";
						}
					}
					else {	# INCOMPLETE or UNASSIGNED
						if ( exists $short_ref->{ $allele } ) {
							my $num = scalar @{$short_ref->{ $allele }};
							if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
								my $residue = 0;
								my $sero = $1 . "-" . $2;
								$sero =~ s/DRB1/DR/;
								$sero =~ s/DQB1/DQ/;
								my $group = "";
								my $test = 0;
								foreach my $short ( sort @{$short_ref->{ $allele }} ) {
									if ( $short =~ /(\S+)_(\d+)/ ) {
										$group = $1;
										$residue = $2;
										$type = $group;
										if ( $group =~ /$sero/ ) {		# allele name and sero type matches
											$test = 1;
											if ( exists $antigen{ $group } ) {
												$type = $antigen{ $group };
											}
											last;
										}
									}
								}
								if ( $test != 1 ) {
									my @tmp = sort @{$short_ref->{ $allele }};
									my $first_name = $tmp[0];
									if ( $first_name =~ /(\S+)_(\d+)/ ) {
										$group = $1;
										$residue = $2;
										$type = $group;
										if ( exists $antigen{ $group } ) {
											$type = $antigen{ $group };
										}
									}
								}
								$elements2_ref->{ $twoField }->[0] = $type;	# $type
								$line = "INCOMPLETE," . $line;
								$elements2_ref->{ $twoField }->[1] = $line . "\n";
							}
						}
						else {
							if ( $allele eq "DRB1*04:20" ) {
								$type = "DR-0403";
								$line = "InSilico," . $line;
							}
							else {
								$type = "UNA";
								$line = "UNASSIGNED," . $line;
							}
							$elements2_ref->{ $twoField }->[0] = $type;	# $type
							$elements2_ref->{ $twoField }->[1] = $line . "\n";
						}
					}
				}
			}
		}
	}
	return $elements2_ref;
}

# This is extra to combine all residue positions for A, B and C
sub RESIDUES_ABC {
	my ( $residues_A_ref, $residues_B_ref, $residues_C_ref ) = @_;
	my @combined = ();
	push @combined, @$residues_A_ref; 
	push @combined, @$residues_B_ref; 
	push @combined, @$residues_C_ref; 

	my %seen;
	my @unique;
	my $unique_ref = \@unique;
	foreach my $value ( sort { $a <=> $b } @combined ) {
		unless ( exists $seen{ $value } ) {
			push @unique, $value;
			$seen{ $value } = 0;
		}
	}
	
	return $unique_ref;
}

# This is extra to combine all residue positions for A, B and C
sub RESIDUES_DR {
	my ( $residues_DRB1_ref, $residues_DRB3_ref, $residues_DRB4_ref,$residues_DRB5_ref ) = @_;
	my @combined = ();
	push @combined, @$residues_DRB1_ref; 
	push @combined, @$residues_DRB3_ref; 
	push @combined, @$residues_DRB4_ref; 
	push @combined, @$residues_DRB5_ref; 

	my %seen;
	my @unique;
	my $unique_ref = \@unique;
	foreach my $value ( sort { $a <=> $b } @combined ) {
		unless ( exists $seen{ $value } ) {
			push @unique, $value;
			$seen{ $value } = 0;
		}
	}
	
	return $unique_ref;
}

1;
