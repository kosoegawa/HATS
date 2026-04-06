#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: RESIDUES.pm 
# This module was developed to print key residues
# last reviewed on April 5 2026

package RESIDUES;
use strict;
use COMBINE;
use POSIX qw(strftime);
use READ_MSF;

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character
my $header = "Antigen,Protein,CIWD3.0,CWD2.0,EURCWD";

sub pattern {
	my ( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref ) = @_;
	open(FILE, ">output/target_" . $gene . "_" . $date . ".csv");
	print FILE $header . ",";

	my $antigen_ref = COMBINE::ANTIGEN();

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

	foreach my $type ( sort keys %$base_ref ) {	#type is serotype	
#		my $line = "FULL,";
		my $line = "F,";
		if ( exists $antigen_ref->{ $type } ) {
			print FILE $antigen_ref->{ $type } . ",";
		}
		else {
			print FILE $type . ",";
		}

		foreach my $head ( keys %$fasta_ref ) {	# go through fasta
			my $allele = "";
			my $twoField = "";
			if ( $head =~ /$ref_ref->{ $type }/ ) {		# check accession number
				if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {	# allele name
					$org_fasta{ $type } =  $fasta_ref->{ $head };
					$allele = $1;
					$twoField = $allele;
					if ( $allele =~ /($gene\*\d+:\d+)/ ) {
						$twoField = $1;
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE $twoField . "," . $ciwd_ref->{ $twoField } . ",";
					}
					else {
						print FILE $twoField . ",,";
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

					my $aa = substr($seq, $position, 1);		# AA residue at target position
					if (( $gene eq "B" ) && ( $residues_ref->[ $index ] == 67 )) {
						if (( $aa eq "S" ) || ( $aa eq "C" )) {
							$aa= "S/C";
						}
						elsif (( $aa eq "Y" ) || ( $aa eq "F" )) {
							$aa= "Y/F";
						}
					}
					elsif (( $gene eq "B" ) && ( $residues_ref->[ $index ] == 167) && (( $aa eq "S" ) || ( $aa eq "G" ))) {
						$aa = "S/G";
					}
					elsif (( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 67 ) && (( $aa eq "I" ) || ( $aa eq "L" ))) {
						$aa= "I/L";
					}
					elsif (((( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 71 )) || (( $gene eq "DPB1" ) && ( $residues_ref->[ $index ] == 69 ))) && (( $aa eq "K" ) || ( $aa eq "R" ))) {
						$aa= "R/K";
					}
					elsif (( $gene =~ /DRB/ ) && (!( $type eq "DR0301")) && (!( $type eq "DR0302" )) && ( $residues_ref->[ $index ] == 47 ) && (( $aa eq "F" ) || ( $aa eq "Y" ))) {
						$aa= "F/Y";
					}
					elsif (($gene eq "DQB1") && ($residues_ref->[ $index ] == 57) && (!( $type eq "DQ0302")) && (!( $type eq "DQ0303" ))) {
						if ((( $aa eq "D" ) || ( $aa eq "V" ) || ( $aa eq "S" )) || ( $aa eq "A" )) {
							$aa = "D/V/S/A";
						}
					}
					elsif (( $gene eq "DPB1" ) && ( $residues_ref->[ $index ] == 84 ) && (( $aa eq "G" ) || ( $aa eq "V" ))) {
						$aa= "G/V";
					}

					print FILE $aa;
					$line = $line . substr($seq, $position, 1);	# needs to fix later
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

	close FILE;
	return $elements_full_ref;
}

sub LAX {
	my ( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $base_ref, $group_ref, $ciwd_ref, $cwd_ref, $ecwd_ref ) = @_;
	open(FILE, ">output/target_LAX_" . $gene . "_" . $date . ".csv");
	print FILE $header . ",";

	my $antigen_ref = COMBINE::ANTIGEN();

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

	foreach my $type ( sort keys %$base_ref ) {
		if ( exists $antigen_ref->{ $type } ) {
			print FILE $antigen_ref->{ $type } . ",";
		}
		else {
			print FILE $type . ",";
		}
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

					my $twoField = $allele;
					if ( $allele =~ /($gene\*\d+:\d+)/ ) {
						$twoField = $1;
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE $twoField . "," . $ciwd_ref->{ $twoField } . ",";
					}
					else {
						print FILE $twoField . ",,";
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
							$aa = "S/G";
						}
						elsif (( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 67 ) && (( $aa eq "I" ) || ( $aa eq "L" ))) {
							$aa= "I/L";
						}
						elsif (( $gene =~ /DRB/ ) && ( $residues_ref->[ $index ] == 71 ) && (( $aa eq "K" ) || ( $aa eq "R" ))) {
							$aa= "R/K";
						}
						elsif (( $gene =~ /DRB/ ) && (!( $type eq "DR0301")) && (!( $type eq "DR0302" )) && ( $residues_ref->[ $index ] == 47 ) && (( $aa eq "F" ) || ( $aa eq "Y" ))) {
							$aa= "F/Y";
						}
						elsif (($gene eq "DQB1") && ($residues_ref->[ $index ] == 57) && (!( $type eq "DQ0302")) && (!( $type eq "DQ0303" ))) {
							if ((( $aa eq "D" ) || ( $aa eq "V" ) || ( $aa eq "S" )) || ( $aa eq "A" )) {
								$aa = "D/V/S/A";
							}
						}
						elsif (( $gene eq "DPB1" ) && ( $residues_ref->[ $index ] == 84 ) && (( $aa eq "G" ) || ( $aa eq "V" ))) {
							$aa= "G/V";
						}

						print FILE $aa;
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
	close FILE;

	open(FILE, ">output/target_position_LAX_" . $gene . "_" . $date . ".csv");
	print FILE $header . "\n";

	foreach my $type ( sort keys %$base_ref ) {
		if ( exists $antigen_ref->{ $type } ) {
			print FILE $antigen_ref->{ $type } . ",";
		}
		else {
			print FILE $type . ",";
		}

		my $twoField = $ref{ $type };
		if ( $ref{ $type } =~ /($gene\*\d+:\d+)/ ) {
			$twoField = $1;
		}
		if ( exists $ciwd_ref->{ $twoField } ) {
			print FILE $twoField . "," . $ciwd_ref->{ $twoField } . ",";
		}
		else {
			print FILE $twoField . ",,";
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
	close FILE;
}


# generate DEPs for each protein
sub ELEMENTS {
	my ( $elements_full_ref,$fasta_ref,$gene,$null_ref,$qallele_ref,$residues_ref,$leader,$partial_ref,$assigned_ref,$short_ref,$ref_allele,$msf_ref ) = @_;
	
	my $antigen_ref = COMBINE::ANTIGEN();
	my $protein_ref = STRASSIGN::PROTEIN( $assigned_ref ); 
 	my $protein_allele_ref = STRASSIGN::PROTEIN_ALLELE_REF( $assigned_ref );
	my $missing_seq_ref = READ_MSF::MISSING_SEQ( $gene,$null_ref,$qallele_ref,$ref_allele,$msf_ref );

	my %elements2;
	my $elements2_ref = \%elements2;
	foreach my $head ( sort keys %$fasta_ref ) {	# go through fasta
		my $allele = "";
		my $twoField = "";
		my $line = "";
		my $type = "";
		if ( $head =~ /HLA:\S+ ($gene\*\d+:\d+:*\d*:*\d*\S*) \d+ bp/ ) {	# allele name
			$allele = $1;
			unless (( exists $null_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
				if ( $allele =~ /($gene\*\d+:\d+)/ ) {
					$twoField = $1;
					unless (( exists $assigned_ref->{ $allele } )) {	# DPB1*482:01:01 does nor exist
						if ( exists $protein_ref->{ $twoField } ) {	# resolve exon 3 missing sequence issue
							next;
						}
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
								$seq = $missing_seq_ref->{ $allele } . $seq;

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
					}

					if ( exists $assigned_ref->{ $allele } ) {
						my $serotype =  $assigned_ref->{ $allele };

						if ( $protein_allele_ref->{ $twoField } = $allele ) {	#assigned two filed exists
							if ( exists $elements_full_ref->{ $serotype } ) {	# elements_ref contains FULL serotype
								if ( exists $antigen_ref->{ $serotype } ) {
									$serotype = $antigen_ref->{ $serotype };
								}

								$elements2_ref->{ $twoField }->[0] = $serotype;	# $type
#								$elements2_ref->{ $twoField }->[1] = $elements_full_ref->{ $assigned_ref->{ $allele } };	#residues
								$line = "F," . $line;
								$elements2_ref->{ $twoField }->[1] = $line . "\n";
							}
							elsif ( $serotype =~ /(\S+-*\d+)_(LAX)/ ) {	# LAX
								$type = $1;
								my $lax = $2;
								$elements2_ref->{ $twoField }->[0] = $type;	# $type
								$line = "S," . $line;
								$elements2_ref->{ $twoField }->[1] = $line . "\n";
							}
						}
					}
					else {	# INCOMPLETE or UNASSIGNED
						if ( exists $short_ref->{ $allele } ) {	# INCOMPLETE
							my $num = scalar @{$short_ref->{ $allele }};
#							my $serotype = "";
							my $group = "";
							my $residue = 0;
	
							if ( $num == 1 ) {	#short
								my $short = $short_ref->{ $allele }->[0];
								if ( $short =~ /(\S+)_(\d+)/ ) {
									$group = $1;
									$residue = $2;
									$type = $group;
									if ( exists $antigen_ref->{ $group } ) {
										$type = $antigen_ref->{ $group };
									}
								}
							}
							else {			#short cross-reactive
								if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
									my $test = 0;
									my $sero = $1 . $2;	# need to modify here, problem for DPB
									$sero =~ s/DRB1/DR/;
									$sero =~ s/C/Cw/;	#added here 2/20/26 to fix bug
									$sero =~ s/DPB1/DPB/;
								
									my @sorted_short = sort @{$short_ref->{ $allele }};
									my $short = $sorted_short[0];
									if ( $short =~ /(\S+)_(\d+)/ ) {
										$group = $1;
										$type = $group;
										if ( exists $antigen_ref->{ $group } ) {
											$type = $antigen_ref->{ $group };
										}
										$residue = $2;
									}

									if ( exists $protein_ref->{ $twoField } ) {	# resolve exon 3 missing sequence issue
										$type = $protein_ref->{ $twoField };
										$group = $type;
									}

									unless ( $sero =~ /DPB/ ) {
										foreach my $short ( sort @{$short_ref->{ $allele }} ) {
											if ( $short =~ /(\S+)_(\d+)/ ) {
												$group = $1;		# problematic
												$residue = $2;
												if ( $group =~ /$sero/ ) {		# allele name and sero type matches, this is important
													if ( exists $antigen_ref->{ $group } ) {
														$type = $antigen_ref->{ $group };
													}
													else {
														$type = $group;
													}

													$test = 1;
													last;
												}
											}
										}
										if ( $test == 0 ) {
											my @sorted_short = sort @{$short_ref->{ $allele }};
											my $short = $sorted_short[0];
											if ( $short =~ /(\S+)_(\d+)/ ) {
												$group = $1;
												$residue = $2;
											}
										}
									}
								}
							}
							$elements2_ref->{ $twoField }->[0] = $type;	# $type
							$line = "I," . $line;
							$elements2_ref->{ $twoField }->[1] = $line . "\n";
						}
						else {
							if ( $allele eq "DRB1*04:20" ) {
								$type = "DR0403";
								$line = "InSilico," . $line;
							}
							else {
								$type = "UNA";
								$line = "U," . $line;
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



1;
