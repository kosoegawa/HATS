#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: ASSIGNED_SHORT.pm 
# This module was developed to print table
# last reviewed, modified and documented on December 13 2023

package ASSIGNED_SHORT;
use strict;
use GROUP_SORT;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

my %antigen;
my $antigen_ref = \%antigen;
$antigen{ "A-0101" } = "A1";
$antigen{ "A-0203" } = "A203";
$antigen{ "A-0210" } = "A210";
$antigen{ "A-2403" } = "A2403";
$antigen{ "A-2501" } = "A25";
$antigen{ "A-1101" } = "A11";
$antigen{ "A-2901" } = "A29";
$antigen{ "A-7401" } = "A74";
$antigen{ "A-6901" } = "A69";
$antigen{ "A-3601" } = "A36";
$antigen{ "A-4301" } = "A43";
$antigen{ "A-8001" } = "A80";

$antigen{ "B-5102" } = "B5102";
$antigen{ "B-5103" } = "B5103";
$antigen{ "B-5201" } = "B52";
$antigen{ "B-0703" } = "B703";
$antigen{ "B-4501" } = "B45";
$antigen{ "B-1302" } = "B13";
$antigen{ "B-1401" } = "B64";
$antigen{ "B-1402" } = "B65";
$antigen{ "B-1512" } = "B76";
$antigen{ "B-1514" } = "B76";
$antigen{ "B-1513" } = "B77";
$antigen{ "B-3901" } = "B3901";
$antigen{ "B-3902" } = "B3902";
$antigen{ "B-5701" } = "B57";
$antigen{ "B-5801" } = "B58";
$antigen{ "B-4901" } = "B49";
$antigen{ "B-5001" } = "B50";
$antigen{ "B-5401" } = "B54";
$antigen{ "B-5501" } = "B55";
$antigen{ "B-5601" } = "B56";
$antigen{ "B-2708" } = "B2708";
$antigen{ "B-4005" } = "B4005";
$antigen{ "B-4101" } = "B41";
$antigen{ "B-4201" } = "B42";
$antigen{ "B-4601" } = "B46";
$antigen{ "B-5301" } = "B53";
$antigen{ "B-5901" } = "B59";
$antigen{ "B-6701" } = "B67";
$antigen{ "B-7301" } = "B73";
$antigen{ "B-7801" } = "B78";
$antigen{ "B-8101" } = "B81";
$antigen{ "B-8201" } = "B82";

$antigen{ "C-0303" } = "Cw9";
$antigen{ "C-0304" } = "Cw10";

$antigen{ "DR-0103" } = "DR103";
$antigen{ "DR-1505" } = "DR-1501";
$antigen{ "DR-1605" } = "DR-1602";
$antigen{ "DR-0301" } = "DR17";
$antigen{ "DR-0302" } = "DR18";
$antigen{ "DR-0701" } = "DR7";
$antigen{ "DR-0901" } = "DR9";
$antigen{ "DR-1001" } = "DR10";
$antigen{ "DR-0404" } = "DR-0401";
$antigen{ "DR-1312" } = "DR-1303";	#added
$antigen{ "DR-1320" } = "DR-1301";	#added
$antigen{ "DR-1349" } = "DR-1305";	#added
$antigen{ "DR-1417" } = "DR-1402";
$antigen{ "DR-1419" } = "DR-1402";	#added

$antigen{ "DQ-0201" } = "DQ2";
$antigen{ "DQ-0203" } = "DQ2";
$antigen{ "DQ-0301" } = "DQ7";
$antigen{ "DQ-0304" } = "DQ7";
$antigen{ "DQ-0302" } = "DQ8";
$antigen{ "DQ-0303" } = "DQ9";
$antigen{ "DQ-0401" } = "DQ4";
$antigen{ "DQ-0501" } = "DQ5";
$antigen{ "DQ-0502" } = "DQ5";
$antigen{ "DQ-0503" } = "DQ5";
$antigen{ "DQ-0602" } = "DQ6";	#comment out for SAB
$antigen{ "DQ-0608" } = "DQ6";
$antigen{ "DQ-0610" } = "DQ6";

$antigen{"DP-0101"} = "DP-01";
#$antigen{"DP-0201"} = "DP-0201";
$antigen{"DP-0301"} = "DP-03";
#$antigen{"DP-0401"} = "DP-0401";
$antigen{"DP-1001"} = "DP-10";
$antigen{"DP-1501"} = "DP-15";
$antigen{"DP-1801"} = "DP-18";
$antigen{"DP-4601"} = "DP-46";
#$antigen{"DP-0402"} = "DP-0402";
#$antigen{"DP-0202"} = "DP-0202";
$antigen{"DP-0601"} = "DP-06";
$antigen{"DP-1301"} = "DP-13";
$antigen{"DP-4501"} = "DP-45";
$antigen{"DP-8001"} = "DP-80";

$antigen{"DP-1101"} = "DP-01"; 
$antigen{"DP-3401"} = "DP-15"; 
$antigen{"DP-13601"} = "DP-45";  
$antigen{"DP-6901"} = "DP-03";

my @broad = ("A9","A10","A19","A28","B5","B12","B14","B15","B16","B17","B21","B22","B40","B70","Cw3","DR2","DR3","DR5","DR6","DQ1","DQ3");

# print SHORT
sub PRINT {
	my ( $unassigned_ref, $short_ref ) = @_;
	print "SHORT\n";
	open(FILE, ">output/Assigned_SHORT_" . $date . ".csv");
	foreach my $allele (sort @$unassigned_ref) {
		if ( exists $short_ref->{ $allele } ) {
			print FILE $allele . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			for (my $index = 0; $index < $num; $index++) {
				print FILE $short_ref->{ $allele }->[$index];
				if ( $index < $num - 1) {
					print FILE ",";
				}
				else {
					print FILE "\n";
				}
			}
		}
		else {
			print FILE $allele . "\n";
		}
	}
	close FILE;
}

sub PRINT_RESIDUES {
	my ( $elements_ref,$gene,$residues_all_ref,$database ) = @_;
	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$elements_ref;
	my $elements_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically
	my %residue; 
	open(FILE, ">output/" . $gene . "_Allele_Residues_" . $database . "_" . $date . ".csv");
	print FILE "Allele,Serotype,COMMENT,";
	# print residues
	for (my $index = 0; $index < scalar @$residues_all_ref; $index++) {
		print FILE $residues_all_ref->[ $index ];
		my $limit = scalar @$residues_all_ref - 1;
		if ( $index < $limit ) {
			print FILE ",";
		}
		else {
			print FILE "\n";
		}
	}
	
	foreach my $twoField ( @$elements_sorted_ref ) {		#go through all alleles
		unless ( exists $residue{ $twoField } ) {
			if ( exists $antigen{ $elements_ref->{ $twoField }->[0] } ) {
				print FILE $twoField . ",". $antigen{ $elements_ref->{ $twoField }->[0] } . "," . $elements_ref->{ $twoField }->[1];	# paragrarh mark is included in elements
			}
			else {
				print FILE $twoField . ",". $elements_ref->{ $twoField }->[0] . "," . $elements_ref->{ $twoField }->[1];	# paragrarh mark is included in elements
			}
			$residue{ $twoField } = 0;
		}
	}
	close FILE;
}

# print combined table
sub COMBINED {
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2,$dr_ref) = @_;
	print "COMBINED\n";
	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically

	open(FILE, ">output/" . $gene . "_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
	print FILE "Allele,COMMENT,Serotype,ImputedAntigen,Broad,CIWD3.0,CWD2.0,EURCWD,Bw46C12DR5X\n";	# modified to ImputedAntigen

	foreach my $allele ( @$alleles_sorted_ref ) {		#go through all alleles
		my $twoField = "";
		unless (( exists $nullAllele_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
			if ( $allele =~ /($gene\*\d+:\d+)/ ) {
				$twoField = $1;
			}
		}
		if ( exists $nullAllele_ref->{ $allele } ) {		# Null
			print FILE $allele . ",," . $nullAllele_ref->{ $allele } . "\n";
		}
		elsif ( exists $qallele_ref->{ $allele } ) {		# Questionable
			print FILE $allele . ",," . $qallele_ref->{ $allele } . "\n";
		}
		elsif ( exists $assigned_ref->{ $allele } ) {	# assigned
			if ( $assigned_ref->{ $allele } =~ /(\S+_*\S*)_(LAX)/ ) {	# LAX
				my $group = $1;
				my $lax = $2;
				my $serotype = $group;
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				if ( $lax eq "LAX" ) {
					$lax = "SEROTYPE";
				}
				if ( exists $cross_ref->{ $allele } ) {		# cross-reactivity
					print FILE $allele . "," . $lax . "_C," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group };
					my $cross = "";
					for ( my $index = 0; $index < scalar @{$cross_ref->{ $allele }}; $index++ ) {
						if ( $index == 0 ) {
							$cross = $cross_ref->{ $allele }->[$index];
						}
						else {
							$cross = $cross . "," . $cross_ref->{ $allele }->[$index];
						}
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE "," . $ciwd_ref->{ $twoField};
					}
					else {
						print FILE ",";
					}
					if ( exists $cwd_ref->{ $twoField } ) {
						print FILE "," . $cwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $ecwd_ref->{ $twoField } ) {
						print FILE "," . $ecwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if (( exists $bw_ref->{ $group } ) && ( !exists $bw_ref2->{ $allele } )) {
						print FILE "," . $bw_ref->{ $group };
					}
					elsif ( exists $bw_ref2->{ $allele }) {
						print FILE "," . $bw_ref2->{ $allele } . ",";
					}
					else {
						print FILE ",";
					}
					print FILE "," . $cross
				}
				else {	# no cross-reactive
					if (( $serotype eq "A-0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
						print FILE $allele . "," . $lax . "," . $serotype . ",A11,A11";
					}
					elsif (( $serotype eq "B-4406" ) && ( $allele =~ /B\*51/ )) {	
						print FILE $allele . "," . $lax . "," . $serotype . ",B51,B5";
					}
					else {
						print FILE $allele . "," . $lax . "," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group };
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE "," . $ciwd_ref->{ $twoField};
					}
					else {
						print FILE ",";
					}
					if ( exists $cwd_ref->{ $twoField } ) {
						print FILE "," . $cwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $ecwd_ref->{ $twoField } ) {
						print FILE "," . $ecwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $bw_ref2->{ $allele }) {
						print FILE "," . $bw_ref2->{ $allele };
					}
					elsif ( exists $bw_ref->{ $group }) {
						print FILE "," . $bw_ref->{ $group };
					}
					else {
						print FILE ",";
					}
				}
				print FILE "\n";
			}
			else {		# stringent
				my $group = $assigned_ref->{ $allele };
				my $serotype = $group;
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				my $full = ",FULL,";
				if (( $serotype eq "A-0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
					print FILE $allele . $full . $serotype . ",A11,A11";
				}
				elsif (( $serotype eq "B-4406" ) && ( $allele =~ /B\*51/ )) {	
					print FILE $allele . $full . $serotype . ",B51,B5";
				}
				else {
					print FILE $allele . $full  . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group };
				}
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField};
				}
				else {
					print FILE ",";
				}
				if ( exists $cwd_ref->{ $twoField } ) {
					print FILE "," . $cwd_ref->{ $twoField };
				}
				else {
					print FILE ",";
				}
				if ( exists $ecwd_ref->{ $twoField } ) {
					print FILE "," . $ecwd_ref->{ $twoField };
				}
				else {
					print FILE ",";
				}
				if ( exists $bw_ref->{ $group }) {
					print FILE "," . $bw_ref->{ $group };
				}
				elsif ( exists $bw_ref2->{ $allele }) {
					print FILE "," . $bw_ref2->{ $allele };
				}
				else {
					print FILE ",";
				}
				print FILE "\n";
			}
		}
		elsif ( exists $short_ref->{ $allele } ) {	# short
			print FILE $allele . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
				my $residue = 0;
				my $sero = $1 . "-" . $2;
				$sero =~ s/DRB1/DR/;
				my $serotype = "";
				my $group = "";
				my $test = 0;
				foreach my $short ( sort @{$short_ref->{ $allele }} ) {
					if ( $short =~ /(\S+)_(\d+)/ ) {
						$group = $1;
						$serotype = $group;
						$residue = $2;
						if ( $group =~ /$sero/ ) {		# allele name and sero type matches
							$test = 1;
							last;
						}
					}
				}
				
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				# go through three different conditions: $base_type and $type is identical
				if ( $test == 1 ) {
					if ( $num == 1 ) {	# specific short
						print FILE "S," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE "SC," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
				}
				else {		# no match => capture the first element
					my @tmp = sort @{$short_ref->{ $allele }};
					my $first_name = $tmp[0];
					if ( $first_name =~ /(\S+)_(\d+)/ ) {
						$group = $1;
						$serotype = $group;
						if ( exists $antigen{ $group } ) {
							$serotype = $antigen{ $group };
						}

						$residue = $2;
					}
					if ( $num == 1 ) {	# specific short
						print FILE "S," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE "SC," . $serotype . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
				}

				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",";
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

				if ( $gene eq "B" ) {	# HLA-B
					if ( exists $bw_ref->{ $group } ) {
						if ( !exists $bw_ref2->{ $allele } ) {
							if ( $bw_ref->{ $group } eq "Bw4" ) {
								if (( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Bw6" ) {
								if (( $residue == 76 ) ||( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 76) && ( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Negative" ) {
								print FILE $bw_ref->{ $group } . ",";
							}
							else {
								print FILE "Negative,";
							}
						}
						else {
							print FILE $bw_ref2->{ $allele } . ",";
						}
					}
				}
				else {
					if ( exists $bw_ref->{ $group } ) {
						if ( exists $bw_ref2->{ $allele } ) {
							if ( $bw_ref->{ $group } eq $bw_ref2->{ $allele } ) {
								print FILE $bw_ref->{ $group } . ",";
							}
							else {
								print FILE $bw_ref2->{ $allele } . ",";
							}
						}
						elsif ( !exists $bw_ref2->{ $allele } ) {
							if ( $bw_ref->{ $group } eq "Bw4" ) {
								if (( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Bw6" ) {
								if (( $residue == 76 ) ||( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 76) && ( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Negative" ) {
								print FILE $bw_ref->{ $group } . ",";
							}

							elsif (($bw_ref->{$group} eq "C1") && ($residue != 80)) {
								print FILE $bw_ref->{ $group } . ",";
							}
							elsif (($bw_ref->{$group} eq "C2") && ($residue != 80)) {
								print FILE $bw_ref->{ $group } . ",";
							}
							else {
								print FILE ",";
							}
						}
					}
					elsif ( exists $bw_ref2->{ $allele }) {
						print FILE $bw_ref2->{ $allele } . ",";
					}
					else {
						print FILE ",";
					}
				}
			}
			
			my $index = 0;
			foreach my $short ( sort @{$short_ref->{ $allele }} ) {
				print FILE $short;
				if ( $index < $num - 1) {
					print FILE ",";
				}
				else {
					print FILE "\n";
				}
				$index++;
			}
		}
		else {		# no match
			print FILE $allele;
			if (( $allele =~ /DRB1\*08:04:02/ ) || ( $allele =~ /DRB1\*04:20/ ) ||
			( $allele =~ /DQB1\*05:03:02/ ) || ( $allele =~ /DQB1\*06:01:02/ ) || ( $allele =~ /DQB1\*06:05:02/ ) ||
			( $allele =~ /DQB1\*06:06/ )) {	# missing key residues 9 - 14
				if ( $allele =~ /DRB1\*08:04:02/ ) {
					print FILE ",InSilico,DR-0801,DR8,DR8";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",InSilico,DR-0403,DR4,DR4";
				}
				elsif ( $allele =~ /DQB1\*05:03:02/ ) {
					print FILE ",InSilico,DQ5,DQ5,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:01:02/ ) {
					print FILE ",InSilico,DQ6,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:05:02/ ) {
					print FILE ",InSilico,DQ6,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:06/ ) {
					print FILE ",InSilico,DQ6,DQ6,DQ1";
				}
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",,";
				}
				if ( exists $cwd_ref->{ $twoField } ) {
					print FILE $cwd_ref->{ $twoField } . ",";
				}
				else {
					print FILE ",";
				}
				if ( exists $ecwd_ref->{ $twoField } ) {
					print FILE $ecwd_ref->{ $twoField } . "\n";
				}
				else {
					print FILE "\n";
				}

			}
			else {
				print FILE ",UNA,U,None,None,None";
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",,";
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
				if ( exists $bw_ref2->{ $allele }) {
					print FILE $bw_ref2->{ $allele } . "\n";
				}
				else {
					if ( $gene eq "B" ) {
						print FILE "Negative\n";
					}
					else {
						print FILE "\n";
					}
				}
			}
		}
	}
	close FILE;
}

sub COMBINED_TWO {
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2,$c1c2_ref) = @_;
	print "TWO FIELD COMBINED\n";
	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );

	open(FILE, ">output/" . $gene . "_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
	if (( $gene eq "B" ) || ( $gene eq "C" )) {
		print FILE "Allele,Serotype,Qualifiers,COMMENT,ImputedAntigen,Broad,CIWD3.0,CWD2.0,EURCWD,Bw4/Bw6,C1/C2\n";
	}
	elsif ( $gene eq "A" ) {
		print FILE "Allele,Serotype,Qualifiers,COMMENT,ImputedAntigen,Broad,CIWD3.0,CWD2.0,EURCWD,Bw4/Bw6\n";
	}
	else {
		print FILE "Allele,Serotype,Qualifiers,COMMENT,ImputedAntigen,Broad,CIWD3.0,CWD2.0,EURCWD,DR5X\n";
	}
	my %twoField;

	foreach my $allele ( @$alleles_sorted_ref ) {		#go through all alleles
		my $twoField = "";
		my $qualifier = "I";
		unless (( exists $nullAllele_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
			if ( $allele =~ /(\w+\*\d+:\d+)/ ) {
				$twoField = $1;
			}
		}
		if ( exists $twoField{$twoField} ) {
			next;
		}
		if ( exists $nullAllele_ref->{ $allele } ) {		# Null
			next;
		}
		elsif ( exists $qallele_ref->{ $allele } ) {		# Questionable
			next;
		}
		elsif ( exists $assigned_ref->{ $allele } ) {	# assigned
			if ( $assigned_ref->{ $allele } =~ /(\S+_*\S*)_(LAX)/ ) {	# LAX
				my $group = $1;
				my $lax = $2;
				my $serotype = $group;
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				my $who = "";
				unless ( $base_ref->{ $group } =~ /Cw1[2-8]/ ) {
					$who = $base_ref->{ $group };
				}
				if ( $lax eq "LAX" ) {
					$lax = "SEROTYPE";
					$qualifier = "S";
				}
				if ( exists $cross_ref->{ $allele } ) {		# cross-reactivity
					print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . "_C," . $who . "," . $broad_ref->{ $group };
					my $cross = "";
					for ( my $index = 0; $index < scalar @{$cross_ref->{ $allele }}; $index++ ) {
						if ( $index == 0 ) {
							$cross = $cross_ref->{ $allele }->[$index];
						}
						else {
							$cross = $cross . "," . $cross_ref->{ $allele }->[$index];
						}
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE "," . $ciwd_ref->{ $twoField};
					}
					else {
						print FILE ",";
					}
					if ( exists $cwd_ref->{ $twoField } ) {
						print FILE "," . $cwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $ecwd_ref->{ $twoField } ) {
						print FILE "," . $ecwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if (( exists $bw_ref->{ $group } ) && ( !exists $bw_ref2->{ $allele } )) {
						print FILE "," . $bw_ref->{ $group };
					}
					elsif ( exists $bw_ref2->{ $allele }) {
						print FILE "," . $bw_ref2->{ $allele } . ",";
					}
					else {
						print FILE ",";
					}
					print FILE "," . $cross
				}
				else {	# no cross-reactive
					if (( $serotype eq "A-0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",A11,A11";
					}
					elsif (( $serotype eq "B-4406" ) && ( $allele =~ /B\*51/ )) {	# request from MFV on February 10 2021
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",B51,B5";
					}
					else {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . "," . $who . "," . $broad_ref->{ $group };
					}
					if ( exists $ciwd_ref->{ $twoField } ) {
						print FILE "," . $ciwd_ref->{ $twoField};
					}
					else {
						print FILE ",";
					}
					if ( exists $cwd_ref->{ $twoField } ) {
						print FILE "," . $cwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $ecwd_ref->{ $twoField } ) {
						print FILE "," . $ecwd_ref->{ $twoField };
					}
					else {
						print FILE ",";
					}
					if ( exists $bw_ref2->{ $allele }) {
						print FILE "," . $bw_ref2->{ $allele };
					}
					elsif ( exists $bw_ref->{ $group }) {
						print FILE "," . $bw_ref->{ $group };
					}
					else {
						print FILE ",";
					}
					if (( $gene eq "B" ) || ( $gene eq "C" )) {
						if ( $c1c2_ref->{ $allele } ) {
							print FILE "," . $c1c2_ref->{ $allele };
						}
						else {
							print FILE ",";
						}
					}
				}
				print FILE "\n";
			}
			else {		# stringent
				my $group = $assigned_ref->{ $allele };
				my $serotype = $group;
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				my $who = "";
				unless ( $base_ref->{ $group } =~ /Cw1[2-8]/ ) {
					$who = $base_ref->{ $group };
				}
				my $full = ",FULL,";
				$qualifier = "F";
				if (( $serotype eq "A-0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "A11,A11";
				}
				elsif (( $serotype eq "B-4406" ) && ( $allele =~ /B\*51/ )) {	# request from MFV on February 10 2021
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "B51,B5";
				}
				else {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . $who . "," . $broad_ref->{ $group };
				}
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField};
				}
				else {
					print FILE ",";
				}
				if ( exists $cwd_ref->{ $twoField } ) {
					print FILE "," . $cwd_ref->{ $twoField };
				}
				else {
					print FILE ",";
				}
				if ( exists $ecwd_ref->{ $twoField } ) {
					print FILE "," . $ecwd_ref->{ $twoField };
				}
				else {
					print FILE ",";
				}
				if ( exists $bw_ref->{ $group }) {
					print FILE "," . $bw_ref->{ $group };
				}
				else {
					print FILE ",";
				}
				if (( $gene eq "B" ) || ( $gene eq "C" )) {
					if ( $c1c2_ref->{ $allele } ) {
						print FILE "," . $c1c2_ref->{ $allele };
					}
					else {
						print FILE ",";
					}
				}
				print FILE "\n";
			}
		}
		elsif ( exists $short_ref->{ $allele } ) {	# short
			print FILE $twoField . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
				my $residue = 0;
				my $sero = $1 . "-" . $2;
				$sero =~ s/DRB1/DR/;
				my $serotype = "";
				my $group = "";
				my $test = 0;
				foreach my $short ( sort @{$short_ref->{ $allele }} ) {
					if ( $short =~ /(\S+)_(\d+)/ ) {
						$group = $1;
						$serotype = $group;
						$residue = $2;
						if ( $group =~ /$sero/ ) {		# allele name and sero type matches
							$test = 1;
							last;
						}
					}
				}
				
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				my $who = "";
				unless ( $base_ref->{ $group } =~ /Cw1[2-8]/ ) {
					$who = $base_ref->{ $group };
				}
				if ( $test == 1 ) {
					if ( $num == 1 ) {	# specific short
						print FILE $serotype . "," . $qualifier . ",SHORT," . $who . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE $serotype . "," . $qualifier . ",SHORT-CROSSREACTIVE," . $who . "," . $broad_ref->{ $group } . ",";
					}
				}
				else {		# no match => capture the first element
					my @tmp = sort @{$short_ref->{ $allele }};
					my $first_name = $tmp[0];
					if ( $first_name =~ /(\S+)_(\d+)/ ) {
						$group = $1;
						$residue = $2;
						$serotype = $group;
						if ( exists $antigen{ $group } ) {
							$serotype = $antigen{ $group };
						}
						unless ( $base_ref->{ $group } =~ /Cw1[2-8]/ ) {
							$who = $base_ref->{ $group };
						}
					}
					if ( $num == 1 ) {	# specific short
						print FILE $serotype . "," . $qualifier . ",SHORT," . $who . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE $serotype . "," . $qualifier . ",SHORT-CROSSREACTIVE," . $who . "," . $broad_ref->{ $group } . ",";
					}
				}

				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",";
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

				if (( $gene eq "B" ) || ( $gene eq "C" )) {	# HLA-B
					if ( exists $bw_ref->{ $group } ) {
						if ( !exists $bw_ref2->{ $allele } ) {
							if ( $bw_ref->{ $group } eq "Bw4" ) {
								if (( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Bw6" ) {
								if (( $residue == 76 ) ||( $residue == 82 ) || ( $residue == 83)) {
									print FILE "Negative,";
								}
								elsif (( $residue != 76) && ( $residue != 82 ) && ( $residue != 83)) {
									print FILE $bw_ref->{ $group } . ",";
								}
							}
							elsif ( $bw_ref->{ $group } eq "Negative" ) {
								print FILE $bw_ref->{ $group } . ",";
							}
							else {
								print FILE "Negative,";
							}
						}
						else {
							print FILE $bw_ref2->{ $allele } . ",";
						}
						if ( $c1c2_ref->{ $allele } ) {
							print FILE $c1c2_ref->{ $allele } . ",";
						}
						else {
							print FILE ",";
						}
					}
				}
				else {
					if (( exists $bw_ref->{ $group } ) && ( !exists $bw_ref2->{ $allele } )) {
						if ( $bw_ref->{ $group } eq "Bw4" ) {
							if (( $residue == 82 ) || ( $residue == 83)) {
								print FILE "Negative,";
							}
							elsif (( $residue != 82 ) && ( $residue != 83)) {
								print FILE $bw_ref->{ $group } . ",";
							}
						}
						elsif ( $bw_ref->{ $group } eq "Bw6" ) {
							if (( $residue == 76 ) ||( $residue == 82 ) || ( $residue == 83)) {
								print FILE "Negative,";
							}
							elsif (( $residue != 76) && ( $residue != 82 ) && ( $residue != 83)) {
								print FILE $bw_ref->{ $group } . ",";
							}
						}
						elsif ( $bw_ref->{ $group } eq "Negative" ) {
							print FILE $bw_ref->{ $group } . ",";
						}

						else {
							print FILE ",";
						}
					}
					elsif ( exists $bw_ref2->{ $allele }) {
						print FILE $bw_ref2->{ $allele } . ",";
					}
					else {
						print FILE ",";
					}
				}
			}
			
			my $index = 0;
			foreach my $short ( sort @{$short_ref->{ $allele }} ) {
				print FILE $short;
				if ( $index < $num - 1) {
					print FILE ",";
				}
				else {
					print FILE "\n";
				}
				$index++;
			}
		}
		else {		# no match
			print FILE $twoField;
			if (( $allele =~ /DRB1\*08:04:02/ ) || ( $allele =~ /DRB1\*04:20/ ) ||
			( $allele =~ /DQB1\*05:03:02/ ) || ( $allele =~ /DQB1\*06:01:02/ ) || ( $allele =~ /DQB1\*06:05:02/ ) ||
			( $allele =~ /DQB1\*06:06/ )) {	# missing key residues 9 - 14
				if ( $allele =~ /DRB1\*08:04:02/ ) {
					print FILE ",DR-0801,F,InSilico,DR8,DR8,";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",DR-0403,S,InSilico,DR4,DR4";
				}
				elsif ( $allele =~ /DQB1\*05:03:02/ ) {
					print FILE ",DQ5,S,InSilico,DQ5,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:01:02/ ) {
					print FILE ",DQ6,F,InSilico,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:05:02/ ) {
					print FILE ",DQ6,F,InSilico,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:06/ ) {
					print FILE ",DQ6,S,InSilico,DQ6,DQ1";
				}
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",,";
				}
				if ( exists $cwd_ref->{ $twoField } ) {
					print FILE $cwd_ref->{ $twoField } . ",";
				}
				else {
					print FILE ",";
				}
				if ( exists $ecwd_ref->{ $twoField } ) {
					print FILE $ecwd_ref->{ $twoField } . "\n";
				}
				else {
					print FILE "\n";
				}

			}
			else {
				print FILE ",UNA,U,None,None,None";
				if ( exists $ciwd_ref->{ $twoField } ) {
					print FILE "," . $ciwd_ref->{ $twoField} . ",";
				}
				else {
					print FILE ",,";
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
				if ( exists $bw_ref2->{ $allele }) {
					print FILE $bw_ref2->{ $allele } . ",";
				}
				else {
					if (( $gene eq "B" ) || ( $gene eq "C" )) {
						print FILE "Negative,";
					}
					else {
						print FILE ",";
					}
				}
				if (( $gene eq "B" ) && ( $c1c2_ref->{ $allele } )) {
					print FILE $c1c2_ref->{ $allele } . "\n";
				}
				else {
					print FILE "\n";
				}
			}
		}
		$twoField{$twoField} = 0;
	}
	close FILE;
}

sub ANTIGEN {
	return $antigen_ref;
}

1;
