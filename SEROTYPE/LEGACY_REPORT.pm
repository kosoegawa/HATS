#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: LEGACY_REPORT.pm 
# This module was developed to print lagacy table
# Imputed column:
# Associated Antigen A0203 & A210 => Legacy A203 and A210
# serotype A0305, but allelel name A*11, e.g., A*11:04 => A11
# A2410 & A2423 => A2403 (closest associated antigen)
# serotype B4406, but allele name B*51, e.g., B*51:23 => B51
# last reviewed, modified and documented on February 19 2026

package LEGACY_REPORT;
use strict;
use GROUP_SORT;
use COMBINE;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

my @broad = ("A9","A10","A19","A28","B5","B12","B14","B15","B16","B17","B21","B22","B40","B70","Cw3","DR2","DR3","DR5","DR6","DQ1","DQ3");

# print combined table
sub COMBINED {
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2,$dr_ref) = @_;
	print "COMBINED\n";

	my $antigen_ref = COMBINE::ANTIGEN();

	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically

	open(FILE, ">output/" . $gene . "_Legacy_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
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
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
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
					if (( $serotype eq "A0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
						print FILE $allele . "," . $lax . "," . $serotype . ",A11,A11";
					}
					elsif ( $serotype eq "A0203" ) {
						print FILE $allele . "," . $lax . "," . $serotype . ",A203,A2";
					}
					elsif ( $serotype eq "A0210" ) {
						print FILE $allele . "," . $lax . "," . $serotype . ",A210,A2";
					}
					elsif ( $serotype eq "B0703" ) {
						print FILE $allele . "," . $lax . "," . $serotype . ",B703,B7";
					}
					elsif ( $serotype eq "DR0103" ) {
						print FILE $allele . "," . $lax . "," . $serotype . ",DR103,DR1";
					}
					elsif (( $serotype eq "A2403" ) || ( $serotype eq "A2410" ) || ( $serotype eq "A2423" )) {
						print FILE $allele . "," . $lax . "," . $serotype . ",A2403,A9";
					}
					elsif (( $serotype eq "B4406" ) && ( $allele =~ /B\*51/ )) {	
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
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
				}
				my $full = ",FULL,";
				if (( $serotype eq "A0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
					print FILE $allele . $full . $serotype . ",A11,A11";
				}
				elsif ( $serotype eq "A0203" ) {
					print FILE $allele . $full . $serotype . ",A203,A2";
				}
				elsif ( $serotype eq "A0210" ) {
					print FILE $allele . $full . $serotype . ",A210,A2";
				}
				elsif ( $serotype eq "B0703" ) {
					print FILE $allele . $full . $serotype . ",B703,B7";
				}
				elsif ( $serotype eq "DR0103" ) {
					print FILE $allele . $full . $serotype . ",DR103,DR1";
				}
				elsif (( $serotype eq "A2403" ) || ( $serotype eq "A2410" ) || ( $serotype eq "A2423" )) {
					print FILE $allele . $full . $serotype . ",A2403,A9";
				}
				elsif (( $serotype eq "B4406" ) && ( $allele =~ /B\*51/ )) {	
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
				my $sero = $1 . $2;
				$sero =~ s/DRB1/DR/;
				$sero =~ s/C/Cw/;	#added here 2/20/26 to fix bug
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
				
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
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
						if ( exists $antigen_ref->{ $group } ) {
							$serotype = $antigen_ref->{ $group };
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
					print FILE ",InSilico,DR0801,DR8,DR8";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",InSilico,DR0403,DR4,DR4";
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
				print FILE ",UNA,U,None,None";
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
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2,$c1c2_ref) = @_;
	print "TWO FIELD COMBINED\n";

	my $antigen_ref = COMBINE::ANTIGEN();

	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );

	open(FILE, ">output/" . $gene . "_Legacy_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
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
				if ( $twoField eq "DPB1*26:01" ) {	# DPB1*26:01:01 is partial sequence, but DPB1*26:01:02:01 is full length
					next;
				}
				my $group = $1;
				my $lax = $2;
				my $serotype = $group;
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
				}
				if ( $lax eq "LAX" ) {
					$lax = "SEROTYPE";
					$qualifier = "S";
				}
				if ( exists $cross_ref->{ $allele } ) {		# cross-reactivity
					print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . "_C," . $base_ref->{ $group } . "," . $broad_ref->{ $group };
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
					if (( $serotype eq "A0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",A11,A11";
					}
					elsif ( $serotype eq "A0203" ) {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",A203,A2";
					}
					elsif ( $serotype eq "A0210" ) {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",A210,A2";
					}
					elsif ( $serotype eq "B0703" ) {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",B703,B7";
					}
					elsif ( $serotype eq "DR0103" ) {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",DR103,DR1";
					}
					elsif (( $serotype eq "A2403" ) || ( $serotype eq "A2410" ) || ( $serotype eq "A2423" )) {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",A2403,A9";
					}
					elsif (( $serotype eq "B4406" ) && ( $allele =~ /B\*51/ )) {	# request from MFV on February 10 2021
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . ",B51,B5";
					}
					else {
						print FILE $twoField . "," . $serotype . "," . $qualifier . "," . $lax . "," . $base_ref->{ $group } . "," . $broad_ref->{ $group };
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
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
				}
				my $full = ",FULL,";
				$qualifier = "F";
				if (( $serotype eq "A0305" ) && ( $allele =~ /A\*11/ )) {	# request from MFV on February 10 2021
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "A11,A11";
				}
				elsif ( $serotype eq "A0203" ) {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "A203,A2";
				}
				elsif ( $serotype eq "A0210" ) {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "A210,A2";
				}
				elsif ( $serotype eq "B0703" ) {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "B703,B7";
				}
				elsif ( $serotype eq "DR0103" ) {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "DR103,DR1";
				}
				elsif (( $serotype eq "A2403" ) || ( $serotype eq "A2410" ) || ( $serotype eq "A2423" )) {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "A2403,A9";
				}
				elsif (( $serotype eq "B4406" ) && ( $allele =~ /B\*51/ )) {	# request from MFV on February 10 2021
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . "B51,B5";
				}
				else {
					print FILE $twoField . "," . $serotype . "," . $qualifier . $full . $base_ref->{ $group } . "," . $broad_ref->{ $group };
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
			if ( $twoField eq "DPB1*26:01" ) {	# DPB1*26:01:01 is partial sequence, but DPB1*26:01:02:01 is full length
				next;
			}
			print FILE $twoField . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
				my $residue = 0;
#				my $sero = $1 . "-" . $2;
				my $sero = $1 . $2;
				$sero =~ s/DRB1/DR/;
				$sero =~ s/C/Cw/;	#added here 2/20/26 to fix bug
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
				
				if ( exists $antigen_ref->{ $group } ) {
					$serotype = $antigen_ref->{ $group };
				}
				if ( $test == 1 ) {
					if ( $num == 1 ) {	# specific short
						print FILE $serotype . "," . $qualifier . ",SHORT," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE $serotype . "," . $qualifier . ",SHORT-CROSSREACTIVE," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
				}
				else {		# no match => capture the first element
					my @tmp = sort @{$short_ref->{ $allele }};
					my $first_name = $tmp[0];
					if ( $first_name =~ /(\S+)_(\d+)/ ) {
						$group = $1;
						$residue = $2;
						$serotype = $group;
						if ( exists $antigen_ref->{ $group } ) {
							$serotype = $antigen_ref->{ $group };
						}
					}
					if ( $num == 1 ) {	# specific short
						print FILE $serotype . "," . $qualifier . ",SHORT," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
					}
					else {		# cross-reactive short
						print FILE $serotype . "," . $qualifier . ",SHORT-CROSSREACTIVE," . $base_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
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
				else {	# not B or C
					if (( exists $bw_ref->{ $group } ) && ( !exists $bw_ref2->{ $allele } )) {
						if ( $bw_ref->{ $group } eq "Bw4" ) {
							if (( $residue == 82 ) || ( $residue == 83)) {
								print FILE ",";		# deleted Negative for HLA-A
							}
							elsif (( $residue != 82 ) && ( $residue != 83)) {
								print FILE $bw_ref->{ $group } . ",";
							}
						}
						elsif ( $bw_ref->{ $group } eq "Bw6" ) {
							if (( $residue == 76 ) ||( $residue == 82 ) || ( $residue == 83)) {
								print FILE ",";		# deleted Negative for HLA-A
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
					print FILE ",DR0801,F,InSilico,DR8,DR8,";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",DR0403,S,InSilico,DR4,DR4";
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
				print FILE ",UNA,U,,None,None";
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


1;
