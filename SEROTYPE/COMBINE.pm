#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: COMBINE.pm 
# This module was developed to print table
# last reviewed, modified and documented on February 19 2026

package COMBINE;
use strict;
use GROUP_SORT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

my %antigen;
my $antigen_ref = \%antigen;
$antigen{ "A0101" } = "A1";
$antigen{ "A2501" } = "A25";
$antigen{ "A1101" } = "A11";
$antigen{ "A7401" } = "A74";
$antigen{ "A6901" } = "A69";
$antigen{ "A3601" } = "A36";
$antigen{ "A4301" } = "A43";
$antigen{ "A8001" } = "A80";

$antigen{ "B5201" } = "B52";
$antigen{ "B4501" } = "B45";
$antigen{ "B1302" } = "B13";
$antigen{ "B1401" } = "B64";
$antigen{ "B1402" } = "B65";
$antigen{ "B1512" } = "B76";
$antigen{ "B1513" } = "B77";
$antigen{ "B5701" } = "B57";
$antigen{ "B5801" } = "B58";
$antigen{ "B4901" } = "B49";
$antigen{ "B5001" } = "B50";
$antigen{ "B5401" } = "B54";
$antigen{ "B4101" } = "B41";
$antigen{ "B4201" } = "B42";
$antigen{ "B4601" } = "B46";
$antigen{ "B5301" } = "B53";
$antigen{ "B5901" } = "B59";
$antigen{ "B7301" } = "B73";
$antigen{ "B8101" } = "B81";
$antigen{ "B8201" } = "B82";

$antigen{ "Cw0102" } = "Cw1";
$antigen{ "Cw0202" } = "Cw2";
$antigen{ "Cw0303" } = "Cw9";
$antigen{ "Cw1402" } = "Cw14";
$antigen{ "Cw1701" } = "Cw17";
$antigen{ "Cw1801" } = "Cw18";

$antigen{ "DR0301" } = "DR17";
$antigen{ "DR0302" } = "DR18";
$antigen{ "DR0701" } = "DR7";
$antigen{ "DR0901" } = "DR9";
$antigen{ "DR1001" } = "DR10";

$antigen{ "DQ0201" } = "DQ2";
$antigen{ "DQ0301" } = "DQ7";
$antigen{ "DQ0302" } = "DQ8";
$antigen{ "DQ0303" } = "DQ9";
$antigen{ "DQ0401" } = "DQ4";
$antigen{ "DQ0501" } = "DQ5";
$antigen{ "DQ0602" } = "DQ6";	#comment out for SAB

my @broad = ("A9","A10","A19","A28","B5","B12","B14","B15","B16","B17","B21","B22","B40","B70","Cw3","DR2","DR3","DR5","DR6","DQ1","DQ3");

# print combined table
sub COMBINED {
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref) = @_;
	print "COMBINED\n";
	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );	# sort allele numerically

	open(FILE, ">output/" . $gene . "_Allele_Antigen_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
	print FILE "Allele,Qualifier,Associated,Split,Broad,CIWD3.0,CWD2.0,EURCWD\n";

	foreach my $allele ( @$alleles_sorted_ref ) {		#go through all alleles
		my $twoField = "";
		unless (( exists $nullAllele_ref->{ $allele } ) || ( exists $qallele_ref->{ $allele })) {		# Null
			if ( $allele =~ /($gene\*\d+:\d+)/ ) {
				$twoField = $1;
			}
		}
		if ( exists $nullAllele_ref->{ $allele } ) {		# Null
			print FILE $allele . "," . $nullAllele_ref->{ $allele } . "\n";
		}
		elsif ( exists $qallele_ref->{ $allele } ) {		# Questionable
			print FILE $allele . "," . $qallele_ref->{ $allele } . "\n";
		}
		elsif ( exists $assigned_ref->{ $allele } ) {	# assigned
			if ( $assigned_ref->{ $allele } =~ /(\S+_*\S*)_(LAX)/ ) {	# LAX
				my $group = $1;
#				my $serotype = $1;
				my $lax = $2;
				my $serotype = $group;
				if ( exists $antigen{ $serotype } ) {
					$serotype = $antigen{ $serotype };
				}
				if ( $lax eq "LAX" ) {
					$lax = "S";		############
				}
				if ( exists $cross_ref->{ $allele } ) {		# cross-reactivity
					print FILE $allele . "," . $serotype . "," . $lax . "_C," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
					print FILE $cross
				}
				else {	# no cross-reactive
					if ( $serotype eq $broad_ref->{ $group } ) {	# No ASSOCIATED & SPLIT, e.g., A36 
						print FILE $allele . ",". $lax  . "," . "," . "," . $broad_ref->{ $group };
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {	# No associated antigen, e.g., A25
						if ( $serotype eq "B4005" ) {
							print FILE $allele . "," . $lax . "," . $serotype . "," . "," . $broad_ref->{ $group };
						}
						else {
							print FILE $allele . "," . $lax . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE $allele . "," . $lax . "," . $serotype . "," . "," . $broad_ref->{ $group };
					}
					else {
						print FILE $allele . "," . $lax . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
				}
				print FILE "\n";
			}
			else {		# stringent
				my $group = $assigned_ref->{ $allele };
				my $serotype = $group;
				if ( exists $antigen{ $serotype } ) {
					$serotype = $antigen{ $group };
				}
				my $full = "F";
				if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
					print FILE $allele . "," . $full . "," . "," . "," . $broad_ref->{ $group };
				}
				elsif ( $serotype eq $parent_ref->{ $group } ) {
					if ( $serotype eq "B4005" ) {
						print FILE $allele . "," . $full . "," . $serotype . "," . "," . $broad_ref->{ $group };
					}
					else {
						print FILE $allele . "," . $full . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
					}
				}
				elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
					print FILE $allele . "," . $full . "," . $serotype . "," . "," . $broad_ref->{ $group };
				}
				else {
					print FILE $allele . "," . $full . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
				print FILE "\n";
			}
		}
		elsif ( exists $short_ref->{ $allele } ) {	# short
			print FILE $allele . ",";
			my $num = scalar @{$short_ref->{ $allele }};
			if ( $allele =~ /(\S+)\*(\d+):\d+:*\d*:*\d*/ ) {		#[1-9]+0* was important, B40
				my $residue = 0;
				my $sero = $1 . $2;	# need to modify here
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
				
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				# go through three different conditions: $base_type and $type is identical
				if ( $test == 1 ) {
					if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
						print FILE "I,," . "," . $broad_ref->{ $group } . ",";
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {	#ASSOCIATED BLANK
						if ( $serotype eq "B4005" ) {
							print FILE "I," . $serotype . ",," . $broad_ref->{ $group } . ",";
						}
						else {
							print FILE "I,," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE  "I," . $serotype .",," . $broad_ref->{ $group } . ",";
					}
					else {
						print FILE "I," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
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
					if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
						print FILE "I,," . "," . $broad_ref->{ $group } . ",";
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {	#ASSOCIATED BLANK
						if ( $serotype eq "B4005" ) {
							print FILE "I," . $serotype . ",," . $broad_ref->{ $group } . ",";
						}
						else {
							print FILE "I,," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE  "I," . $serotype .",," . $broad_ref->{ $group } . ",";
					}
					else {
						print FILE "I," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
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
					print FILE ",InSilico,DR0801,,DR8";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",InSilico,DR0403,,DR4";
				}
				elsif ( $allele =~ /DQB1\*05:03:02/ ) {	# sequence extended
					print FILE ",InSilico,,DQ5,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:01:02/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:05:02/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:06/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
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
				print FILE ",U,UNA,None,None";
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
				print FILE "\n";

			}
		}
	}
	close FILE;
}

sub COMBINED_TWO {
	my ($database,$nullAllele_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
	$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref) = @_;
	print "TWO FIELD COMBINED\n";
	my @combined;
	my $combined_ref = \@combined;
	push @combined, keys %$nullAllele_ref;
	push @combined, keys %$qallele_ref;
	push @combined, keys %$assigned_ref;
	push @combined, @$unassigned_ref;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $combined_ref );

	open(FILE, ">output/" . $gene . "_Protein_Antigen_Table_IMGT_HLA_" . $database . "_" . $date . ".csv");
	print FILE "Protein,Qualifier,Associated,Split,Broad,CIWD3.0,CWD2.0,EURCWD\n";
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
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				if ( $lax eq "LAX" ) {
					$qualifier = "S";
				}
				if ( exists $cross_ref->{ $allele } ) {		# cross-reactivity
					print FILE $twoField . "," . $serotype . "," . $qualifier . "_C," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
					print FILE "," . $cross
				}
				else {	# no cross-reactive
					if ( $serotype eq $broad_ref->{ $group } ) {	# No ASSOCIATED & SPLIT, e.g., A36 
						print FILE $twoField . "," . $qualifier . "," . "," . "," . $broad_ref->{ $group };
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {	# No associated antigen, e.g., A25
						if ( $serotype eq "B4005" ) {
							print FILE $twoField . "," . $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group };
						}
						else {
							print FILE $twoField . "," . $qualifier . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE $twoField . "," . $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group };
					}
					else {
						print FILE $twoField . "," . $qualifier . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
				}
				print FILE "\n";
			}
			else {		# stringent
				my $group = $assigned_ref->{ $allele };
				my $serotype = $group;
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				$qualifier = "F";
				if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
					print FILE $twoField . "," . $qualifier . "," . "," . "," . $broad_ref->{ $group };
				}
				elsif ( $serotype eq $parent_ref->{ $group } ) {
					if ( $serotype eq "B4005" ) {
						print FILE $twoField . "," . $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group };
					}
					else {
						print FILE $twoField . "," . $qualifier . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
					}
				}
				elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
					print FILE $twoField . "," . $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group };
				}
				else {
					print FILE $twoField . "," . $qualifier . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group };
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
				
				if ( exists $antigen{ $group } ) {
					$serotype = $antigen{ $group };
				}
				if ( $test == 1 ) {
					if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
						print FILE $qualifier . "," . "," . "," . $broad_ref->{ $group } . ",";
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {
						if ( $serotype eq "B4005" ) {
							print FILE $qualifier . "," .$serotype . "," .  "," . $broad_ref->{ $group } . ",";
						}
						else {
							print FILE $qualifier . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group } . ",";
					}
					else {
						print FILE $qualifier . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
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
					}

					if ( $serotype eq $broad_ref->{ $group } ) {	# ASSOCIATED & SPLIT: BLANK
						print FILE $qualifier . "," . "," . "," . $broad_ref->{ $group } . ",";
					}
					elsif ( $serotype eq $parent_ref->{ $group } ) {
						if ( $serotype eq "B4005" ) {
							print FILE $qualifier . "," .$serotype . "," .  "," . $broad_ref->{ $group } . ",";
						}
						else {
							print FILE $qualifier . "," . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
						}
					}
					elsif ( $parent_ref->{ $group } eq $broad_ref->{ $group } ) {	# No split, e.g., A2
						print FILE $qualifier . "," . $serotype . "," . "," . $broad_ref->{ $group } . ",";
					}
					else {
						print FILE $qualifier . "," . $serotype . "," . $parent_ref->{ $group } . "," . $broad_ref->{ $group } . ",";
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
					print FILE ",InSilico,DR0801,,DR8";
				}
				elsif ( $allele =~ /DRB1\*04:20/ ) {	# missing key residues 9 - 14
					print FILE ",InSilico,DR0403,,DR4";
				}
				elsif ( $allele =~ /DQB1\*05:03:02/ ) {
					print FILE ",InSilico,,DQ5,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:01:02/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:05:02/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
				}
				elsif ( $allele =~ /DQB1\*06:06/ ) {
					print FILE ",InSilico,,DQ6,DQ1";
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
				print FILE ",U,UNA,None,None";
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
				print FILE "\n";
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
