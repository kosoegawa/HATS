#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: COUNT.pm 
# This module was developed to generate Summary table
# last modified and documented on December 7 2024

#no warnings 'experimental::smartmatch';
#eliminated smartmatch
package COUNT;
use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
use Openfile;
use GROUP_SORT;
use ASSIGNED_SHORT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
#my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

sub COUNT {		# deal with remaining serotypes with strict mode
	my ($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref) = @_;

	my %unique;
	open (FILE, ">output/Summary_" . $gene . "_" . $date . ".csv");
	print FILE "WHOAssigned,Serotype,STRINGENT,STR_OUT,LAX,LAX_OUT,UNASSIGNED\n";
	my $str_total = 0;
	my $str_out_total = 0;
	my $lax_total = 0;
	my $lax_out_total = 0;
	my $unassign_total = 0;
	my %outliers;

	foreach my $basetype ( @$basetype_ref ) {	# e.g., @basetype = ("DQ5","DQ6","DQ2","DQ7","DQ8","DQ9","DQ4");
		foreach my $type ( @{$sero_ref->[0]} ) {	# go through serotype
			unless ( exists $base_ref->{ $type } ) {
				print $type . "PROBLEM: PAY ATTENTION!!!\n";
			}
			if ( $basetype eq $base_ref->{ $type } ) {
				print FILE $basetype . ",";
#				print $type . "\n";
				my $str = 0;
				my $str_out = 0;
				my $lax = 0;
				my $lax_out = 0;
				my $unassign = 0;
				foreach my $csv (@$csv_ref) {	# go through CSV files
					my @list = Openfile::open_file_from_list($csv);
					my $name = $type . "_";
					if ( $csv =~ /$name/ ) {	# target file, DR4.2 takes DR402 and DR412: this is causing problem
						if ( $csv =~ /LAX/ ) {	# LAX file
							$lax = scalar @list;	# number of lines
#							print $lax . "\n";
							my $out = 0;
							my @matched;
							foreach my $list ( @list ) {	# go through allele list
								if ( exists $unique{ $list } ) {
									print $type ."\tDuplicated: " . $csv . "\t" . $list . "\n";
								}
								$unique{ $list } = 0;
								if ( $list =~ /$key_ref->{$type}/ ) {
									$out++;		# allele group matches
								#	print $key_ref->{$type} . "\n";
									push @matched, $list;
								}
							}
							$lax_out = $lax - $out;	# outlier
							if ( $lax_out != 0 ) {	# deal with outlier
								foreach my $list ( @list ) {
									my $matched_count = 0;
									foreach my $matched ( @matched ) {
										if ( $list eq $matched ) {
											$matched_count++;
										}
									}
									if ( $matched_count == 0 ) {
										$outliers{ $list } = $type;
									}
								}
							}
						}
						else {	# strict
							$str = scalar @list;
							my $out = 0;
							my @matched;
							foreach my $list ( @list ) {
								if ( exists $unique{ $list } ) {
									print $type . "\tDuplicated: " . $csv . "\t" . $list . "\n";
								}
								$unique{ $list } = 0;
#								print $type . "\n";
								if ( $list =~ /$key_ref->{$type}/ ) {
									$out++;
								#	print $key_ref->{$type} . "\n";
									push @matched, $list;
								}
							}
							$str_out = $str - $out;
							if ( $str_out != 0 ) {	# deal with outlier
								foreach my $list ( @list ) {
									my $matched_count = 0;
									foreach my $matched ( @matched ) {
										if ( $list eq $matched ) {
											$matched_count++;
										}
									}
									if ( $matched_count == 0 ) {
										$outliers{ $list } = $type;
									}
								}
							}
						}
					}
					elsif ( $csv =~ /unassigned_$gene/ ) {
						my $count = 0;
						foreach my $subtype ( @{$sero_ref->[1]} ) {
							if ( $subtype eq $type) {
								$count++;
							}
						}
						if ( $count == 0 ) {
							foreach my $list ( @list ) {
								if ( $list =~ /$key_ref->{$type}/ ) {
									$unassign++;
								}
							}
						}
					}
				}
				print FILE  $type . "," . $str . "," . $str_out . "," . $lax . "," . $lax_out . "," . $unassign . "\n";
				$str_total += $str;
				$str_out_total += $str_out;
				$lax_total += $lax;
				$lax_out_total += $lax_out;
				$unassign_total += $unassign;
			}
		}
	}

	print FILE  "Total,," . $str_total . "," . $str_out_total . "," . $lax_total . "," . $lax_out_total . "," . $unassign_total . "\n";
	my $null_count = scalar keys %$null_ref;
	my $qallele_count = scalar keys %$qallele_ref;
	my $expressed = $str_total + $lax_total + $unassign_total;
	my $sum = $null_count + $qallele_count + $expressed;
	print FILE "Total alleles," . $sum . "\n";
	print FILE "Null alleles," . $null_count . "\n";
	print FILE "Q alleles," . $qallele_count . "\n";
	print FILE "Expressed alleles," . $expressed . "\n";
	my $unique_num = scalar (keys %unique);
	if ($unique_num != ($str_total + $lax_total)) {
		my $diff = $str_total + $lax_total - $unique_num;
		print $diff . " Duplicated allele exists!\n";
	}
	close FILE;

	open ( FILE, ">output/Outliers_" . $gene . "_" . $date . ".csv" );
	print FILE "Allele,Type\n";
	my @alleles = keys %outliers;
	my $alleles_ref = \@alleles;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $alleles_ref );
	foreach my $outliers ( @$alleles_sorted_ref ) {
		print FILE $outliers . "," . $outliers{ $outliers } . "\n";
	}
	close FILE;

}

sub TWOFIELD {		# deal with remaining serotypes with strict mode
	my ($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref) = @_;

	my %unique;
	open (FILE, ">output/TwoField_Summary_" . $gene . "_" . $date . ".csv");
	print FILE "BASE,TYPE,STRINGENT,STR_OUT,LAX,LAX_OUT,UNASSIGNED\n";
	my $str_total = 0;
	my $str_out_total = 0;
	my $lax_total = 0;
	my $lax_out_total = 0;
	my $unassign_total = 0;
	my %outliers;

	foreach my $basetype ( @$basetype_ref ) {
		foreach my $type ( @{$sero_ref->[0]} ) {	# go through serotype
			if ( $basetype eq $base_ref->{ $type } ) {
				print FILE $basetype . ",";
#				print $type . "\n";
				my $str = 0;
				my $str_out = 0;
				my $lax = 0;
				my $lax_out = 0;
				my $unassign = 0;
				foreach my $csv (@$csv_ref) {	# go through CSV files
					my @list = Openfile::open_file_from_list($csv);

					my @twofield;
					my %twofield;
					foreach my $list ( @list ) {
						my $twofield = "";
						if ( $list =~ /(\w+\*\d+:\d+)/ ) {
							$twofield = $1;
						}
						unless ( exists $twofield{ $twofield } ) {
							unless ($twofield eq "") {
								$twofield{ $twofield } = 0;
								push @twofield, $twofield;
							}
						}
					}

					my $name = $type . "_";
					if ( $csv =~ /$name/ ) {	# target file
						if ( $csv =~ /LAX/ ) {	# LAX file
							$lax = scalar @twofield;	# number of lines
							my $out = 0;
							my @matched;
							foreach my $twofield ( @twofield ) {
								if ( exists $unique{ $twofield } ) {
									print "TwoField Duplicated: " . $csv . "\t" . $twofield . "\n";
								}
								$unique{ $twofield } = 0;
								if ( $twofield =~ /$key_ref->{$type}/ ) {
									$out++;		# allele group matches
								#	print $key_ref->{$type} . "\n";
									push @matched, $twofield;
								}
							}
							$lax_out = $lax - $out;	# outlier
							if ( $lax_out != 0 ) {	# deal with outlier
								foreach my $twofield ( @twofield ) {
									my $matched_count = 0;
									foreach my $matched ( @matched ) {
										if ( $twofield eq $matched ) {
											$matched_count++;
										}
									}
									if ( $matched_count == 0 ) {
										$outliers{ $twofield } = $type;
									}
								}
							}
						}
						else {	# strict
							$str = scalar @twofield;
							my $out = 0;
							my @matched;
							foreach my $twofield ( @twofield ) {
								if ( exists $unique{ $twofield } ) {
									print "TwoField Duplicated: " . $csv . "\t" . $twofield . "\n";
								}
								$unique{ $twofield } = 0;
								if ( $twofield =~ /$key_ref->{$type}/ ) {
									$out++;
								#	print $key_ref->{$type} . "\n";
									push @matched, $twofield;
								}
							}
							$str_out = $str - $out;
							if ( $str_out != 0 ) {	# deal with outlier
								foreach my $twofield ( @twofield ) {
									my $matched_count = 0;
									foreach my $matched ( @matched ) {
										if ( $twofield eq $matched ) {
											$matched_count++;
										}
									}
									if ( $matched_count == 0 ) {
										$outliers{ $twofield } = $type;
									}
								}
							}
						}
					}
					elsif ( $csv =~ /unassigned_$gene/ ) {
						my $count = 0;
						foreach my $subtype ( @{$sero_ref->[1]} ) {
							if ( $subtype eq $type) {
								$count++;
							}
						}
						if ( $count == 0 ) {
							foreach my $twofield ( @twofield ) {
								if ( $twofield =~ /$key_ref->{$type}/ ) {
									$unassign++;
								}
							}
						}
					}
				}
				print FILE  $type . "," . $str . "," . $str_out . "," . $lax . "," . $lax_out . "," . $unassign . "\n";
				$str_total += $str;
				$str_out_total += $str_out;
				$lax_total += $lax;
				$lax_out_total += $lax_out;
				$unassign_total += $unassign;
			}
		}
	}

	print FILE  "Total,," . $str_total . "," . $str_out_total . "," . $lax_total . "," . $lax_out_total . "," . $unassign_total . "\n";
	my $null_count = scalar keys %$null_ref;
	my $qallele_count = scalar keys %$qallele_ref;
	my $expressed = $str_total + $lax_total + $unassign_total;
	my $sum = $null_count + $qallele_count + $expressed;
	print FILE "Total two-field alleles," . $sum . "\n";
	print FILE "Null alleles," . $null_count . "\n";
	print FILE "Q alleles," . $qallele_count . "\n";
	print FILE "Expressed alleles," . $expressed . "\n";
	my $unique_num = scalar (keys %unique);
	if ($unique_num != ($str_total + $lax_total)) {
		my $diff = $str_total + $lax_total - $unique_num;
		print $diff . " Duplicated allele exists!\n";
	}
	close FILE;

	open ( FILE, ">output/Two_Field_Outliers_" . $gene . "_" . $date . ".csv" );
	print FILE "Allele,Type\n";
	my @alleles = keys %outliers;
	my $alleles_ref = \@alleles;
	my $alleles_sorted_ref = GROUP_SORT::SORT( $alleles_ref );
	foreach my $outliers ( @$alleles_sorted_ref ) {
		print FILE $outliers . "," . $outliers{ $outliers } . "\n";
	}
	close FILE;
}

sub SUMMARY {
	my ( $csv_ref, $gene, $null_ref, $qallele_ref, $whotype_ref ) = @_;
	my @list = Openfile::open_file_from_list($csv_ref);
	my @who_values = values %$whotype_ref;

	my $header = shift @list;
	my %full;
	my $full_ref = \%full;

	my %who;
	my $who_ref = \%who;

	my $common_full = 0;
	my $inter_full = 0;
	my $wd_full = 0;
		
	my $common_sero = 0;
	my $inter_sero = 0;
	my $wd_sero = 0;
		
	my $common_short = 0;
	my $inter_short = 0;
	my $wd_short = 0;
		
	my $common_sc = 0;
	my $inter_sc = 0;
	my $wd_sc = 0;
		
	my $common_ins = 0;	# in silico
	my $inter_ins = 0;
	my $wd_ins = 0;

	foreach my $line ( @list ) {
		my @elements = split( ",", $line );
		my $type = $elements[2];
		if (( $type eq "NULL" ) || ( $type eq "Questionable" )) {
			next;
		}
		unless ( exists $full_ref->{ $type } ) {
			$full_ref->{ $type }->[0] = 0;		# FULL
			$full_ref->{ $type }->[1] = 0;		# Serotype
			$full_ref->{ $type }->[2] = 0;		# S
			$full_ref->{ $type }->[3] = 0;		# SC
			$full_ref->{ $type }->[4] = 0;		# InSilico

			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[0] = 0;		# FULL
				$who_ref->{ $type }->[1] = 0;		# Serotype
				$who_ref->{ $type }->[2] = 0;		# S
				$who_ref->{ $type }->[3] = 0;		# SC
				$who_ref->{ $type }->[4] = 0;		# InSilico
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[0] = 0;		# FULL
				$who_ref->{ $type }->[1] = 0;		# Serotype
				$who_ref->{ $type }->[2] = 0;		# S
				$who_ref->{ $type }->[3] = 0;		# SC
				$who_ref->{ $type }->[4] = 0;		# InSilico
			}
		}

		if ( $elements[1] eq "UNA" ) {
			$full_ref->{ $type }->[0] = $full_ref->{ $type }->[0] + 1;
			next;
		}
		elsif ( $elements[1] eq "FULL" ) {
			$full_ref->{ $type }->[0] = $full_ref->{ $type }->[0] + 1;
			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[0] = $who_ref->{ $type }->[0] + 1;
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[0] = $who_ref->{ $type }->[0] + 1;
			}
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_full++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_full++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_full++;
			}
			next;
		}
		elsif ( $elements[1] eq "SEROTYPE" ) {
			$full_ref->{ $type }->[1] = $full_ref->{ $type }->[1] + 1;
			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[1] = $who_ref->{ $type }->[1] + 1;
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[1] = $who_ref->{ $type }->[1] + 1;
			}
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_sero++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_sero++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_sero++;
			}
			next;
		}
		elsif ( $elements[1] eq "S" ) {
			$full_ref->{ $type }->[2] = $full_ref->{ $type }->[2] + 1;
			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[2] = $who_ref->{ $type }->[2] + 1;
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[2] = $who_ref->{ $type }->[2] + 1;
			}
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_short++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_short++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_short++;
			}
			next;
		}
		elsif ( $elements[1] eq "SC" ) {
			$full_ref->{ $type }->[3] = $full_ref->{ $type }->[3] + 1;
			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[3] = $who_ref->{ $type }->[3] + 1;
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[3] = $who_ref->{ $type }->[3] + 1;
			}
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_sc++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_sc++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_sc++;
			}
			next;
		}
		elsif ( $elements[1] eq "InSilico" ) {
			$full_ref->{ $type }->[4] = $full_ref->{ $type }->[4] + 1;
			if ( grep ( /^$type$/, @who_values ) ) {
				$who_ref->{ $type }->[4] = $who_ref->{ $type }->[4] + 1;
			}
			elsif ( exists $whotype_ref->{ $type } ) {
				$who_ref->{ $type }->[4] = $who_ref->{ $type }->[4] + 1;
			}
			if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
				$common_ins++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
				$inter_ins++;
			}
			elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
				$wd_ins++;
			}
			next;
		}
	}

	open(FILE, ">output/" . $gene . "_AlleleSerotypeCount_" . $date . ".csv"); 
	print FILE "SERO,FULL,SEROTYPE,S,SC,InSilico\n";
	my $full = 0;
	my $sero = 0;
	my $short = 0;
	my $sc = 0;
	my $ins = 0;
	my $who_full = 0;
	my $who_sero = 0;
	my $who_short = 0;
	my $who_sc = 0;
	my $who_ins = 0;
	foreach my $type ( sort keys %full ) {
		print FILE $type . ",";
		print FILE $full_ref->{ $type }->[0] . ",";		# FULL
		print FILE $full_ref->{ $type }->[1] . ",";		# Serotype
		print FILE $full_ref->{ $type }->[2] . ",";		# S
		print FILE $full_ref->{ $type }->[3] . ",";		# SC
		print FILE $full_ref->{ $type }->[4] . "\n";		# InSilico

		$full = $full + $full_ref->{ $type }->[0];
		$sero = $sero + $full_ref->{ $type }->[1];
		$short = $short + $full_ref->{ $type }->[2];
		$sc = $sc + $full_ref->{ $type }->[3];
		$ins = $ins + $full_ref->{ $type }->[4];

		if ( exists $who_ref->{ $type } ) {
			$who_full = $who_full + $who_ref->{ $type }->[0];
			$who_sero = $who_sero + $who_ref->{ $type }->[1];
			$who_short = $who_short + $who_ref->{ $type }->[2];
			$who_sc = $who_sc + $who_ref->{ $type }->[3];
			$who_ins = $who_ins + $who_ref->{ $type }->[4];
		}
	}
	print FILE "\n";

	print FILE "AssignedTotal," . $full . "," . $sero . "," . $short . "," . $sc . "," . $ins . "\n";
	print FILE "WHOAcceptedTotal," . $who_full . "," . $who_sero . "," . $who_short . "," . $who_sc . "," . $who_ins . "\n\n";
	print FILE "CommonTotal," . $common_full . "," . $common_sero . "," . $common_short . "," . $common_sc . "," . $common_ins . "\n";
	print FILE "IntermediateTotal," . $inter_full . "," . $inter_sero . "," . $inter_short . "," . $inter_sc . "," . $inter_ins . "\n";
	print FILE "WellDocumentedTotal," . $wd_full . "," . $wd_sero . "," . $wd_short . "," . $wd_sc . "," . $wd_ins . "\n";

	my $rare_full = $full - $common_full - $inter_full - $wd_full;
	my $rare_sero = $sero - $common_sero - $inter_sero - $wd_sero;
	my $rare_short = $short - $common_short - $inter_short - $wd_short;
	my $rare_sc = $sc - $common_sc - $inter_sc - $wd_sc;
	my $rare_ins = $ins - $common_ins - $inter_ins - $wd_ins; 
	print FILE "RareTotal," . $rare_full . "," . $rare_sero . "," . $rare_short . "," . $rare_sc . "," . $rare_ins . "\n\n";
	
	my $null_count = scalar keys %$null_ref;
	my $qallele_count = scalar keys %$qallele_ref;
	my $subtotal = $full + $sero + $short + $sc + $ins;
	print FILE "TOTAL EXPRESSED ALLELES," . $subtotal . "\n";
	print FILE "Null alleles," . $null_count . "\n";
	print FILE "Q alleles," . $qallele_count . "\n";
	my $total = $full + $sero + $short + $sc + $ins + $null_count + $qallele_count;
	print FILE "TOTAL ALLELES," . $total . "\n";

	close FILE;
}


sub SUMMARY_TWO {
	my ($csv_ref, $gene, $sero_ref, $null_ref, $qallele_ref, $basetype_ref, $whotype_ref) = @_;
	my @list = Openfile::open_file_from_list($csv_ref);
	my $header = shift @list;
	my $common_full = 0;
	my $inter_full = 0;
	my $wd_full = 0;
		
	my $common_sero = 0;
	my $inter_sero = 0;
	my $wd_sero = 0;
		
	my $common_short = 0;
	my $inter_short = 0;
	my $wd_short = 0;
		
	my $common_sc = 0;
	my $inter_sc = 0;
	my $wd_sc = 0;
		
	my $common_ins = 0;
	my $inter_ins = 0;
	my $wd_ins = 0;
		
	my %full;
	my $full_ref = \%full;
	$full_ref->{ "None" }->[0] = 0;		# FULL
	$full_ref->{ "None" }->[1] = 0;		# Serotype
	$full_ref->{ "None" }->[2] = 0;		# S
	$full_ref->{ "None" }->[3] = 0;		# SC
	$full_ref->{ "None" }->[4] = 0;		# InSilico

	my %who;
	my $who_ref = \%who;

	my @serotype;
	my %unique;
	my $antigen_ref = ASSIGNED_SHORT::ANTIGEN();
	foreach my $type ( @{$sero_ref->[0]} ) {	# go through serotype
		my $official = "";
		if ( exists $antigen_ref->{ $type } ) {
			$official = $antigen_ref->{ $type };	# convert to official serotype, B-1501 => B62
		}
		else {
			$official = $type;
		}
		$type = $official;	#reassign $type value
		unless ( exists $unique{ $type } ) {
			push @serotype, $type;
			$unique{ $type } = 0;
		}
	}

	foreach my $type ( @serotype ) {	# go through serotype
		unless ( exists $full_ref->{ $type } ) {
			$full_ref->{ $type }->[0] = 0;		# FULL
			$full_ref->{ $type }->[1] = 0;		# Serotype
			$full_ref->{ $type }->[2] = 0;		# S
			$full_ref->{ $type }->[3] = 0;		# SC
			$full_ref->{ $type }->[4] = 0;		# InSilico

			$who_ref->{ $type }->[0] = 0;		# FULL
			$who_ref->{ $type }->[1] = 0;		# Serotype
			$who_ref->{ $type }->[2] = 0;		# S
			$who_ref->{ $type }->[3] = 0;		# SC
			$who_ref->{ $type }->[4] = 0;		# InSilico
		}
	}

	my %two;
	foreach my $line ( @list ) {
		my @elements = split( ",", $line );
		if ( $elements[1] eq "UNA" ) {
			if ( $elements[0] =~ /($gene\*\d+:\d+)/ ) {
				my $two_allele = $1;
				if ( exists $two{ $two_allele } ) {
					next;
				}
				else {
					$full_ref->{ "None" }->[0] = $full_ref->{ "None" }->[0] + 1;
					$two{ $two_allele } = 0;
					next;
				}
			}
		}

		if (( exists $null_ref->{ $elements[0] } ) || ( exists $qallele_ref->{ $elements[0] } )) {
			next;
		}

		if ( $elements[0] =~ /($gene\*\d+:\d+)/ ) {
			my $two_allele = $1;
			if ( exists $two{ $two_allele } ) {
				next;
			}
			else {
				foreach my $type ( @serotype ) {	# go through serotype
					if ( $elements[2] eq $type ) {
						if ( $elements[1] eq "FULL" ) {
							$full_ref->{ $type }->[0] = $full_ref->{ $type }->[0] + 1;
							if ( exists $whotype_ref->{$type} ) {	
								if ( grep( /^$whotype_ref->{$type}$/, @$basetype_ref )) {	#literal value lookup
									$who_ref->{ $type }->[0] = $who_ref->{ $type }->[0] + 1;
								}
							}
							# CIWD
							if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
								$common_full++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
								$inter_full++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
								$wd_full++;
							}
						}
						elsif ( $elements[1] eq "SEROTYPE" ) {
							$full_ref->{ $type }->[1] = $full_ref->{ $type }->[1] + 1;
							if ( exists $whotype_ref->{$type} ) {	
								if ( grep( /^$whotype_ref->{$type}$/, @$basetype_ref )) {	#literal value lookup
									$who_ref->{ $type }->[1] = $who_ref->{ $type }->[1] + 1;
								}
							}
							if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
								$common_sero++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
								$inter_sero++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
								$wd_sero++;
							}
						}
						elsif ( $elements[1] eq "S" ) {
							$full_ref->{ $type }->[2] = $full_ref->{ $type }->[2] + 1;
							if ( exists $whotype_ref->{$type} ) {	
								if ( grep( /^$whotype_ref->{$type}$/, @$basetype_ref )) {	#literal value lookup
									$who_ref->{ $type }->[2] = $who_ref->{ $type }->[2] + 1;
								}
							}
							if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
								$common_short++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
								$inter_short++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
								$wd_short++;
							}
						}
						elsif ( $elements[1] eq "SC" ) {
							$full_ref->{ $type }->[3] = $full_ref->{ $type }->[3] + 1;
							if ( exists $whotype_ref->{$type} ) {	
								if ( grep( /^$whotype_ref->{$type}$/, @$basetype_ref )) {	#literal value lookup
									$who_ref->{ $type }->[3] = $who_ref->{ $type }->[3] + 1;
								}
							}
							if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
								$common_sc++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
								$inter_sc++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
								$wd_sc++;
							}
						}
						else {
							$full_ref->{ $type }->[4] = $full_ref->{ $type }->[4] + 1;
							if ( exists $whotype_ref->{$type} ) {	
								if ( grep( /^$whotype_ref->{$type}$/, @$basetype_ref )) {	#literal value lookup
									$who_ref->{ $type }->[4] = $who_ref->{ $type }->[4] + 1;
								}
							}
							if (( exists $elements[5] ) && ( $elements[5] eq "C" )) {
								$common_ins++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "I" )) {
								$inter_ins++;
							}
							elsif (( exists $elements[5] ) && ( $elements[5] eq "WD" )) {
								$wd_ins++;
							}
						}
					}
				}
				$two{ $two_allele } = 0;
			}
		}

	}

	open(FILE, ">output/" . $gene . "_ProteinSerotypeCount_" . $date . ".csv"); 
	print FILE "SERO,FULL,SEROTYPE,S,SC,InSilico\n";
	my $full = 0;
	my $sero = 0;
	my $short = 0;
	my $sc = 0;
	my $ins = 0;
	my $who_full = 0;
	my $who_sero = 0;
	my $who_short = 0;
	my $who_sc = 0;
	my $who_ins = 0;
	foreach my $type ( sort keys %full ) {
		unless ( $type eq "None" ) {
			print FILE $type . ",";
			print FILE $full_ref->{ $type }->[0] . ",";		# FULL
			print FILE $full_ref->{ $type }->[1] . ",";		# Serotype
			print FILE $full_ref->{ $type }->[2] . ",";		# S
			print FILE $full_ref->{ $type }->[3] . ",";		# SC
			print FILE $full_ref->{ $type }->[4] . "\n";		# InSilico
			$full = $full + $full_ref->{ $type }->[0];
			$sero = $sero + $full_ref->{ $type }->[1];
			$short = $short + $full_ref->{ $type }->[2];
			$sc = $sc + $full_ref->{ $type }->[3];
			$ins = $ins + $full_ref->{ $type }->[4];
			if ( exists $who_ref->{ $type } ) {
				$who_full = $who_full + $who_ref->{ $type }->[0];
				$who_sero = $who_sero + $who_ref->{ $type }->[1];
				$who_short = $who_short + $who_ref->{ $type }->[2];
				$who_sc = $who_sc + $who_ref->{ $type }->[3];
				$who_ins = $who_ins + $who_ref->{ $type }->[4];
			}
		}
	}
	
	print FILE "None,";
	print FILE $full_ref->{ "None" }->[0] . ",";		# FULL
	print FILE $full_ref->{ "None" }->[1] . ",";		# Serotype
	print FILE $full_ref->{ "None" }->[2] . ",";		# S
	print FILE $full_ref->{ "None" }->[3] . ",";		# SC
	print FILE $full_ref->{ "None" }->[4] . "\n\n";		# InSilico
	print FILE "AssignedTotal," . $full . "," . $sero . "," . $short . "," . $sc . "," . $ins . "\n";
	print FILE "WHOAcceptedTotal," . $who_full . "," . $who_sero . "," . $who_short . "," . $who_sc . "," . $who_ins . "\n\n";
	print FILE "CommonTotal," . $common_full . "," . $common_sero . "," . $common_short . "," . $common_sc . "," . $common_ins . "\n";
	print FILE "IntermediateTotal," . $inter_full . "," . $inter_sero . "," . $inter_short . "," . $inter_sc . "," . $inter_ins . "\n";
	print FILE "Well-DocumentedTotal," . $wd_full . "," . $wd_sero . "," . $wd_short . "," . $wd_sc . "," . $wd_ins . "\n";

	my $rare_full = $full - $common_full - $inter_full - $wd_full;
	my $rare_sero = $sero - $common_sero - $inter_sero - $wd_sero;
	my $rare_short = $short - $common_short - $inter_short - $wd_short;
	my $rare_sc = $sc - $common_sc - $inter_sc - $wd_sc;
	my $rare_ins = $ins - $common_ins - $inter_ins - $wd_ins; 
	print FILE "RareTotal," . $rare_full . "," . $rare_sero . "," . $rare_short . "," . $rare_sc . "," . $rare_ins . "\n\n";
	
	my $null_count = scalar keys %$null_ref;
	my $qallele_count = scalar keys %$qallele_ref;
	my $subtotal = $full + $sero + $short + $sc + $ins + $full_ref->{ "None" }->[0];
	print FILE "TOTAL EXPRESSED TWO-FILED ALLELES," . $subtotal . "\n";
	print FILE "Null alleles," . $null_count . "\n";
	print FILE "Q alleles," . $qallele_count . "\n";
	my $total = $full + $sero + $short + $sc + $ins + $full_ref->{ "None" }->[0] + $null_count + $qallele_count;
	print FILE "TOTAL TWO-FILED ALLELES," . $total . "\n";

	close FILE;
}

1;
