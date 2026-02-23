#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: COUNT.pm 
# This module was developed to generate Summary table
# Deleted SUMMARY and SUMMARY_TWO, and moved them to SUMCOUNT.pm
# last modified and documented on February 20 2026

package COUNT;
use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
use Openfile;
use GROUP_SORT;
use ASSIGNED_SHORT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
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
								#print $type . "\n";
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
								#	print $type . "\n";
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


1;
