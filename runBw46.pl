#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module:runBw46.pl 
# Driver for Bw4 and Bw6
# last modified and documented on June 19 2026

use strict;
use lib 'SEROTYPE';
use HLAA_INFO;
use HLAB_INFO;
use CBw;
use POSIX qw(strftime);
use COMBINE;
use HATS_VERSION;


my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

my $output = "output/";
my $combined = "COMBINED/";

my @bw4bw6 = glob('COMBINED/Antigen_Bw4Bw6_*');
my $bw4bw6_count = scalar @bw4bw6;
if ( $bw4bw6_count > 0 ) {
	unlink @bw4bw6;
}

#capture input file
my @file = glob('input/*');
my $database =  HATS_VERSION::IMGT_HLA_VERSION();	# IPD-IMGT/HLA database version
my $hats = HATS_VERSION::VERSION();	# HATS version
my $file = "";
foreach my $tmp ( @file ) {
	print $tmp . "\n";
	if ( $tmp =~ /hla_prot\.fasta/ ) {
		$file = $tmp;
	}
}	

my $a_bw_ref = HLAA_INFO::BW();
my $b_bw_ref = HLAB_INFO::BW();
my $c_bw_ref = CBw::BW();
my $antigen_ref = COMBINE::ANTIGEN();

open ( FILE, ">" . $output . $combined . "Antigen_Bw4Bw6_" . $database . "_" . $date . ".csv" );	#create an empty file to tage database version	
print FILE "B Antigen,Bw4/Bw6 Inclusion\n";
foreach my $antigen ( sort keys %$b_bw_ref ) {
	if ( exists $antigen_ref->{$antigen} ) {
		print FILE $antigen_ref->{$antigen} . "," . $b_bw_ref->{$antigen} . "\n";
	}
	else {
		print FILE $antigen . "," . $b_bw_ref->{$antigen} . "\n";
	}
}

print FILE "\nA Antigen,Bw4 Associated Reactivity\n";
foreach my $antigen ( sort keys %$a_bw_ref ) {
	if ( exists $antigen_ref->{$antigen} ) {
		print FILE $antigen_ref->{$antigen} . "," . $a_bw_ref->{$antigen} . "\n";
	}
	else {
		print FILE $antigen . "," . $a_bw_ref->{$antigen} . "\n";
	}
}

print FILE "\nCw Antigen,Potential Bw6 Associated Reactivity\n";
foreach my $antigen ( sort keys %$c_bw_ref ) {
	if ( $antigen eq "Cw1212" ) {
		if ( exists $antigen_ref->{$antigen} ) {
			print FILE $antigen_ref->{$antigen} . "," . $c_bw_ref->{$antigen} . "\n";
		}
		else {
			print FILE $antigen . "," . $c_bw_ref->{$antigen} . "\n";
		}
	}
}


close FILE;


