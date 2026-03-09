#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: runHATSversion.pl
# last modified and documented on March 8 2026

use strict;
use lib 'SEROTYPE';
use POSIX qw(strftime);
use HATS_VERSION;

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

#capture input file
my @file = glob('input/*');
my $database = "3.39.0";	# IPD-IMGT/HLA database version
my $hats = HATS_VERSION::VERSION();	# HATS version
my $file = "";
foreach my $tmp ( @file ) {
	print $tmp . "\n";
	if ( $tmp =~ /hla_prot\.fasta\.(.*+)/ ) {
	# capture database version
		$database = $1;
		$file = $tmp;
	}
}	

open ( FILE, ">COMBINED/" . $hats . "_IMGT_" . $database . ".csv" );	#create an empty file to tag database version	
close FILE;

