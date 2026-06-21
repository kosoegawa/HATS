#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: combine.pl
# This is to merge all serotype tables in RESULTS directory
# last modified and documented on April 5 2026

use strict;
use lib 'SEROTYPE';
use Openfile;
use GROUP_SORT;
use RESULT_COMBINE;
use PRACTICAL_COMBINE;


my $output = "output/";
my $combined_dir = $output . "COMBINED/";
my @csv = glob($combined_dir . "*.csv");
foreach my $csv ( @csv ) {
	unlink $csv;
}

my $results_dir = $output . "RESULTS/";
RESULT_COMBINE::COMBINE( $results_dir, $combined_dir );

$results_dir = $output . "TWORESULTS/";
RESULT_COMBINE::COMBINE( $results_dir, $combined_dir  );

$results_dir = $output . "PRACTICAL/";
PRACTICAL_COMBINE::COMBINE( $results_dir, $combined_dir  );

$results_dir = $output . "TWOPRACTICAL/";
PRACTICAL_COMBINE::COMBINE( $results_dir, $combined_dir  );



