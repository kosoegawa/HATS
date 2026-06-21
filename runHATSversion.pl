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

my $output = "output/";
my $combined = "COMBINED/";

my $database =  HATS_VERSION::IMGT_HLA_VERSION();	# IPD-IMGT/HLA database version
my $hats = HATS_VERSION::VERSION();	# HATS version

open ( FILE, ">" . $output . $combined . $hats . "_IMGT_" . $database . ".csv" );	#create an empty file to tag database version	
close FILE;

