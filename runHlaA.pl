#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: runHlaA.pl
# Driver for HLA-A
# last modified and documented on October 26 2023

use strict;
use lib 'SEROTYPE';
use ORGANIZE;
use STRASSIGN;
use RESIDUES;
use NullAllele;
use QAllele;
use ASSIGN;
use HLAA_INFO;
use COUNT;
use ASSIGNED_SHORT;
use ABw;
use Bw46ASSIGN;
use File::Copy;

my $date = `date +%F`;          # invoke bash date command
chomp $date;    # remove newline character

#capture input file
my @file = glob('input/hla_prot.fasta*');
my $file = "";
foreach my $tmp ( @file ) {
	$file = $tmp;
	print $file . "\n";
}	
# capture database version
my $database = "3.39.0";
if ( $file =~ /hla_prot\.fasta\.(.*+)/ ) {
	$database = $1;
}

#remove all csv files
my @csv = glob('output/*.csv');
my $csv_ref = \@csv;
my $csvs = scalar @csv;
if ( $csvs > 0 ) {
	unlink @csv;
}

my @fasta = glob('output/*.fasta');
my $fasta_count = scalar @fasta;
if ( $fasta_count > 0 ) {
	unlink @fasta;
}

open ( FILE, ">output/" . $database . ".csv" );	#create an empty file to tag database version	
close FILE;

my $fasta_ref = ORGANIZE::fasta( $file );	# organize fasta
my $gene = HLAA_INFO::HLAA();	# capture gene
my $ciwd_ref = ORGANIZE::CIWD( $gene );		#CIWD3.0
my $cwd_ref = ORGANIZE::CWD( $gene );		#CWD2.0
my $ecwd_ref = ORGANIZE::EURCWD( $gene );	#EUR CWD

my $leader = HLAA_INFO::HLAA_LEADER();
my $ref_ref = HLAA_INFO::REF("ALL");
my $residues_all_ref = HLAA_INFO::RESIDUES("ALL");
my $partial_ref = HLAA_INFO::PARTIAL();	# added this to handle partial sequence in a better way

my $group_ref = HLAA_INFO::GROUP();
my $base_ref = HLAA_INFO::BASE();
my $basetype_ref = HLAA_INFO::BASETYPE();

#print target residues
my $elements_ref = RESIDUES::pattern( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref, $partial_ref, $basetype_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );
#print relax target residues
RESIDUES::LAX( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref, $partial_ref, $basetype_ref, $base_ref, $group_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );

#print null alleles
my $null_ref = NullAllele::all( $fasta_ref, $gene );

#print Q alleles
my $qallele_ref = QAllele::all( $fasta_ref, $gene );

# partial sequences are present
# Stringent condition
# added $partial_ref as an argument to handle partial reference sequences
my $assigned_ref = STRASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref, $partial_ref );

# capture group IDs
my @group;
my %element;
foreach my $element ( sort values %$group_ref ) {
	unless ( exists $element{ $element } ) {
		push @group, $element;
		$element{ $element } = 0;
	}
}

my $known_cross_ref = HLAA_INFO::KNOWN_CROSS();
my @known_cross = keys %$known_cross_ref;

# assign LAX condition
my %cross;
my $cross_ref;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = HLAA_INFO::REF( $group[ $index ] );
	foreach my $known ( @known_cross ) {
		if (exists $ref_ref->{ $known } ) {
			delete( $ref_ref->{ $known } );
		}
	}
	my $residues_ref = HLAA_INFO::RESIDUES( $group[ $index ] );
	$cross_ref = ASSIGN::CROSS( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref );
	$assigned_ref = ASSIGN::ASSIGN($fasta_ref, $assigned_ref, $gene, $leader,$ref_ref, $residues_ref, $partial_ref, $known_cross_ref );
}

my $unassigned_ref = ASSIGN::UNASSIGNED( $fasta_ref, $assigned_ref, $gene );

@csv = glob('output/*.csv');
my $sero_ref = HLAA_INFO::SERO();
my $key_ref = HLAA_INFO::KEY();
#generate Summary table
COUNT::COUNT($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);
# generate two-field summary table
COUNT::TWOFIELD($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);

my %short;
my $short_ref = \%short;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = HLAA_INFO::REF( $group[ $index ] );
	my $residues_ref = HLAA_INFO::RESIDUES( $group[ $index ] );
	$short_ref = ASSIGN::SHORT($fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref  );
}

# generates residues for all two-field alleles
my $elements2_ref = RESIDUES::ELEMENTS ( $elements_ref,$fasta_ref,$gene,$null_ref,$qallele_ref,$residues_all_ref,$leader,$partial_ref,$assigned_ref,$short_ref );
ASSIGNED_SHORT::PRINT_RESIDUES( $elements2_ref,$gene,$residues_all_ref,$database );

# assign SHORT
ASSIGNED_SHORT::PRINT( $unassigned_ref, $short_ref );

# generate final table
my $bw_ref = HLAA_INFO::BW();
my $broad_ref = HLAA_INFO::BROAD();

$ref_ref = ABw::REF();
my $residues_ref = ABw::RESIDUES();
my $bw = ABw::BW();
my $bw_ref2 = Bw46ASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $bw );

ASSIGNED_SHORT::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2 );
ASSIGNED_SHORT::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2 );

@csv = glob("output/" . $gene . "_Serotype_Table_IMGT_HLA_*");
foreach my $csv ( @csv ) {
	my $whotype_ref = HLAA_INFO::WHO();
	COUNT::SUMMARY($csv, $gene, $null_ref, $qallele_ref, $whotype_ref);
	COUNT::SUMMARY_TWO($csv, $gene, $sero_ref, $null_ref, $qallele_ref, $basetype_ref, $whotype_ref);
}

@csv = glob("RESULTS/" . $gene . "_Serotype_Table_IMGT_HLA_*");
$csvs = scalar @csv;
if ( $csvs > 0 ) {
	unlink @csv;
}

copy("output/" . $gene . "_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "RESULTS/") or die "Copy failed: $!";
copy("output/" . $gene . "_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "TWORESULTS/") or die "Copy failed: $!";
