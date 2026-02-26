#!/usr/bin/perl -w
#
# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: runDRB3.pl 
# Driver for HLA-DRB3
# If partial sequences are used as a reference, add the optional argument
# last modified and documented on February 26 2026

use strict;
use lib 'SEROTYPE';
use ORGANIZE;
use STRASSIGN;
use RESIDUES;
use NullAllele;
use QAllele;
use ASSIGN;
use DRB3_INFO;
use COUNT;
use ASSIGNED_SHORT;
use COPYRESULT;
use POSIX qw(strftime);
use SUMCOUNT;
use COMBINE;
use LEGACY_REPORT;
use PRACTICAL;

my $date = strftime "%Y-%m-%d", localtime;
chomp $date;    # remove newline character

#capture input file
my @file = glob('input/*');
my $database = "3.39.0";	# IPD-IMGT/HLA database version
my $hats = "HATSv3.0.0";	# HATS version
my $file = "";
foreach my $tmp ( @file ) {
	print $tmp . "\n";
	if ( $tmp =~ /hla_prot\.fasta\.(.*+)/ ) {
	# capture database version
		$database = $1;
		$file = $tmp;
	}
	elsif ( $tmp =~ /(HATSv.*)/ ) {
		$hats = $1;
	}
}	

#remove all csv files
my @csv = glob('output/*.csv');
my $csv_ref = \@csv;
my $csvs = scalar @csv;
if ( $csvs > 0 ) {
	unlink @csv;
}

open ( FILE, ">output/" . $hats . "_IMGT_" . $database . ".csv" );	#create an empty file to tag database version	
close FILE;

my $fasta_ref = ORGANIZE::fasta( $file );	# organize fasta

my $gene = DRB3_INFO::DRB3();
my $ciwd_ref = ORGANIZE::CIWD( $gene );
my $cwd_ref = ORGANIZE::CWD( $gene );
my $ecwd_ref = ORGANIZE::EURCWD( $gene );

my $leader = DRB3_INFO::DRB3_LEADER();
my $ref_ref = DRB3_INFO::REF("ALL");
my $residues_all_ref = DRB3_INFO::RESIDUES("ALL");
my $partial_ref = DRB3_INFO::PARTIAL();

my $group_ref = DRB3_INFO::GROUP();
my $base_ref = DRB3_INFO::BASE();
my $basetype_ref = DRB3_INFO::BASETYPE();

#print target residues
my $elements_ref = RESIDUES::pattern( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref, $partial_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );
#print relax target residues
RESIDUES::LAX( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref, $partial_ref, $base_ref, $group_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );

#print null alleles
my $null_ref = NullAllele::all( $fasta_ref, $gene );

#print Q alleles
my $qallele_ref = QAllele::all( $fasta_ref, $gene );

# Stringent condition
my $assigned_ref = STRASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_all_ref );

# capture group IDs
my @group;
my %element;
foreach my $element ( sort values %$group_ref ) {
	unless ( exists $element{ $element } ) {
		push @group, $element;
		$element{ $element } = 0;
	}
}

my $known_cross_ref = DRB3_INFO::KNOWN_CROSS();
my @known_cross = keys %$known_cross_ref;

# assign LAX condition
my %cross;
my $cross_ref;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = DRB3_INFO::REF( $group[ $index ] );
	foreach my $known ( @known_cross ) {
		if (exists $ref_ref->{ $known } ) {
			delete( $ref_ref->{ $known } );
			print $known . " deleted\n";
		}
	}
	my $residues_ref = DRB3_INFO::RESIDUES( $group[ $index ] );
	$cross_ref = ASSIGN::CROSS( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref );
	$assigned_ref = ASSIGN::ASSIGN($fasta_ref, $assigned_ref, $gene, $leader,$ref_ref, $residues_ref, $partial_ref, $known_cross_ref );
}

my $unassigned_ref = ASSIGN::UNASSIGNED( $fasta_ref, $assigned_ref, $gene );

@csv = glob('output/*.csv');
my $sero_ref = DRB3_INFO::SERO();
my $key_ref = DRB3_INFO::KEY();
#generate Summary table
COUNT::COUNT($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);
# generate two-field summary table
COUNT::TWOFIELD($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);

my %short;
my $short_ref = \%short;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = DRB3_INFO::REF( $group[ $index ] );
	my $residues_ref = DRB3_INFO::RESIDUES( $group[ $index ] );
	$short_ref = ASSIGN::SHORT($fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref  );
}

# generates residues for all two-field alleles
my $elements2_ref = RESIDUES::ELEMENTS ( $elements_ref,$fasta_ref,$gene,$null_ref,$qallele_ref,$residues_all_ref,$leader,$partial_ref,$assigned_ref,$short_ref );
ASSIGNED_SHORT::PRINT_RESIDUES( $elements2_ref,$gene,$residues_all_ref,$database );

# assign SHORT
ASSIGNED_SHORT::PRINT( $unassigned_ref, $short_ref );
# generate final table
my $broad_ref = DRB3_INFO::BROAD();

LEGACY_REPORT::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );
LEGACY_REPORT::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );
my $parent_ref = DRB3_INFO::PARENT();

COMBINE::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );
COMBINE::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );

PRACTICAL::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );
PRACTICAL::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$parent_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref );

@csv = glob("output/" . $gene . "_Allele_Antigen_Practical_Table_IMGT_HLA_*");
foreach my $csv ( @csv ) {
	my $out_name = "_Allele_Antigen_Count_";
	SUMCOUNT::SUMMARY($csv, $gene, $null_ref, $qallele_ref, $out_name, $residues_all_ref);
}

@csv = glob("output/" . $gene . "_Protein_Antigen_Practical_Table_IMGT_HLA_*");
foreach my $csv ( @csv ) {
	my $out_name = "_Protein_Antigen_Count_";
	SUMCOUNT::SUMMARY($csv, $gene, $null_ref, $qallele_ref, $out_name, $residues_all_ref);
}

COPYRESULT::COPYRESULT( $gene, $database, $date );
COPYRESULT::COPYTWORESULT( $gene, $database, $date );
COPYRESULT::COPYPRACTICAL( $gene, $database, $date );
COPYRESULT::COPYRESIDUE( $gene, $database, $date );
COPYRESULT::COPYLEGACY( $gene, $database, $date );


