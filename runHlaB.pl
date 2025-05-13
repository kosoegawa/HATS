#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: runHlaB.pl 
# Driver for HLA-B
# If partial sequences are used as a reference, add the optional argument
# last modified and documented on May 13 2025

use strict;
use lib 'SEROTYPE';
use ORGANIZE;
use STRASSIGN;
use RESIDUES;
use NullAllele;
use QAllele;
use ASSIGN;
use HLAB_INFO;
use COUNT;
use ASSIGNED_SHORT;
use BBw;
use BC1C2;
use Bw46ASSIGN;
use COPYRESULT;
use POSIX qw(strftime);

my $date = strftime "%Y-%m-%d", localtime;
#my $date = `date +%F`;          # invoke bash date command
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

my @fasta = glob('output/*.fasta');
my $fasta_count = scalar @fasta;
if ( $fasta_count > 0 ) {
	unlink @fasta;
}

open ( FILE, ">output/" . $hats . "_IMGT_" . $database );	#create an empty file to tag database version	
close FILE;

my $fasta_ref = ORGANIZE::fasta( $file );	# organize fasta
my $gene = HLAB_INFO::HLAB();
my $ciwd_ref = ORGANIZE::CIWD( $gene );
my $cwd_ref = ORGANIZE::CWD( $gene );
my $ecwd_ref = ORGANIZE::EURCWD( $gene );

my $leader = HLAB_INFO::HLAB_LEADER();
my $ref_ref = HLAB_INFO::REF("ALL");
my $residues_all_ref = HLAB_INFO::RESIDUES("ALL");
my $partial_ref = HLAB_INFO::PARTIAL();	# added this to handle partial sequence in a better way

my $group_ref = HLAB_INFO::GROUP();
my $base_ref = HLAB_INFO::BASE();
my $basetype_ref = HLAB_INFO::BASETYPE();

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
foreach my $element ( sort values %$group_ref ) {	#values
	unless ( exists $element{ $element } ) {
		push @group, $element;
		$element{ $element } = 0;
	}
}

my $known_cross_ref = HLAB_INFO::KNOWN_CROSS();
my @known_cross = keys %$known_cross_ref;

# assign LAX condition
my %cross;
my $cross_ref;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = HLAB_INFO::REF( $group[ $index ] );
	foreach my $known ( @known_cross ) {
		if (exists $ref_ref->{ $known } ) {
			delete( $ref_ref->{ $known } );
		}
	}
	my $residues_ref = HLAB_INFO::RESIDUES( $group[ $index ] );
	$cross_ref = ASSIGN::CROSS( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref );
	$assigned_ref = ASSIGN::ASSIGN($fasta_ref, $assigned_ref, $gene, $leader,$ref_ref, $residues_ref, $partial_ref, $known_cross_ref );
}

my $unassigned_ref = ASSIGN::UNASSIGNED( $fasta_ref, $assigned_ref, $gene );

@csv = glob('output/*.csv');
my $sero_ref = HLAB_INFO::SERO();
my $key_ref = HLAB_INFO::KEY();
#generate Summary table
COUNT::COUNT($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);
# generate two-field summary table
COUNT::TWOFIELD($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);

my %short;
my $short_ref = \%short;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	my $ref_ref = HLAB_INFO::REF( $group[ $index ] );
	my $residues_ref = HLAB_INFO::RESIDUES( $group[ $index ] );
	$short_ref = ASSIGN::SHORT($fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref  );
}

my $elements2_ref = RESIDUES::ELEMENTS ( $elements_ref,$fasta_ref,$gene,$null_ref,$qallele_ref,$residues_all_ref,$leader,$partial_ref,$assigned_ref,$short_ref );
ASSIGNED_SHORT::PRINT_RESIDUES( $elements2_ref,$gene,$residues_all_ref,$database );

# assign SHORT
ASSIGNED_SHORT::PRINT( $unassigned_ref, $short_ref );

# generate final table
my $bw_ref = HLAB_INFO::BW();
my $broad_ref = HLAB_INFO::BROAD();

$ref_ref = BBw::REF();
my $residues_ref = BBw::RESIDUES();
my $bw = BBw::BW();
my $bw_ref2 = Bw46ASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $bw );

$ref_ref = BC1C2::REF();
$residues_ref = BC1C2::RESIDUES();
my $c1c2 = BC1C2::C1C2();
my $c1c2_ref = Bw46ASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $c1c2 );

ASSIGNED_SHORT::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2 );
ASSIGNED_SHORT::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$bw_ref,$bw_ref2,$c1c2_ref );

@csv = glob("output/" . $gene . "_Serotype_Table_IMGT_HLA_*");
foreach my $csv ( @csv ) {
	my $whotype_ref = HLAB_INFO::WHO();
	COUNT::SUMMARY($csv, $gene, $null_ref, $qallele_ref, $whotype_ref);
	COUNT::SUMMARY_TWO($csv, $gene, $sero_ref, $null_ref, $qallele_ref, $basetype_ref, $whotype_ref);
}

COPYRESULT::COPYRESULT( $gene, $database, $date );
COPYRESULT::COPYTWORESULT( $gene, $database, $date );
COPYRESULT::COPYRESIDUE( $gene, $database, $date );

