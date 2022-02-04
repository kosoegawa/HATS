#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: runDRB1.pl 
# Driver for HLA-DRB1
# If partial sequences are used as a reference, add the optional argument
# last modified and documented on February 3 2022

use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
use ORGANIZE;
use STRASSIGN;
use RESIDUES;
use NullAllele;
use QAllele;
use ASSIGN;
use DRB1_INFO;
use COUNT;
use ASSIGNED_SHORT;
use DR5231;
use DR5231ASSIGN;
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

open ( FILE, ">output/" . $database . ".csv" );	#create an empty file to tage database version	
close FILE;

my $fasta_ref = ORGANIZE::fasta( $file );	# organize fasta

my $gene = DRB1_INFO::DRB1();
my $ciwd_ref = ORGANIZE::CIWD( $gene );
my $cwd_ref = ORGANIZE::CWD( $gene );
my $ecwd_ref = ORGANIZE::EURCWD( $gene );

my $leader = DRB1_INFO::DRB1_LEADER();
my $ref_ref = DRB1_INFO::REF("ALL");
my $residues_ref = DRB1_INFO::RESIDUES("ALL");
my $partial_ref = DRB1_INFO::PARTIAL();

my $group_ref = DRB1_INFO::GROUP();
my $base_ref = DRB1_INFO::BASE();
my $basetype_ref = DRB1_INFO::BASETYPE();

#print target residues
RESIDUES::pattern( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $basetype_ref, $base_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );
#print relax target residues
RESIDUES::LAX( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $basetype_ref, $base_ref, $group_ref, $ciwd_ref, $cwd_ref, $ecwd_ref );

#print null alleles
my $null_ref = NullAllele::all( $fasta_ref, $gene );

#print Q alleles
my $qallele_ref = QAllele::all( $fasta_ref, $gene );

# partial sequences are present
# Stringent condition
# added $partial_ref as an argument to handle partial reference sequences
my $assigned_ref = STRASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref );

# capture group IDs
my @group;
my %element;
foreach my $element ( sort values %$group_ref ) {
	unless ( exists $element{ $element } ) {
		push @group, $element;
		$element{ $element } = 0;
	}
}
# assign LAX condition
my %cross;
my $cross_ref;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = DRB1_INFO::REF( $group[ $index ] );
	$residues_ref = DRB1_INFO::RESIDUES( $group[ $index ] );
	$cross_ref = ASSIGN::CROSS( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref );
	$assigned_ref = ASSIGN::ASSIGN($fasta_ref, $assigned_ref, $gene, $leader,$ref_ref, $residues_ref, $partial_ref );
}

my $unassigned_ref = ASSIGN::UNASSIGNED( $fasta_ref, $assigned_ref, $gene );

@csv = glob('output/*.csv');
my $sero_ref = DRB1_INFO::SERO();
my $key_ref = DRB1_INFO::KEY();
#generate Summary table
COUNT::COUNT($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);
# generate two-field summary table
COUNT::TWOFIELD($csv_ref, $gene, $sero_ref, $key_ref, $null_ref, $base_ref, $basetype_ref, $qallele_ref);

my %short;
my $short_ref = \%short;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = DRB1_INFO::REF( $group[ $index ] );
	$residues_ref = DRB1_INFO::RESIDUES( $group[ $index ] );
	$short_ref = ASSIGN::SHORT($fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref  );
}

# assign SHORT
ASSIGNED_SHORT::PRINT( $unassigned_ref, $short_ref );

# generate final table
my $dr5231_ref = DRB1_INFO::DR5231();	# this is empty place holder
my $broad_ref = DRB1_INFO::BROAD();

$ref_ref = DR5231::REF();
my $dr5231 = DR5231::DR5231();
my %dr5231;
my $dr5231_ref2 = \%dr5231;
foreach my $key ( keys %$ref_ref ) {
	$residues_ref = DR5231::RESIDUES( $key );
	my $tmp_ref =  DR5231ASSIGN::all( $fasta_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $dr5231 );
	%dr5231 = ( %dr5231, %$tmp_ref );
}

#ASSIGNED_SHORT::AC_COMBINED( $database,$null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
ASSIGNED_SHORT::COMBINED( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$dr5231_ref,$dr5231_ref2 );
ASSIGNED_SHORT::COMBINED_TWO( $database, $null_ref,$qallele_ref,$assigned_ref,$unassigned_ref,$short_ref,$gene,$base_ref,$basetype_ref,$cross_ref,
$broad_ref,$ciwd_ref,$cwd_ref,$ecwd_ref,$dr5231_ref,$dr5231_ref2 );

@csv = glob("output/" . $gene . "_Serotype_Table_IMGT_HLA_*");
foreach my $csv ( @csv ) {
	my $whotype_ref = DRB1_INFO::WHO();
	COUNT::SUMMARY($csv, $gene, $sero_ref, $null_ref, $qallele_ref, $basetype_ref, $whotype_ref);
	COUNT::SUMMARY_TWO($csv, $gene, $sero_ref, $null_ref, $qallele_ref, $basetype_ref, $whotype_ref);
}

@csv = glob("RESULTS/" . $gene . "_Serotype_Table_IMGT_HLA_*");
$csvs = scalar @csv;
if ( $csvs > 0 ) {
	unlink @csv;
}

copy("output/" . $gene . "_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "RESULTS/") or die "Copy failed: $!";
copy("output/" . $gene . "_TwoField_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "TWORESULTS/") or die "Copy failed: $!";
