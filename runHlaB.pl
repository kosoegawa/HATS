#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: runHlaB.pl 
# Driver for HLA-B
# If partial sequences are used as a reference, add the optional argument
# last modified and documented on October 18 2022

use strict;
use lib '/data/kazu/workplace/serotype/SEROTYPE';
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

my $gene = HLAB_INFO::HLAB();
my $ciwd_ref = ORGANIZE::CIWD( $gene );
my $cwd_ref = ORGANIZE::CWD( $gene );
my $ecwd_ref = ORGANIZE::EURCWD( $gene );

my $leader = HLAB_INFO::HLAB_LEADER();
my $ref_ref = HLAB_INFO::REF("ALL");
my $residues_ref = HLAB_INFO::RESIDUES("ALL");
my $partial_ref = HLAB_INFO::PARTIAL();	# added this to handle partial sequence in a better way

my $group_ref = HLAB_INFO::GROUP();
my $base_ref = HLAB_INFO::BASE();
my $basetype_ref = HLAB_INFO::BASETYPE();

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
foreach my $element ( sort values %$group_ref ) {	#values
	unless ( exists $element{ $element } ) {
		push @group, $element;
		$element{ $element } = 0;
	}
}
# assign LAX condition
my %cross;
my $cross_ref;
for ( my $index = 0; $index < scalar @group; $index++ ) {
	$ref_ref = HLAB_INFO::REF( $group[ $index ] );
	$residues_ref = HLAB_INFO::RESIDUES( $group[ $index ] );
	$cross_ref = ASSIGN::CROSS( $fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $partial_ref, $cross_ref );
	$assigned_ref = ASSIGN::ASSIGN($fasta_ref, $assigned_ref, $gene, $leader,$ref_ref, $residues_ref, $partial_ref );
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
	$ref_ref = HLAB_INFO::REF( $group[ $index ] );
	$residues_ref = HLAB_INFO::RESIDUES( $group[ $index ] );
	$short_ref = ASSIGN::SHORT($fasta_ref, $assigned_ref, $gene, $leader, $ref_ref, $residues_ref, $short_ref, $partial_ref  );
}

# assign SHORT
ASSIGNED_SHORT::PRINT( $unassigned_ref, $short_ref );

# generate final table
my $bw_ref = HLAB_INFO::BW();
my $broad_ref = HLAB_INFO::BROAD();

$ref_ref = BBw::REF();
$residues_ref = BBw::RESIDUES();
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

@csv = glob("RESULTS/" . $gene . "_Serotype_Table_IMGT_HLA_*");
$csvs = scalar @csv;
if ( $csvs > 0 ) {
	unlink @csv;
}

copy("output/" . $gene . "_Serotype_Table_IMGT_HLA_" . $database . "_" . $date . ".csv", "RESULTS/") or die "Copy failed: $!";

