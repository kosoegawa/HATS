#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: DRB1_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last reviewed, modified and documented on April 28 2025

package DRB1_INFO;
use strict;

my @dr1 = (9, 10, 11, 12, 13, 70, 71);	#70,71 D,E => 103; changed 
my @dr3 = (9, 10, 11, 12, 13, 47, 58, 71, 74);
my @dr4 = (9, 10, 11, 12, 13, 70, 71, 74);
my @dr7 = (9, 10, 11, 12, 13, 60, 74);	# added residue 60S that is immunogenic  14,
my @dr8 = (9, 10, 11, 12, 13, 16, 58, 60, 67, 71, 74);
my @dr9 = (9, 10, 11, 12, 13, 60, 74);	#60S is included. See DR7
my @dr10 = (9, 10, 11, 12, 13);
my @dr11 = (9, 10, 11, 12, 13, 58, 67, 70, 71);	# 70, 71 were added to distinguish DR1102
my @dr12 = (9, 10, 11, 12, 13, 58, 60, 67, 71, 74);	# residue 71 is required to distinguish from DR1317, removed 47
my @dr13 = (9, 10, 11, 12, 13, 58, 60, 67, 70, 71);	# residue 58 is required to distinguish from DR11 & DR1102
my @dr14 = (9, 10, 11, 12, 13, 16, 58, 60, 67,70, 71, 74);
my @dr15 = (9, 10, 11, 12, 13, 67, 70);
my @dr16 = (9, 10, 11, 12, 13, 67, 70);

my %dr5231;
my %dr1;
my %group;
my %base;
my @subtype;
$dr1{"DR-0101"} = "HLA00664";		# DRB1*01:01:01:01
$dr1{"DR-0103"} = "HLA00667";		# DRB1*01:03:01
$group{"DR-0101"} = "DR1"; $group{"DR-0103"} = "DR1";
$base{"DR-0101"} = "DR1"; $base{"DR-0103"} = "DR103";
push @subtype, ("DR-103");

my %dr3;
$dr3{"DR-0301"} = "HLA00671";		# DRB1*03:01:01:01: DR17
$dr3{"DR-0302"} = "HLA00673";		# DRB1*03:02:01: DR18
$group{"DR-0301"} = "DR3"; $group{"DR-0302"} = "DR3";
$base{"DR-0301"} = "DR17"; $base{"DR-0302"} = "DR18";
push @subtype, ("DR-0302");

my %dr4;
$dr4{"DR-0401"} = "HLA00685";		# DRB1*04:01:01:01
$dr4{"DR-0404"} = "HLA00689";		# DRB1*04:04:01
$dr4{"DR-0402"} = "HLA00687";		# DRB1*04:02:01
$dr4{"DR-0403"} = "HLA00688";		# DRB1*04:03:01:01
$dr4{"DR-0412"} = "HLA00698";		# DRB1*04:12
$dr4{"DR-0415"} ="HLA00701";		# DRB1*04:15
$dr4{"DR-1410"} = "HLA00842";		# DRB1*14:10
$group{"DR-0401"} = "DR4"; $group{"DR-0404"} = "DR4"; $group{"DR-0402"} = "DR4"; $group{"DR-0403"} = "DR4"; $group{"DR-0412"} = "DR4"; $group{"DR-0415"} = "DR4";
$group{"DR-1410"} = "DR4";
$base{"DR-0401"} = "DR4"; $base{"DR-0404"} = "DR4"; $base{"DR-0402"} = "DR4"; $base{"DR-0403"} = "DR4"; $base{"DR-0412"} = "DR4"; $base{"DR-0415"} = "DR4";
$base{"DR-1410"} = "DR4";
push @subtype, ("DR-0404","DR-0402","DR403","DR-0412","DR-0415","DR-1410");

my %dr7;
$dr7{"DR-0701"} = "HLA00719";		# DRB1*07:01:01:01
$group{"DR-0701"} = "DR7";
$base{"DR-0701"} = "DR7";

my %dr8;
$dr8{"DR-0801"} = "HLA00723";		# DRB1*08:01:01
$dr8{"DR-0803"} = "HLA00727";		# DRB1*08:03:02:01 266 bp
$dr8{"DR-0808"} = "HLA00734";		# DRB1*08:08 266 bp
$dr8{"DR-0818"} = "HLA00744";		# DRB1*08:18, used to replace partial DRB1*08:05
$group{"DR-0801"} = "DR8"; $group{"DR-0803"} = "DR8"; $group{"DR-0808"} = "DR8"; $group{"DR-0818"} = "DR8";
$base{"DR-0801"} = "DR8"; $base{"DR-0803"} = "DR8"; $base{"DR-0808"} = "DR8"; $base{"DR-0818"} = "DR8";
push @subtype, ("DR-0818","DR-0803","DR-0808");

my %dr9;
$dr9{"DR-0901"} = "HLA00749";		# DRB1*09:01:02:01
#$dr9{"DR-0902"} = "HLA01513";		# DRB1*09:02:01
$group{"DR-0901"} = "DR9";# $group{"DR-0902"} = "DR9";
$base{"DR-0901"} = "DR9";# $base{"DR-0902"} = "DR9";
#push @subtype, ("DR-0902");

my %dr10;
$dr10{"DR-1001"} = "HLA00750";		# DRB1*10:01:01:01
$group{"DR-1001"} = "DR10";
$base{"DR-1001"} = "DR10";
my %dr11;
$dr11{"DR-1101"} = "HLA00751";		# DRB1*11:01:01:01
$dr11{"DR-1102"} = "HLA00754";		# DRB1*11:02:01:01
$dr11{"DR-1103"} = "HLA00755";		# DRB1*11:03:01 266 bp
$dr11{"DR-1105"} = "HLA00758";		# DRB1*11:05
$dr11{"DR-1107"} = "HLA00760";		# DRB1*11:07:01
$dr11{"DR-1108"} = "HLA00761";		# DRB1*11:08:01 266 bp
$dr11{"DR-1117"} = "HLA00771";		# DRB1*11:17
$group{"DR-1101"} = "DR11"; $group{"DR-1102"} = "DR11"; $group{"DR-1103"} = "DR11"; $group{"DR-1105"} = "DR11"; $group{"DR-1107"} = "DR11"; $group{"DR-1108"} = "DR11"; $group{"DR-1117"} = "DR11";
$base{"DR-1101"} = "DR11"; $base{"DR-1102"} = "DR11"; $base{"DR-1103"} = "DR11"; $base{"DR-1105"} = "DR11"; $base{"DR-1107"} = "DR11"; $base{"DR-1108"} = "DR11"; $base{"DR-1117"} = "DR11";
push @subtype, ("DR-1102","DR-1103","DR-1105","DR-1107","DR-1108","DR-1117");

my %dr12;
$dr12{"DR-1201"} = "HLA00789";		# DRB1*12:01:01:01
$dr12{"DR-1202"} = "HLA00790";		# DRB1*12:02:01:01 266 bp
$group{"DR-1201"} = "DR12";$group{"DR-1202"} = "DR12";
$base{"DR-1201"} = "DR12";$base{"DR-1202"} = "DR12";

my %dr13;
$dr13{"DR-1301"} = "HLA00797";		# DRB1*13:01:01:01
$dr13{"DR-1303"} = "HLA00799";		# DRB1*13:03:01
$dr13{"DR-1305"} = "HLA00802";		# DRB1*13:05:01
$dr13{"DR-1312"} = "HLA00810";		# DRB1*13:12:01 266 bp
$dr13{"DR-1317"} = "HLA00815";		# DRB1*13:17
$dr13{"DR-1320"} = "HLA00818";		# DRB1*13:20:01:01 266 bp
$dr13{"DR-1339"} = "HLA01166";		# DRB1*13:39
$dr13{"DR-1343"} = "HLA01237";		# DRB1*13:43
$dr13{"DR-1349"} = "HLA30961";		# DRB1*13:49:02 183 bp
$group{"DR-1301"} = "DR13"; $group{"DR-1303"} = "DR13"; $group{"DR-1305"} = "DR13"; $group{"DR-1312"} = "DR13"; $group{"DR-1317"} = "DR13";
$group{"DR-1320"} = "DR13"; $group{"DR-1339"} = "DR13"; $group{"DR-1343"} = "DR13"; $group{"DR-1349"} = "DR13";
$base{"DR-1301"} = "DR13"; $base{"DR-1303"} = "DR13"; $base{"DR-1305"} = "DR13"; $base{"DR-1312"} = "DR13"; $base{"DR-1317"} = "DR13";
$base{"DR-1320"} = "DR13"; $base{"DR-1339"} = "DR13"; $base{"DR-1343"} = "DR13"; $base{"DR-1349"} = "DR13";

push @subtype, ("DR-1202","DR-1303","DR-1305","DR-1312","DR-1317","DR-1320","DR-1339","DR-1343","DR-1349");

my %dr14;
$dr14{"DR-1401"} = "HLA02371";		# DRB1*14:54:01:01
$dr14{"DR-1402"} = "HLA00834";		# DRB1*14:02:01:01
$dr14{"DR-1403"} = "HLA00835";		# DRB1*14:03:01
$dr14{"DR-1404"} = "HLA00836";		# DRB1*14:04:01:01
$dr14{"DR-1405"} = "HLA00837";		# DRB1*14:05:01:01
$dr14{"DR-1411"} = "HLA00843";		# DRB1*14:11
$dr14{"DR-1414"} = "HLA00846";		# DRB1*14:14
$dr14{"DR-1417"} = "HLA00849";		# DRB1*14:17 266 bp
$dr14{"DR-1419"} = "HLA00851";		# DRB1*14:19
$dr14{"DR-1422"} = "HLA00854";		# DRB1*14:22
$dr14{"DR-1424"} = "HLA00856";		# DRB1*14:24
$dr14{"DR-1448"} = "HLA01734";		# DRB1*14:48
$group{"DR-1401"} = "DR14"; $group{"DR-1402"} = "DR14"; $group{"DR-1403"} = "DR14"; $group{"DR-1404"} = "DR14"; $group{"DR-1405"} = "DR14"; $group{"DR-1411"} = "DR14";
$group{"DR-1414"} = "DR14"; $group{"DR-1417"} = "DR14"; $group{"DR-1419"} = "DR14"; $group{"DR-1422"} = "DR14"; $group{"DR-1424"} = "DR14"; $group{"DR-1448"} = "DR14";
$base{"DR-1401"} = "DR14"; $base{"DR-1402"} = "DR14"; $base{"DR-1403"} = "DR1403"; $base{"DR-1404"} = "DR1404"; $base{"DR-1405"} = "DR14"; $base{"DR-1411"} = "DR14";
$base{"DR-1414"} = "DR14"; $base{"DR-1417"} = "DR14"; $base{"DR-1419"} = "DR14"; $base{"DR-1422"} = "DR14"; $base{"DR-1424"} = "DR14"; $base{"DR-1448"} = "DR14";

push @subtype, ("DR-1402","DR-1403","DR-1404","DR-1405","DR-1411","DR-1414","DR-1417","DR-1419","DR-1422","DR-1424","DR-1448");	# modify here if serotype modified

my %dr15;
$dr15{"DR-1501"} = "HLA00865";		# DRB1*15:01:01:01
$dr15{"DR-1504"} = "HLA00871";		# DRB1*15:04 266 bp
$dr15{"DR-1505"} = "HLA00872";		# DRB1*15:05
$group{"DR-1501"} = "DR15";$group{"DR-1504"} = "DR15";$group{"DR-1505"} = "DR15";
$base{"DR-1501"} = "DR15";$base{"DR-1504"} = "DR15";$base{"DR-1505"} = "DR15";
my %dr16;
$dr16{"DR-1601"} = "HLA00876";		# DRB1*16:01:01
$dr16{"DR-1602"} = "HLA00878";		# DRB1*16:02:01:01 266 bp
$dr16{"DR-1605"} = "HLA00882";		# DRB1*16:05:01 266 bp
$group{"DR-1601"} = "DR16";$group{"DR-1602"} = "DR16";$group{"DR-1605"} = "DR16";
$base{"DR-1601"} = "DR16";$base{"DR-1602"} = "DR16";$base{"DR-1605"} = "DR16";
push @subtype, ("DR-1504","DR-1602","DR-1605","DR-1505");


sub DRB1 {
	my $gene = "DRB1";
	return $gene;
}

sub DRB1_LEADER {
	my $leader = 28;		# DRB1 specific
	return $leader;
}

sub GROUP {
	my $group_ref = \%group;
	return $group_ref;
}

sub BASE {
	my $base_ref = \%base;
	return $base_ref;
}

sub BASETYPE {
	my @basetype = ("DR1","DR103","DR15","DR16","DR17","DR18","DR4","DR7","DR8","DR9","DR10",
		"DR11","DR12","DR13","DR14","DR1403","DR1404");
	my $basetype_ref = \@basetype;
}

sub BROAD {
	my %broad;
	foreach my $base ( keys %base ) {
		if (( $base{ $base } eq "DR1" ) || ( $base{ $base } eq "DR103" )) {
			$broad{ $base } = "DR1";
		}
		elsif (( $base{ $base } eq "DR15" ) || ( $base{ $base } eq "DR16" )) {
			$broad{ $base } = "DR2";
		}
		elsif (( $base{ $base } eq "DR17" ) || ( $base{ $base } eq "DR18" )) {
			$broad{ $base } = "DR3";
		}
		elsif (( $base{ $base } eq "DR11" ) || ( $base{ $base } eq "DR12" )) {
			$broad{ $base } = "DR5";
		}
		elsif (( $base{ $base } eq "DR13" ) || ( $base{ $base } eq "DR14" ) || ( $base{ $base } eq "DR1403" ) || ( $base{ $base } eq "DR1404" )) {
			$broad{ $base } = "DR6";
		}
		else {
			$broad{ $base } = $base{ $base };
		}
	}
	my $broad_ref = \%broad;
	return $broad_ref;
}

sub RESIDUES {
	my ( $serotype ) = @_;
	my @combined = ();
	push @combined, @dr1;
	push @combined, @dr3;
	push @combined, @dr4;
	push @combined, @dr7;
	push @combined, @dr8;
	push @combined, @dr9;
	push @combined, @dr10;
	push @combined, @dr11;
	push @combined, @dr12;
	push @combined, @dr13;
	push @combined, @dr14;
	push @combined, @dr15;
	push @combined, @dr16;
	my %seen;
	my @unique;
	foreach my $value ( sort { $a <=> $b } @combined ) {
		unless ( exists $seen{ $value } ) {
			push @unique, $value;
			$seen{ $value } = 0;
		}
	}
	
	my @residues = ();
	my $residues_ref = \@residues;
	if ( $serotype eq "DR1" ) {
		@residues = @dr1; 
	}
	elsif ( $serotype eq "DR3" ) {
		@residues = @dr3;
	}
	elsif ( $serotype eq "DR4" ) {
		@residues = @dr4;
	}
	elsif ( $serotype eq "DR7" ) {
		@residues = @dr7;
	}
	elsif ( $serotype eq "DR8" ) {
		@residues = @dr8; 
	}
	elsif ( $serotype eq "DR9" ) {
		@residues = @dr9;
	}
	elsif ( $serotype eq "DR10" ) {
		@residues = @dr10;
	}
	elsif ( $serotype eq "DR11" ) {
		@residues = @dr11;
	}
	elsif ( $serotype eq "DR12" ) {
		@residues = @dr12;
	}
	elsif ( $serotype eq "DR13" ) {
		@residues = @dr13;
	}
	elsif ( $serotype eq "DR14" ) {
		@residues = @dr14;
	}
	elsif ( $serotype eq "DR15" ) {
		@residues = @dr15;
	}
	elsif ( $serotype eq "DR16" ) {
		@residues = @dr16;
	}
	else {
		@residues = @unique;
	}
	return $residues_ref;
}

sub REF {
	my ( $serotype ) = @_;
	my %ref;
	my $ref_ref = \%ref;

	if ( $serotype eq "DR1" ) {
		%ref = %dr1; 
	}
	elsif ( $serotype eq "DR3" ) {
		%ref = %dr3;
	}
	elsif ( $serotype eq "DR4" ) {
		%ref = %dr4;
	}
	elsif ( $serotype eq "DR7" ) {
		%ref = %dr7;
	}
	elsif ( $serotype eq "DR8" ) {
		%ref = %dr8;
	}
	elsif ( $serotype eq "DR9" ) {
		%ref = %dr9;
	}
	elsif ( $serotype eq "DR10" ) {
		%ref = %dr10;
	}
	elsif ( $serotype eq "DR11" ) {
		%ref = %dr11;
	}
	elsif ( $serotype eq "DR12" ) {
		%ref = %dr12;
	}
	elsif ( $serotype eq "DR13" ) {
		%ref = %dr13;
	}
	elsif ( $serotype eq "DR14" ) {
		%ref = %dr14;
	}
	elsif ( $serotype eq "DR15" ) {
		%ref = %dr15;
	}
	elsif ( $serotype eq "DR16" ) {
		%ref = %dr16;
	}
	else {
		%ref = (%dr1,%dr3,%dr4,%dr7,%dr8,%dr9,%dr10,%dr11,%dr12,%dr13,%dr14,%dr15,%dr16);
	}
	
	return $ref_ref;
}

sub SERO {
	my @sero;
	my %ref = (%dr1,%dr3,%dr4,%dr7,%dr8,%dr9,%dr10,%dr11,%dr12,%dr13,%dr14,%dr15,%dr16);
	my @tmp = sort keys %ref;	# all keys
	for ( my $index = 0; $index < scalar @tmp; $index++ ) {
		$sero[0][$index] = $tmp[$index];
	}
	for ( my $index = 0; $index < scalar @subtype; $index++ ) {
		$sero[1][$index] = $subtype[$index];
	}
	my $sero_ref = \@sero;
	return $sero_ref;
}

sub KEY {
	my %tmp = (%dr1,%dr3,%dr4,%dr7,%dr8,%dr9,%dr10,%dr11,%dr12,%dr13,%dr14,%dr15,%dr16);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if (( $key eq "DR-0101" ) || ( $key eq "DR-0103" )) {
			$ref{$key} = "DRB1\\*01";
		}
		elsif ( $key =~ /DR-03/ ) {
			$ref{$key} = "DRB1\\*03";
		}
		elsif ( $key =~ /DR-04/ ) {
			$ref{$key} = "DRB1\\*04";
		}
		elsif ( $key =~ /DR-07/ ) {
			$ref{$key} = "DRB1\\*07";
		}
		elsif ( $key =~ /DR-08/ ) {
			$ref{$key} = "DRB1\\*08";
		}
		elsif ( $key =~ /DR-09/ ) {
			$ref{$key} = "DRB1\\*09";
		}
		elsif ( $key =~ /DR-10/ ) {
			$ref{$key} = "DRB1\\*10";
		}
		elsif ( $key =~ /DR-11/ ) {
			$ref{$key} = "DRB1\\*11";
		}
		elsif ( $key =~ /DR-12/ ) {	# DR12, DR1209, DR1221
			$ref{$key} = "DRB1\\*12";
		}
		elsif ( $key =~ /DR-13/ ) {
			$ref{$key} = "DRB1\\*13";
		}
		elsif ( $key =~ /DR-14/ ) {
			$ref{$key} = "DRB1\\*14";
		}
		elsif ( $key =~ /DR-15/ ) {
			$ref{$key} = "DRB1\\*15";
		}
		else {
			$ref{$key} = "DRB1\\*16";
		}
	}
	return $key_ref;

}

sub DR5231 {	# this is a place holder for consistency
	my $dr5231_ref = \%dr5231;
	return $dr5231_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "X" x 34;

#	$partial{ "DR-0902" } = $seq;	# still partial
	$partial{ "DR-1339" } = $seq;	# still partial
	$partial{ "DR-1343" } = $seq;	# still partial
	$partial{ "DR-1349" } = $seq;	# still partial
	$partial{ "DR-1422" } = $seq;	# still partial
	$partial{ "DR-1107" } = $seq;	# still partial
	$partial{ "DR-1448" } = $seq;	# still partial
	$partial{ "DR-1505" } = $seq;	# still partial
	$partial{ "DRB1*03:10" } = "X" x 29;	# still partial
	$partial{ "DRB1*04:20" } = "X" x 42;	# unusual partial sequence
	$partial{ "DRB1*04:22" } = "X" x 37;	# unusual partial sequence
	$partial{ "DRB1*11:26" } = "X" x 36;	# unusual partial sequence
	$partial{ "DRB1*13:30" } = "X" x 35;	# unusual partial sequence
	$partial{ "DRB1*13:55" } = "X" x 35;	# unusual partial sequence
	$partial{ "DRB1*14:31" } = "X" x 35;	# unusual partial sequence
	$partial{ "general" } = $seq;
		
	return $partial_ref;
}

sub WHO {
	my %who;
	my $whotype_ref = \%who;
	$who{"DR-0101"} = "DR1"; $who{"DR-0103"} = "DR103"; $who{"DR-1501"} = "DR15"; $who{"DR-1601"} = "DR16"; $who{"DR-0301"} = "DR17";
	$who{"DR-0302"} = "DR18"; $who{"DR-0401"} = "DR4"; $who{"DR-0701"} = "DR7"; $who{"DR-0801"} = "DR8"; $who{"DR-0901"} = "DR9";
	$who{"DR-1001"} = "DR10"; $who{"DR-1101"} = "DR11"; $who{"DR-1201"} = "DR12"; $who{"DR-1301"} = "DR13"; $who{"DR-1401"} = "DR14";
	$who{"DR-1403"} = "DR1403"; $who{"DR-1404"} = "DR1404";

	return $whotype_ref;
}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "DR-1417" } = 0;
	$known_cross{ "DR-1349" } = 0;
	return $known_cross_ref;
}

1;
