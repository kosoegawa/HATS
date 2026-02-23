#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: HLAC_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on February 19 2026

package HLAC_INFO;
use strict;

# residue 80 for C1 & C2
# residue 76 to distinguish from HLA-B
# removed residue 173 completely on January 24 2025
my @cw1 = (66,73,76,77,80,82,83,99,138,163,177);	# added 138, 99 is requied to separate from Cw14 and Cw16
my @cw2 = (45,66,73,76,77,80,82,83,143,163,177);	#143 to distinguish from Cw17
my @cw3 = (45,49,66,76,77,80,82,83,91,103,147,163,177);	# 91 to separate Cw10 from Cw9
my @cw4 = (24,49,66,73,76,77,80,90,91,163);	# 49 to separate Cw4 from Cw403, removed 99, added 90
my @cw5 = (49,73,76,77,80,138,163,177);		# residue 138 is required
my @cw6 = (24,49,66,73,76,77,80,82,83,99,143,163,175);	# 143: to distinguish from Cw17; 99 is required to distinguish from Cw18
my @cw7 = (66,73,76,77,80,90,147,163,177);	#added 90
my @cw8 = (73,76,77,80,138,163,175,177);
my @cw12 = (73,76,77,80,82,83,90,147,158,163,177);		# 158 was included to remove B*67:02
my @cw1212 = (45,73,76,77,80,82,83,90,147,158,163);
my @cw1204 = (24,49,66,73,76,77,80,82,83,90,147,158,163);	# removed 99, added 90
my @cw14 = (73,76,77,80,99,163,177);	# 99 is required to separate from Cw1 and Cw16
my @cw15 = (66,73,76,77,80,82,83,163,177);
my @cw16 = (45,66,73,76,77,80,99,138,163,177);	# 99 is required to separate from Cw1 and Cw14
my @cw17 = (76,77,80,143,147,163);
my @cw18 = (24,49,73,76,77,80,99);
my @extra = (109);	# FULL only

my %c1c2;
my %group;	# conventional group to chose specific residues
my %base;	# recognized serotype
my @subtype;

my %cw1;
$cw1{"Cw0102"} = "HLA00401";	# C*01:02:01:01
#$cw1{"Cw0114"} = "HLA02790";	# C*01:14 Dropped on 1/24/25
$group{"Cw0102"} = "Cw1";# $group{"Cw0114"} = "Cw1";
$base{"Cw0102"} = "Cw1";# $base{"Cw0114"} = "Cw1";
$c1c2{"Cw0102"} = "Negative";# $c1c2{"Cw0114"} = "Negative";
#push @subtype, ("Cw0114");		#,"Cw195"

my %cw2;
$cw2{"Cw0202"} = "HLA00404";	# C*02:02:01
$group{"Cw0202"} = "Cw2";
$base{"Cw0202"} = "Cw2";
$c1c2{"Cw0202"} = "Negative";

#$cw2{"Cw0227"} = "HLA03705";	# C*02:27:01 Dropped on 1/24/25
#$group{"Cw0227"} = "Cw2";
#$base{"Cw0227"} = "Cw2";
#$c1c2{"Cw0227"} = "Negative";
#push @subtype, ("Cw0227");	#,"Cw2131"

my %cw3;
$cw3{"Cw0303"} = "HLA00411";	# C*03:03:01:01
$cw3{"Cw0304"} = "HLA00413";	# C*03:04:01:01
$cw3{"Cw0307"} = "HLA00417";	# C*03:07:01:01
$cw3{"Cw0308"} = "HLA00418";	# C*03:08 366 bp
$group{"Cw0303"} = "Cw3"; $group{"Cw0304"} = "Cw3"; $group{"Cw0307"} = "Cw3"; $group{"Cw0308"} = "Cw3";
$base{"Cw0303"} = "Cw9"; $base{"Cw0304"} = "Cw10"; $base{"Cw0307"} = "Cw10"; $base{"Cw0308"} = "Cw10";

$c1c2{"Cw0303"} = "Negative"; $c1c2{"Cw0304"} = "Negative"; $c1c2{"Cw0307"} = "Negative"; $c1c2{"Cw0308"} = "Negative";
push @subtype, ("Cw0307","Cw0304","Cw0308");	#"Cw370","Cw3127","Cw3450",

my %cw4;
$cw4{"Cw0401"} = "HLA00420";	# C*04:01:01:01
$cw4{"Cw0403"} = "HLA00423";	# C*04:03:01:01
$cw4{"Cw0408"} = "HLA01311";	# C*04:08
$cw4{"Cw0410"} = "HLA01645";	# C*04:10
$cw4{"Cw0427"} = "HLA02831";	# C*04:27 366 bp
$group{"Cw0401"} = "Cw4"; $group{"Cw0403"} = "Cw4"; $group{"Cw0408"} = "Cw4"; $group{"Cw0410"} = "Cw4"; $group{"Cw0427"} = "Cw4";
$base{"Cw0401"} = "Cw4"; $base{"Cw0403"} = "Cw4"; $base{"Cw0408"} = "Cw4"; $base{"Cw0410"} = "Cw4"; $base{"Cw0427"} = "Cw4";
$c1c2{"Cw0401"} = "Negative"; $c1c2{"Cw0403"} = "Negative"; $c1c2{"Cw0408"} = "Negative"; $c1c2{"Cw0410"} = "Negative"; $c1c2{"Cw0427"} = "Negative";
push @subtype, ("Cw0403","Cw0408","Cw0410","Cw0427");

my %cw5;
$cw5{"Cw0501"} = "HLA00427";	#C*05:01:01:01
$cw5{"Cw0509"} = "HLA01856";	#C*05:09:01
$cw5{"Cw0810"} = "HLA01835";	# C*08:10
$group{"Cw0501"} = "Cw5"; $group{"Cw0509"} = "Cw5"; $group{"Cw0810"} = "Cw5";
$base{"Cw0501"} = "Cw5"; $base{"Cw0509"} = "Cw5"; $base{"Cw0810"} = "Cw5";
$c1c2{"Cw0501"} = "Negative"; $c1c2{"Cw0509"} = "Negative"; $c1c2{"Cw0810"} = "Negative";
push @subtype, ("Cw0509","Cw0810");

my %cw6;
$cw6{"Cw0602"} = "HLA00430";	# C*06:02:01:01
$cw6{"Cw0608"} = "HLA01628";	# C*06:08
$cw6{"Cw0627"} = "HLA04299";	# C*06:27 366 bp
$group{"Cw0602"} = "Cw6"; $group{"Cw0608"} = "Cw6"; $group{"Cw0627"} = "Cw6";	
$base{"Cw0602"} = "Cw6"; $base{"Cw0608"} = "Cw6"; $base{"Cw0627"} = "Cw6";
$c1c2{"Cw0602"} = "Negative"; $c1c2{"Cw0608"} = "Negative"; $c1c2{"Cw0627"} = "Negative";
push @subtype, ("Cw0608","Cw0627");	

my %cw7;
$cw7{"Cw0701"} = "HLA00433";	# C*07:01:01:01
$cw7{"Cw0702"} = "HLA00434";	# C*07:02:01:01
$cw7{"Cw0704"} = "HLA00436";	# C*07:04:01:01
$cw7{"Cw0707"} = "HLA00439";	# C*07:07 
$cw7{"Cw0717"} = "HLA01662";	# C*07:17:01:01
$group{"Cw0701"} = "Cw7"; $group{"Cw0702"} = "Cw7"; $group{"Cw0704"} = "Cw7"; $group{"Cw0717"} = "Cw7"; $group{"Cw0707"} = "Cw7";
$base{"Cw0701"} = "Cw7"; $base{"Cw0702"} = "Cw7"; $base{"Cw0704"} = "Cw7"; $base{"Cw0717"} = "Cw7"; $base{"Cw0707"} = "Cw7";
$c1c2{"Cw0701"} = "Negative";  $c1c2{"Cw0702"} = "Negative"; $c1c2{"Cw0704"} = "Negative"; $c1c2{"Cw0717"} = "Negative"; $c1c2{"Cw0707"} = "Negative";
push @subtype, ("Cw0702","Cw0704","Cw0717","Cw0707");	#

my %cw8;
$cw8{"Cw0801"} = "HLA00445";	# C*08:01:01:01
$cw8{"Cw0802"} = "HLA00446";	# C*08:02:01:01
$cw8{"Cw0803"} = "HLA00447";	# C*08:03:01
$cw8{"Cw0806"} = "HLA00450";	# C*08:06 366 bp
$group{"Cw0801"} = "Cw8"; $group{"Cw0802"} = "Cw8"; $group{"Cw0803"} = "Cw8"; $group{"Cw0806"} = "Cw8";
$base{"Cw0801"} = "Cw8"; $base{"Cw0802"} = "Cw8"; $base{"Cw0803"} = "Cw8"; $base{"Cw0806"} = "Cw8";
$c1c2{"Cw0801"} = "Negative"; $c1c2{"Cw0802"} = "Negative"; $c1c2{"Cw0803"} = "Negative"; $c1c2{"Cw0806"} = "Negative";
push @subtype, ("Cw0802","Cw0803","Cw0806");

my %cw12;
$cw12{"Cw1202"} = "HLA00453";	# C*12:02:01
$group{"Cw1202"} = "Cw12";
$base{"Cw1202"} = "Cw12"; 
$c1c2{"Cw1202"} = "Negative"; 

my %cw1212;
$cw1212{"Cw1212"} = "HLA01878";	# C*12:12, C1 Bw6
$group{"Cw1212"} = "Cw1212";
$base{"Cw1212"} = "Cw12";
$c1c2{"Cw1212"} = "Bw6";

my %cw1204;
$cw1204{"Cw1204"} = "HLA00456";	# C*12:04:01
$group{"Cw1204"} = "Cw1204";
$base{"Cw1204"} = "Cw12";
$c1c2{"Cw1204"} = "Negative";
push @subtype, ("Cw1212", "Cw1204");

my %cw14;
$cw14{"Cw1402"} = "HLA00462";	# C*14:02:01:01
$group{"Cw1402"} = "Cw14";
$base{"Cw1402"} = "Cw14";
$c1c2{"Cw1402"} = "Negative";

my %cw15;
$cw15{"Cw1502"} = "HLA00467";	# C*15:02:01:01
$cw15{"Cw1507"} = "HLA00473";	# C*15:07:01:01
$group{"Cw1502"} = "Cw15"; $group{"Cw1507"} = "Cw15";
$base{"Cw1502"} = "Cw15"; $base{"Cw1507"} = "Cw15";
$c1c2{"Cw1502"} = "Negative"; $c1c2{"Cw1507"} = "Negative";
push @subtype, ("Cw1507");	#"Cw1545",

my %cw16;
$cw16{"Cw1601"} = "HLA00475";	# C*16:01:01:01
$cw16{"Cw1602"} = "HLA00476";	# C*16:02:01:01
$group{"Cw1601"} = "Cw16"; $group{"Cw1602"} = "Cw16";
$base{"Cw1601"} = "Cw16"; $base{"Cw1602"} = "Cw16";
$c1c2{"Cw1601"} = "Negative"; $c1c2{"Cw1602"} = "Negative";
push @subtype, ("Cw1602");

my %cw17;
$cw17{"Cw1701"} = "HLA04311";	# C*17:01:01:02
$group{"Cw1701"} = "Cw17";
$base{"Cw1701"} = "Cw17";
$c1c2{"Cw1701"} = "Negative";

my %cw18;
$cw18{"Cw1801"} = "HLA00483";	# C*18:01:01:01
$group{"Cw1801"} = "Cw18";
$base{"Cw1801"} = "Cw18";
$c1c2{"Cw1801"} = "Negative";


sub HLAC {
	my $gene = "C";
	return $gene;
}

sub HLAC_LEADER {
	my $leader = 23;		# C specific
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
	my @basetype = ("Cw1","Cw2","Cw9","Cw10","Cw4","Cw5","Cw6","Cw7","Cw8","Cw12","Cw14","Cw15","Cw16","Cw17","Cw18");
	my $basetype_ref = \@basetype;
}


sub PARENT {
	my %parent;
	my $parent_ref = \%parent;
	foreach my $key ( keys %base ) {	# $key = A0101
		$parent{ $key } = $base{ $key };	# $base( $key } = "A1";
	}
	return $parent_ref;
}


sub BROAD {
	my %broad;
	foreach my $base ( keys %base ) {
		unless (( $base{ $base } eq "Cw9" ) || ( $base{ $base } eq "Cw10" )) {
#			if (( $base{ $base } eq "Cw12" ) || ( $base{ $base } eq "Cw14" ) || ( $base{ $base } eq "Cw15" ) ||
#			( $base{ $base } eq "Cw16" ) || ( $base{ $base } eq "Cw17" ) || ( $base{ $base } eq "Cw18" )) {
#				$broad{ $base } = "";
#			}
#			else {
#				$broad{ $base } = $base{ $base };
#			}
			$broad{ $base } = $base{ $base };
		}
		else {
			$broad{ $base } = "Cw3";
		}
	}
	my $broad_ref = \%broad;
	return $broad_ref;
}

sub RESIDUES {
	my ( $serotype ) = @_;
	my @combined = ();
	push @combined, @cw1; 
	push @combined, @cw2; 
	push @combined, @cw3; 
	push @combined, @cw4; 
	push @combined, @cw5;
	push @combined, @cw6; 
	push @combined, @cw7; 
	push @combined, @cw8; 
	push @combined, @cw12; 
	push @combined, @cw1204; 
	push @combined, @cw1212;
	push @combined, @cw14; 
	push @combined, @cw15; 
	push @combined, @cw16; 
	push @combined, @cw17; 
	push @combined, @cw18; 
	push @combined, @extra; 

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
	if ( $serotype eq "Cw1" ) {
		@residues = @cw1; 
	}
	elsif ( $serotype eq "Cw2" ) {
		@residues = @cw2;
	}
	elsif ( $serotype eq "Cw3" ) {
		@residues = @cw3;
	}
	elsif ( $serotype eq "Cw4" ) {
		@residues = @cw4;
	}
	elsif ( $serotype eq "Cw5" ) {
		@residues = @cw5;
	}
	elsif ( $serotype eq "Cw6" ) {
		@residues = @cw6;
	}
	elsif ( $serotype eq "Cw7" ) {
		@residues = @cw7;
	}
	elsif ( $serotype eq "Cw8" ) {
		@residues = @cw8;
	}
	elsif ( $serotype eq "Cw12" ) {
		@residues = @cw12;
	}
	elsif ( $serotype eq "Cw1204" ) {
		@residues = @cw1204;
	}
	elsif ( $serotype eq "Cw1212" ) {
		@residues = @cw1212;
	}
	elsif ( $serotype eq "Cw14" ) {
		@residues = @cw14;
	}
	elsif ( $serotype eq "Cw15" ) {
		@residues = @cw15;
	}
	elsif ( $serotype eq "Cw16" ) {
		@residues = @cw16;
	}
	elsif ( $serotype eq "Cw17" ) {
		@residues = @cw17;
	}
	elsif ( $serotype eq "Cw17" ) {
		@residues = @cw17;
	}
	elsif ( $serotype eq "Cw18" ) {
		@residues = @cw18;
	}
	else {
		@residues = @unique;
	}
	return $residues_ref;
}

sub RESIDUES_ABC {
	my @a = (43,44,45,56,62,63,65,66,67,73,74,76,82,83,107,127,144,145,149,151,161,163,166,167,171);
	my @b = (11,45,46,62,63,67,69,70,71,76,82,83,103,127,143,145,147,158,163,167,171,177,178,180);
	my @c = (24,45,49,66,73,76,77,80,82,83,91,99,138,143,147,158,163,173,175,177);

	my ( $serotype ) = @_;
	my @combined = ();
	push @combined, @a; 
	push @combined, @b; 
	push @combined, @c; 

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
	if ( $serotype eq "ALL" ) {
		@residues = @unique;
	}
	return $residues_ref;
}

sub REF {
	my ( $serotype ) = @_;
	my %ref;
	my $ref_ref = \%ref;

	if ( $serotype eq "Cw1" ) {
		%ref = %cw1; 
	}
	elsif ( $serotype eq "Cw2" ) {
		%ref = %cw2;
	}
	elsif ( $serotype eq "Cw3" ) {
		%ref = %cw3;
	}
	elsif ( $serotype eq "Cw4" ) {
		%ref = %cw4;
	}
	elsif ( $serotype eq "Cw5" ) {
		%ref = %cw5;
	}
	elsif ( $serotype eq "Cw6" ) {
		%ref = %cw6;
	}
	elsif ( $serotype eq "Cw7" ) {
		%ref = %cw7;
	}
	elsif ( $serotype eq "Cw8" ) {
		%ref = %cw8;
	}
	elsif ( $serotype eq "Cw12" ) {
		%ref = %cw12;
	}
	elsif ( $serotype eq "Cw1204" ) {
		%ref = %cw1204;
	}
	elsif ( $serotype eq "Cw1212" ) {
		%ref = %cw1212;
	}
	elsif ( $serotype eq "Cw14" ) {
		%ref = %cw14;
	}
	elsif ( $serotype eq "Cw15" ) {
		%ref = %cw15;
	}
	elsif ( $serotype eq "Cw16" ) {
		%ref = %cw16;
	}
	elsif ( $serotype eq "Cw17" ) {
		%ref = %cw17;
	}
	elsif ( $serotype eq "Cw18" ) {
		%ref = %cw18;
	}
	else {
		%ref = (%cw1,%cw2,%cw3,%cw4,%cw5,%cw6,%cw7,%cw8,%cw12,%cw1212,%cw1204,%cw14,%cw15,%cw16,%cw17,%cw18);
	}		
	
	return $ref_ref;
}


sub SERO {
	my @sero;
	my %ref = (%cw1,%cw2,%cw3,%cw4,%cw5,%cw6,%cw7,%cw8,%cw12,%cw1212,%cw1204,%cw14,%cw15,%cw16,%cw17,%cw18);
	my @tmp = sort keys %ref;
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
	my %tmp = (%cw1,%cw2,%cw3,%cw4,%cw5,%cw6,%cw7,%cw8,%cw12,%cw1212,%cw1204,%cw14,%cw15,%cw16,%cw17,%cw18);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
#		if (( $key eq "C-0102" ) || ( $key eq "C-0114" )) {
		if ( $key eq "Cw0102" ) {
			$ref{$key} = "C\\*01";
		}
		elsif ( $key =~ /Cw02/ ) {
			$ref{$key} = "C\\*02";
		}
		elsif (( $key =~ /Cw03/ )) {
			$ref{$key} = "C\\*03";
		}
		elsif ( $key =~/Cw04/ ) {
			$ref{$key} = "C\\*04";
		}
		elsif ( $key =~ /Cw05/ ) {
			$ref{$key} = "C\\*05";
		}
		elsif ( $key =~ /Cw06/ ) {
			$ref{$key} = "C\\*06";
		}
		elsif ( $key =~ /Cw07/ ) {
			$ref{$key} = "C\\*07";
		}
		elsif ( $key =~ /Cw08/ ) {
			$ref{$key} = "C\\*08";
		}
		elsif ( $key =~ /Cw12/ ) {
			$ref{$key} = "C\\*12";
		}
		elsif ( $key =~ /Cw14/ ) {
			$ref{$key} = "C\\*14";
		}
		elsif ( $key =~ /Cw15/ ) {
			$ref{$key} = "C\\*15";
		}
		elsif ( $key =~ /Cw16/ ) {
			$ref{$key} = "C\\*16";
		}
		elsif ( $key =~ /Cw17/ ) {
			$ref{$key} = "C\\*17";
		}
		elsif ( $key =~ /Cw18/ ) {
			$ref{$key} = "C\\*18";
		}
	}
	return $key_ref;

}

sub BW {
	my $bw_ref = \%c1c2;
	return $bw_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "N" x 25;	#change the number of missing nucleotide
	$partial{ "general" } = $seq;
	# C*03:46 contains 1 AA deletion, so the AA alignment is not correct
		
	return $partial_ref;
}

#sub WHO {
#	my %who;
#	my $whotype_ref = \%who;
#	$who{"C-0102"} = "Cw1"; $who{"C-0202"} = "Cw2"; $who{"C-0303"} = "Cw9"; $who{"C-0304"} = "Cw10"; $who{"C-0401"} = "Cw4";
#	$who{"C-0501"} = "Cw5"; $who{"C-0602"} = "Cw6"; $who{"C-0701"} = "Cw7"; $who{"C-0801"} = "Cw8";

#	return $whotype_ref;
#}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "NOTHING" } = 0;
	return $known_cross_ref;
}

1;
