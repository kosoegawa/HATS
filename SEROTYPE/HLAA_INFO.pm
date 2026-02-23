#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# © 2022 Stanford Blood Center L.L.C.
# SPDX-License-Identifier: BSD-3-Clause

# module: HLAA_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on January 22 2026

package HLAA_INFO;
use strict;


# It appears that 44 snd 67 are important for A1 & A36
# residue 163R was added 163RG is well established eplet
my @a1 = (44,67,76,163,166,167);

# A2 62G is mandatory, residues 43 & 73 were added to separate A202 and A211, 144 were added to separate A265
my @a2 = (43,62, 66, 73, 107,138, 144,145, 149,163);	#138 was added for A*02:18 
my @a219 = (62, 65, 82, 83, 144, 145,161,163, 166, 167);	#removed 151
my @a3 = (62, 76,144,161, 163);	# 161D is mandatory
my @a305 = (62,76, 144,145, 161, 163); 	# added 144 and 145 to exclude the other alleles
my @a11 = (62, 76, 144, 161, 163,167);
my @a9 = (62, 65, 82, 83, 107, 144, 151, 163, 166, 167);	#Bw4, 65G for A9, both 144 and 151 are important
my @a25 = (76, 82,83, 144, 149,163,166,167);	#Bw4
my @a26a43 = (62,63,74,76,144, 149,163,166);	 #
my @a29 = (44, 62, 63,102,144,151);		#residue 44 is included to exclude noise for SHORT
my @a30 = (56, 62, 73,76, 82,83,144);	# included 62 to separate A3007, removed 151 to include A*30:04	 
my @a31 = (56,62, 66, 73, 76, 82, 83, 144);	# 62 is included to eliminate A*33 alleles 
my @a32 = (62, 82, 83,109,144,161, 163,167);	#Bw4, 62 is used to separate from A2403, 109 was added
my @a33 = (62, 63, 73, 76, 144, 151, 171);	# residue 171 is unique for A33
my @a3313 = (62, 63, 73, 76, 144, 151,163, 171);	# residue 163 was added
my @a34 = (62, 66, 67, 144, 149, 163);
my @a36 = (44,67,76, 166, 167);	#67 is used to show cross-reactivity with A1, 76 was added to exclude noise
my @a66 = (62, 66, 74, 76, 149, 163);
my @a68 = (62,63, 76, 107, 144, 145);	# added 63
my @a6836 = (62, 76,82,83, 107, 144,145); # 63 was originally included to remove B*15:67 allele, but was able to replace with 144
my @a69 = (62, 76, 107, 144, 145);	#145 was added to exclude noise
my @a74 = (45,56,62, 63, 73,76,109,144, 161, 163);	#45 is used to eliminate HLA-C, residue 56 is used to distinguish from A30, 109 was added
my @a80 = (62,66,144,145, 151, 166, 167);
my @extra = (90,103,127);	# FULL only, added residue 90

my %bw;
my %a1;
my %group;	# conventional group to chose specific residues
my %base;	# recognized serotype
my @subtype;
$a1{"A0101"} = "HLA00001";		# A*01:01:01:01
$group{"A0101"} = "A1";
$base{"A0101"} = "A1";
my %a2;
$a2{"A0201"} = "HLA00005";
$a2{"A0202"} = "HLA00007";	# A*02:02:01:01
$a2{"A0203"} = "HLA00008";
$a2{"A0208"} = "HLA00013";	# A*02:08:01 365 bp
$a2{"A0211"} = "HLA00016";	# A*02:11:01:01
$a2{"A0216"} = "HLA00021";	# A*02:16 365 bp
$a2{"A0218"} = "HLA00024";	# A*02:18 365 bp, WD but in supplemental bead
$a2{"A0220"} = "HLA00026";	# A*02:20:01 365 bp
$a2{"A0285"} = "HLA02235";	# A*02:85
$a2{"A0256"} = "HLA01575";		# A*02:56:01
$group{"A0201"} = "A2"; $group{"A0202"} = "A2"; $group{"A0203"} = "A2"; $group{"A0208"} = "A2"; $group{"A0211"} = "A2";
$group{"A0216"} = "A2"; 
$group{"A0220"} = "A2"; $group{"A0285"} = "A2"; $group{"A0256"} = "A2";$group{"A0218"} = "A2";
$base{"A0201"} = "A2"; $base{"A0202"} = "A2"; $base{"A0203"} = "A2"; $base{"A0208"} = "A2"; $base{"A0211"} = "A2";
$base{"A0216"} = "A2"; 
$base{"A0220"} = "A2"; $base{"A0285"} = "A2"; $base{"A0256"} = "A2"; $base{"A0218"} = "A2";
#my %a210;
$a2{"A0210"} = "HLA00015";
$group{"A0210"} = "A2";
$base{"A0210"} = "A2";
push @subtype, ("A0202","A0203","A0208","A0210","A0211","A0216","A0220","A0285","A0256","A0218");#	
my%a219;
$a219{"A0219"} = "HLA00025";	# A*02:19, I group
$a219{"A0244"} = "HLA01222";	# A*02:44, I group
$a219{"A0323"} = "HLA02528";	# A*03:23:01
$group{"A0219"} = "A219"; $group{"A0244"} = "A219"; $group{"A0323"} = "A219";
$base{"A0219"} = "A2"; $base{"A0244"} = "A2"; $base{"A0323"} = "A3";
push @subtype, ("A0219","A0244","A0323");

my %a3;
$a3{"A0301"} = "HLA00037";		# A*03:01:01:01
$group{"A0301"} = "A3";
$base{"A0301"} = "A3";
my %a305;
$a305{"A0305"} = "HLA01107";		# A*03:05:01
$group{"A0305"} = "A305";
$base{"A0305"} = "A3";		#A305 and A1104 are same serotype
push @subtype, ("A0305");
my %a11;
$a11{"A1101"} = "HLA00043";		# A*11:01:01:01
$group{"A1101"} = "A11";
$base{"A1101"} = "A11";
my %a9;
$a9{"A2301"} = "HLA00048";		#A*23:01:01:01, Bw4
$a9{"A2304"} = "HLA01108";		# A*23:04, Bw4
$a9{"A2402"} = "HLA00050";		# A*24:02:01:01, Bw4
$a9{"A2403"} = "HLA00053";		# A*24:03:01:01, Bw4
$a9{"A2404"} = "HLA00054";		# A*24:04, non-Bw4
$a9{"A2405"} = "HLA00055";		# A*24:05:01, Bw4
$a9{"A0246"} = "HLA01274";		# A*02:46, Bw6Neg
$a9{"A2408"} = "HLA00058";		# A*24:08
$a9{"A2410"} = "HLA00060";		# A*24:10:01:01 365 bp
$a9{"A2414"} = "HLA00064";		# A*24:14:01:01 365 bp
$a9{"A2423"} = "HLA01041";		# A*24:23, included for common allele
$a9{"A2424"} = "HLA01042";		# A*24:24
$bw{"A2301"} = "Bw4"; $bw{"A2304"} = "Bw4"; 	# nothing for A246
$bw{"A2402"} = "Bw4"; $bw{"A2403"} = "Bw4"; $bw{"A2405"} = "Bw4"; $bw{"A2408"} = "Bw4"; $bw{"A2410"} = "Bw4"; $bw{"A2414"} = "Bw4";
$bw{"A2423"} = "Bw4"; $bw{"A2424"} = "Bw4";
$group{"A2301"} = "A9"; $group{"A2304"} = "A9"; $group{"A2408"} = "A9"; $group{"A2410"} = "A9"; $group{"A2414"} = "A9";
$group{"A2402"} = "A9"; $group{"A2403"} = "A9"; $group{"A2404"} = "A9"; $group{"A2405"} = "A9"; $group{"A2424"} = "A9"; $group{"A2423"} = "A9";
$group{"A0246"} = "A9";	# use residue 76 to distinguish from A2403
$base{"A2301"} = "A23"; $base{"A2304"} = "A23"; $base{"A2408"} = "A24";# $base{"A2410"} = "A2403"; $base{"A2414"} = "A24";
$base{"A2410"} = "A24"; $base{"A2414"} = "A24";
$base{"A2402"} = "A24"; # $base{"A2403"} = "A2403"; $base{"A2404"} = "A24"; $base{"A2405"} = "A24"; $base{"A2424"} = "A23"; $base{"A2423"} = "A2403";
$base{"A2403"} = "A24"; $base{"A2404"} = "A24"; $base{"A2405"} = "A24"; $base{"A2424"} = "A23"; $base{"A2423"} = "A24";
$base{"A0246"} = "A2";	# was None
push @subtype, ("A2304","A2403","A2404","A2405","A2408","A2410","A2414","A2423","A2424");
push @subtype, ("A0246");
my %a25;
$a25{"A2501"} = "HLA00071";		# A*25:01:01:01
$bw{"A2501"} = "Bw4";
$group{"A2501"} = "A25";
$base{"A2501"} = "A25";
my %a26a43;
$a26a43{"A2601"} = "HLA00073";		# A*26:01:01:01
$a26a43{"A2603"} = "HLA00075";		# A*26:03:01:01
$a26a43{"A2607"} = "HLA00079";		# A*26:07:01
$a26a43{"A2614"} = "HLA01120";		# A*26:14, maybe added
$a26a43{"A4301"} = "HLA00111";		# A*43:01
$group{"A2601"} = "A26A43"; $group{"A2603"} = "A26A43"; $group{"A2607"} = "A26A43"; $group{"A4301"} = "A26A43";
$group{"A2614"} = "A26A43";
$base{"A2601"} = "A26"; $base{"A2603"} = "A26"; $base{"A2607"} = "A26"; $base{"A4301"} = "A43";
$base{"A2614"} = "A26";
push @subtype, ("A2603","A2607","A2614");
my %a29;
$a29{"A2901"} = "HLA00085";		# A*29:01:01:01
$a29{"A2902"} = "HLA00086";		# A*29:02:01:01 365 bp
$group{"A2901"} = "A29"; $group{"A2902"} = "A29";
$base{"A2901"} = "A29"; $base{"A2902"} = "A29";
push @subtype, ("A2902");
my %a30;
$a30{"A3001"} = "HLA00089";		# A*30:01:01:01
$a30{"A3002"} = "HLA00090";		# A*30:02:01:01
$a30{"A3007"} = "HLA00095";		# A*30:07, included for common allele
#$a30{"A2914"} = "HLA02256";		# A*29:14
$group{"A3001"} = "A30"; $group{"A3002"} = "A30"; $group{"A3007"} = "A30";# $group{"A2914"} = "A30";
$base{"A3001"} = "A30"; $base{"A3002"} = "A30"; $base{"A3007"} = "A30";# $base{"A2914"} = "A30";
push @subtype, ("A3002", "A3007");
my %a31;
$a31{"A3101"} = "HLA00097";		# A*31:01:02:01
$a31{"A3102"} = "HLA00098";		# A*31:02:01 365 bp
$group{"A3101"} = "A31";$group{"A3102"} = "A31";
$base{"A3101"} = "A31";$base{"A3102"} = "A31";
push @subtype, ("A3102");
my %a32;
$a32{"A3201"} = "HLA00101";		# A*32:01:01:01
# removed A2309, but added A3204 on Nov 30 2020
$a32{"A3204"} = "HLA01045";		# A*32:04
#$a32{"A2309"} = "HLA01571";		# A*23:09, similar to A32:13
$bw{"A3201"} = "Bw4";
$bw{"A3204"} = "Bw4";
$group{"A3201"} = "A32";
$group{"A3204"} = "A32";
$base{"A3201"} = "A32";
$base{"A3204"} = "A32";
push @subtype,("A3204");
my %a33;
$a33{"A3301"} = "HLA00104";		# A*33:01:01:01
$a33{"A3303"} = "HLA00106";		# A*33:03:01:01
$group{"A3301"} = "A33"; $group{"A3303"} = "A33";
$base{"A3301"} = "A33"; $base{"A3303"} = "A33";
push @subtype, ("A3303");
my %a3313;	# created to avoid artificial outlier
$a3313{"A3313"} = "HLA02981";		# A*33:13
#$a3313{"A3311"} = "HLA02918";		# A*33:11
$group{"A3313"} = "A3313"; #$group{"A3311"} = "A3313";
$base{"A3313"} = "A33"; #$base{"A3311"} = "A33";
push @subtype, ("A3313");
my %a34;
$a34{"A3401"} = "HLA00108";		# A*34:01:01:01
$a34{"A3402"} = "HLA00109";		# A*34:02:01:01
$group{"A3401"} = "A34"; $group{"A3402"} = "A34";
$base{"A3401"} = "A34"; $base{"A3402"} = "A34";
push @subtype, ("A3402");
my %a36;
$a36{"A3601"} = "HLA00110";		# A*36:01
$group{"A3601"} = "A36";
$base{"A3601"} = "A36";
my %a66;
$a66{"A6601"} = "HLA00112";		# A*66:01:01:01
$a66{"A6602"} = "HLA00113";		# A*66:02
$group{"A6601"} = "A66"; $group{"A6602"} = "A66";
$base{"A6601"} = "A66"; $base{"A6602"} = "A66";
push @subtype, ("A6602");
my %a68;
$a68{"A6801"} = "HLA00115";		# A*68:01:01:01
$a68{"A6810"} = "HLA00972";		# A*68:10, included for common allele A*68:10,S,A246,None,None,C,WD,,A246_65
$a68{"A6813"} = "HLA01047";		# A*68:13:01
$group{"A6801"} = "A68"; $group{"A6810"} = "A68"; $group{"A6813"} = "A68";
$base{"A6801"} = "A68"; $base{"A6810"} = "A68"; $base{"A6813"} = "A68";
my %a6836;
$a6836{"A6836"} = "HLA02701";		# A*68:36
$group{"A6836"} = "A6836";
$base{"A6836"} = "A68";
$bw{"A6836"} = "Bw4"; 
push @subtype, ("A6810","A6813","A6836");

my %a69;
$a69{"A6901"} = "HLA00126";		# A*69:01:01:01
$group{"A6901"} = "A69";
$base{"A6901"} = "A69";
my %a74;
$a74{"A7401"} = "HLA00127";		# A*74:01:01:01
$a74{"A0265"} = "HLA01778";		# A*02:65
$a74{"A3308"} = "HLA02250";		# A*33:08
$group{"A7401"} = "A74"; $group{"A0265"} = "A74"; $group{"A3308"} = "A74";
$base{"A7401"} = "A74"; $base{"A0265"} = "A2"; $base{"A3308"} = "A33";
push @subtype, ("A0265", "A3308");
my %a80;
$a80{"A8001"} = "HLA00130";		# A*80:01:01:01
$group{"A8001"} = "A80";
$base{"A8001"} = "A80";


sub HLAA {
	my $gene = "A";
	return $gene;
}

sub HLAA_LEADER {
	my $leader = 23;		# A specific
	return $leader;
}

sub GROUP {
	my $group_ref = \%group;
	return $group_ref;
}

sub BASE {	# closest WHO accepted antigen
	my $base_ref = \%base;
	return $base_ref;
}

sub BASETYPE {		# populate WHO accepted type
	my %unique;
	my @basetype;
	foreach my $value (sort values %base ) {	# value
		unless ( exists $unique{ $value } ) {
			push @basetype, $value;
			$unique{ $value } = 0;
		}
	}
	my $basetype_ref = \@basetype;
	return $basetype_ref;		#######
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
		if (( $base{ $base } eq "A23" ) || ( $base{ $base } eq "A24" )) {
			$broad{ $base } = "A9";
		}
		elsif (( $base{ $base } eq "A68" ) || ( $base{ $base } eq "A69" )) {
			$broad{ $base } = "A28";
		}
		elsif (( $base{ $base } eq "A25" ) || ( $base{ $base } eq "A26" ) || ( $base{ $base } eq "A34" ) || ( $base{ $base } eq "A66" )) {
			$broad{ $base } = "A10";
		}
		elsif (( $base{ $base } eq "A29" ) || ( $base{ $base } eq "A30" ) || ( $base{ $base } eq "A31" ) || ( $base{ $base } eq "A32" )
		 || ( $base{ $base } eq "A33" ) || ( $base{ $base } eq "A74" )) {
			$broad{ $base } = "A19";
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
	push @combined, @a1; 
	push @combined, @a2; 
	push @combined, @a219;
	push @combined, @a3; 
	push @combined, @a305; 
	push @combined, @a11;
	push @combined, @a9; 
	push @combined, @a25; 
	push @combined, @a26a43; 
	push @combined, @a29; 
	push @combined, @a30; 
	push @combined, @a31; 
	push @combined, @a32; 
	push @combined, @a33; 
	push @combined, @a3313; 
	push @combined, @a34; 
	push @combined, @a36; 
	push @combined, @a66; 
	push @combined, @a68; 
	push @combined, @a6836; 
	push @combined, @a69; 
	push @combined, @a74; 
	push @combined, @a80; 
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
	if ( $serotype eq "A1" ) {
		@residues = @a1; 
	}
	elsif ( $serotype eq "A2" ) {
		@residues = @a2;
	}
	elsif ( $serotype eq "A219" ) {
		@residues = @a219;
	}
	elsif ( $serotype eq "A3" ) {
		@residues = @a3;
	}
	elsif ( $serotype eq "A305" ) {
		@residues = @a305;
	}
	elsif ( $serotype eq "A11" ) {
		@residues = @a11;
	}
	elsif ( $serotype eq "A9" ) {
		@residues = @a9;
	}
	elsif ( $serotype eq "A25" ) {
		@residues = @a25;
	}
	elsif ( $serotype eq "A26A43" ) {
		@residues = @a26a43;
	}
	elsif ( $serotype eq "A29" ) {
		@residues = @a29;
	}
	elsif ( $serotype eq "A30" ) {
		@residues = @a30;
	}
	elsif ( $serotype eq "A31" ) {
		@residues = @a31;
	}
	elsif ( $serotype eq "A32" ) {
		@residues = @a32;
	}
	elsif ( $serotype eq "A33" ) {
		@residues = @a33;
	}
	elsif ( $serotype eq "A3313" ) {
		@residues = @a3313;
	}
	elsif ( $serotype eq "A34" ) {
		@residues = @a34;
	}
	elsif ( $serotype eq "A36" ) {
		@residues = @a36;
	}
	elsif ( $serotype eq "A66" ) {
		@residues = @a66;
	}
	elsif ( $serotype eq "A68" ) {
		@residues = @a68;
	}
	elsif ( $serotype eq "A6836" ) {
		@residues = @a6836;
	}
	elsif ( $serotype eq "A69" ) {
		@residues = @a69;
	}
	elsif ( $serotype eq "A74" ) {
		@residues = @a74;
	}
	elsif ( $serotype eq "A80" ) {
		@residues = @a80;
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

	if ( $serotype eq "A1" ) {
		%ref = %a1; 
	}
	elsif ( $serotype eq "A2" ) {
		%ref = %a2;
	}
	elsif ( $serotype eq "A219" ) {
		%ref = %a219;
	}
	elsif ( $serotype eq "A3" ) {
		%ref = %a3;
	}
	elsif ( $serotype eq "A305" ) {
		%ref = %a305;
	}
	elsif ( $serotype eq "A11" ) {
		%ref = %a11;
	}
	elsif ( $serotype eq "A9" ) {
		%ref = %a9;
	}
	elsif ( $serotype eq "A25" ) {
		%ref = %a25;
	}
	elsif ( $serotype eq "A26A43" ) {
		%ref = %a26a43;
	}
	elsif ( $serotype eq "A29" ) {
		%ref = %a29;
	}
	elsif ( $serotype eq "A30" ) {
		%ref = %a30;
	}
	elsif ( $serotype eq "A31" ) {
		%ref = %a31;
	}
	elsif ( $serotype eq "A32" ) {
		%ref = %a32;
	}
	elsif ( $serotype eq "A33" ) {
		%ref = %a33;
	}
	elsif ( $serotype eq "A3313" ) {
		%ref = %a3313;
	}
	elsif ( $serotype eq "A34" ) {
		%ref = %a34;
	}
	elsif ( $serotype eq "A36" ) {
		%ref = %a36;
	}
	elsif ( $serotype eq "A66" ) {
		%ref = %a66;
	}
	elsif ( $serotype eq "A68" ) {
		%ref = %a68;
	}
	elsif ( $serotype eq "A6836" ) {
		%ref = %a6836;
	}
	elsif ( $serotype eq "A69" ) {
		%ref = %a69;
	}
	elsif ( $serotype eq "A74" ) {
		%ref = %a74;
	}
	elsif ( $serotype eq "A80" ) {
		%ref = %a80;
	}
	else {	# all together
		%ref = (%a1,%a2,%a219,%a3,%a305,%a11,%a9,%a25,%a26a43,%a29,%a30,%a31,%a32,%a33,%a34,%a36,%a66,%a68,%a6836,%a69,%a74,%a80,%a219,%a3313);
	}
	
	return $ref_ref;
}

# 
sub SERO {
	my @sero;
	my %ref = (%a1,%a2,%a219,%a3,%a305,%a11,%a9,%a25,%a26a43,%a29,%a30,%a31,%a32,%a33,%a34,%a36,%a66,%a68,%a6836,%a69,%a74,%a80,%a219,%a3313);
	my @tmp = sort keys %ref;
	for ( my $index = 0; $index < scalar @tmp; $index++ ) {
		$sero[0][$index] = $tmp[$index];	# populate serotype
	}
	for ( my $index = 0; $index < scalar @subtype; $index++ ) {
		$sero[1][$index] = $subtype[$index];	# populate subtype
	}
	my $sero_ref = \@sero;
	return $sero_ref;
}

sub KEY {
	my %tmp = (%a1,%a2,%a219,%a3,%a305,%a11,%a9,%a25,%a26a43,%a29,%a30,%a31,%a32,%a33,%a34,%a36,%a66,%a68,%a6836,%a69,%a74,%a80,%a219,%a3313);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key eq "A0101" ) {
			$ref{$key} = "A\\*01";
		}
		elsif (( $key eq "A0201" ) || ( $key eq "A0202" ) || ( $key eq "A0203" ) || ( $key eq "A0208" ) || ( $key eq "A0210" ) ||
		( $key eq "A0211" ) || ( $key eq "A0218" ) || ( $key eq "A0265" ) || ( $key eq "A0246" ) || ( $key eq "A0285" ) || ( $key eq "A0256" )
		|| ($key eq "A0219") || ($key eq "A0244")) {
			$ref{$key} = "A\\*02";
		}
		elsif (( $key eq "A0301" ) || ( $key eq "A0305" ) || ( $key eq "A0323" )) {
			$ref{$key} = "A\\*03";
		}
		elsif ( $key eq "A1101" ) {
			$ref{$key} = "A\\*11";
		}
		elsif (( $key eq "A2301" ) || ( $key eq "A2304" ) || ( $key eq "A2309" )) {
			$ref{$key} = "A\\*23";
		}
		elsif (( $key eq "A2402" ) || ( $key eq "A2403" ) || ( $key eq "A2404" ) || ( $key eq "A2405" ) || ( $key eq "A2408" ) ||
			( $key eq "A2423" ) || ( $key eq "A2424" )) {
			$ref{$key} = "A\\*24";
		}
		elsif ( $key eq "A2501" ) {
			$ref{$key} = "A\\*25";
		}
		elsif (( $key eq "A2601" ) || ( $key eq "A2603" ) || ( $key eq "A2607" ) || ( $key eq "A2614")) {
			$ref{$key} = "A\\*26";
		}
		elsif (( $key eq "A2901" ) || ( $key eq "A2902" )) {
			$ref{$key} = "A\\*29";
		}
		elsif (( $key eq "A3001" ) || ( $key eq "A3002" ) || ( $key eq "A3007" )) {
			$ref{$key} = "A\\*30";
		}
		elsif ( $key eq "A3101" ) {
			$ref{$key} = "A\\*31";
		}
		elsif (( $key eq "A3201" ) || ( $key eq "A3204" )) {
			$ref{$key} = "A\\*32";
		}
		elsif (( $key eq "A3301" ) || ( $key eq "A3303" ) || ( $key eq "A3308") || ( $key eq "A3313" )) {
			$ref{$key} = "A\\*33";
		}
		elsif (( $key eq "A3401" ) || ( $key eq "A3402" )) {
			$ref{$key} = "A\\*34";
		}
		elsif ( $key eq "A3601" ) {
			$ref{$key} = "A\\*36";
		}
		elsif ( $key eq "A4301" ) {
			$ref{$key} = "A\\*43";
		}
		elsif (( $key eq "A6601" ) || ( $key eq "A6602" )) {
			$ref{$key} = "A\\*66";
		}
		elsif (( $key eq "A6801" ) || ( $key eq "A6810" ) || ( $key eq "A6836")) {
			$ref{$key} = "A\\*68";
		}
		elsif ( $key eq "A6901" ) {
			$ref{$key} = "A\\*69";
		}
		elsif ( $key eq "A7401" ) {
			$ref{$key} = "A\\*74";
		}
		else {
			$ref{$key} = "A\\*80";
		}
	}
	return $key_ref;

}

sub BW {
	my $bw_ref = \%bw;
	return $bw_ref;
}

sub PARTIAL {		# partial sequence
	my %partial;
	my $partial_ref = \%partial;
	my $seq = "X" x 25;	#change the number of missing nucleotide
#	$partial{ "A-2424" } = $seq;	# partial 3.54.0, full 3.62.0
	$partial{ "A0323" } = $seq;	# partial 3.54.0
	$partial{ "A3313" } = $seq;	# partial 3.54.0
	$partial{ "A*11:14" } = "X" x 16;	# partial 3.54.0
	$partial{ "general" } = $seq;
		
	return $partial_ref;
}

# Do I need this?
#sub WHO {
#	my %who;
#	my $whotype_ref = \%who;
#	$who{"A0101"} = "A1"; $who{"A0201"} = "A2"; $who{"A0203"} = "A203"; $who{"A0210"} = "A210"; $who{"A0301"} = "A3";
#	$who{"A1101"} = "A11"; $who{"A2301"} = "A23"; $who{"A2402"} = "A24"; $who{"A2403"} = "A2403"; $who{"A2501"} = "A25";
#	$who{"A2601"} = "A26"; $who{"A2901"} = "A29"; $who{"A3001"} = "A30"; $who{"A3101"} = "A31"; $who{"A3201"} = "A32";
#	$who{"A3301"} = "A33"; $who{"A3401"} = "A34"; $who{"A3601"} = "A36"; $who{"A4301"} = "A43";
#	$who{"A6601"} = "A66"; $who{"A6801"} = "A68"; $who{"A6901"} = "A69"; $who{"A7401"} = "A74"; $who{"A8001"} = "A80";

#	return $whotype_ref;
#}

sub KNOWN_CROSS {	# trick to make SEROTYPE to FULL
	my %known_cross;
	my $known_cross_ref = \%known_cross;
	$known_cross{ "NOTHING" } = 0;
	return $known_cross_ref;
}

1;
