#!/usr/bin/perl -w

# Author: Kazutoyo Osoegawa, Ph.D.
# Developed at Stanford Blood Center
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# module: HLAB_INFO.pm 
# This module was developed to convert HLA allele to HLA serotype
# last modified and documented on October 18 2022
# Included Bw6, Negative and Bw4

package HLAB_INFO;
use strict;

# define key residues
my @b5 = (45,62,63,82,83,103,163,167,171);	#Bw4: no 76  residue 103 is included to distinguish B5102 from B53
my @b7 = (45,63,69,70,71,76,82,83,147,163,171,178);	#Bw6 76,82,83
my @b713 = (45,63,69,70,71,76,82,83,103,147,163,171,178);	#Bw6 76,82,83, 103 was added to eliminate HLA-C alleles
my @b8 = (45,63,69,71,76,82,83,163,177,180);	#Bw6 76,82,83
my @b802 = (45,63,69,71,82,83,163,177,180);	#Bw4 82,83
my @b44 = (45,62,63,82,83,163,167);	#Bw4, added 62 to eliminate B5707
my @b45 = (45,63,76,82,83,163,167);	#Bw6
my @b13 = (45,63,82,83,145,163,178);	#Bw4
my @b14 = (11,45,67,76,82,83,158,163,171);	#Bw6, 158 was added to eliminate B*39:48, 158 did not affect B14
my @b15 = (45,46,62,63,69,70,76,82,83,163,167,171,178,180);	#Bw6, 178 was added to separate B1552 from B703
my @b63 = (45,46,62,63,69,70,76,82,83,103,163,167,171,178,180);	#Bw6, 178 was added to separate B1552 from B703
my @b1523 = (45,46,62,63,67,69,70,76,82,83,163,167,178);	
my @b17 = (45,46,62,63,69,70,83,163);	#removed 82
my @b38 = (45,63,67,71,82,83,158,163,171);	#Bw4, 171 was added to distinguish from B-1809
my @b39 = (45,63,67,71,76,82,83,158,163);	#Bw6
my @b18 = (45,63,76,82,83,163,171);	#Bw6
my @b1809 = (45,63,82,83,163,171);	#Bw4
my @b49 = (45,46,63,82,83,103,163,167,178);	#Bw4
my @b50 = (45,46,63,76,82,83,103,163,167,178);	#Bw6, use 103 instead of 131
my @b22 = (45,63,67,69,76,83,158,163,167,177,178,180);	#Bw6, removed 82
my @b27 = (45,63,70,71,82,83,178);		#Bw4, 178 to separate B48
my @b2708 = (45,63,70,71,76,82,83,143,163,178);	#Bw6, 143 was included to separete B2712 from B4804, 178 for B48
my @b35 = (45,63,70,71,76,82,83,163,171);	#Bw6
my @b37 = (45,62,63,70,76,82,83,163);	#Bw4
my @b3705 = (45,62,63,70,76,82,83,163);	#Bw6
my @b40 = (45,46,76,82,83,143,163,178);	#Bw6
my @b48 = (45,63,70,71,76,82,83,143,163,167,177,178,180);	#Bw6
my @b41 = (45,63,70,71,76,83,163,167,177,178,180);		# remoed 82
my @b42 = (45,67,69,70,71,76,82,83,163,177,178,180);	#Bw6
my @b46 = (45,46,62,63,69,76,82,83,163,167);	#KC1, Negative, 45,46 is required to separate from Cw3
my @b47 = (45,63,70,71,82,83,143,163,177);	#Bw4
my @b53 = (45,63,70,71,82,83,103,163,171);	# residue 103 is included to distinguish B5102 from B53
my @b59 = (45,63,67,71,82,83,163,171,177);	#Bw4
my @b67 = (45,67,71,76,82,83,103,158,163);	#Residue 103 was required to be HLA-B specific, Bw6
my @b73 = (45,63,67,70,76,82,83,163,178);	#KC1, Negative
my @b78 = (45,76,82,83,163,171);	#Bw6
my @b81 = (45,69,70,71,76,82,83,147,163,178);	#Bw6
my @b82b83 = (45,63,67,69,70,71,76,82,83,167,178);	#Bw6, used residue 103 instead of 162 to eliminate residue 162 from FULL
my @extra = (127);

my %bw;
my %group;	# same residue group
my %base;	# convert to known serotype
my @subtype;
# capture reference allele
my %b5;	#Bw4
$b5{"B-5101"} = "HLA00344";	#B*51:01:01:01
$b5{"B-5102"} = "HLA00346";	#B*51:02:01:01
$b5{"B-5103"} = "HLA00348";	#B*51:03
$b5{"B-5201"} = "HLA00362";	#B*52:01:01:01
#$b5{"B-5203"} = "HLA01295";	# B*52:03
$b5{"B-4406"} = "HLA00323";	# B*44:06, added on March 10 2022
$bw{"B-5101"} = "Bw4"; $bw{"B-5102"} = "Bw4"; $bw{"B-5103"} = "Bw4"; $bw{"B-5201"} = "Bw4";# $bw{"B-5203"} = "Bw4";
$bw{"B-4406"} = "Bw4";
$group{"B-5101"} = "B5"; $group{"B-5102"} = "B5"; $group{"B-5103"} = "B5"; $group{"B-5201"} = "B5";# $group{"B-5203"} = "B5";
$group{"B-4406"} = "B5";
$base{"B-5101"} = "B51"; $base{"B-5102"} = "B5102"; $base{"B-5103"} = "B5103"; $base{"B-5201"} = "B52";# $base{"B-5203"} = "B52";
$base{"B-4406"} = "B44";
push @subtype, ("B-5102","B-5103");
push @subtype, ("B-4406");

my %b7;
$b7{"B-0702"} = "HLA00132";		#B*07:02:01:01, Bw6
$b7{"B-0703"} = "HLA00135";	#B*07:03, Bw6
$b7{"B-0715"} = "HLA01050";	#B*07:15, Bw6 negative
$b7{"B-0736"} = "HLA01808";	# B*07:36, Bw4
$bw{"B-0702"} = "Bw6"; $bw{"B-0703"} = "Bw6"; $bw{"B-0715"} = "Negative"; $bw{"B-0736"} = "Bw4"; 
$group{"B-0702"} = "B7"; $group{"B-0703"} = "B7"; $group{"B-0715"} = "B7"; $group{"B-0736"} = "B7"; 
$base{"B-0702"} = "B7"; $base{"B-0703"} = "B703"; $base{"B-0715"} = "B7"; $base{"B-0736"} = "B7"; 

my %b713;
$b713{"B-0713"} = "HLA00145";	# B*07:13, ref seq was partial 3.38.0, but full length becomes available 3.39.0
$bw{"B-0713"} = "Negative";
$group{"B-0713"} = "B713";
$base{"B-0713"} = "None";
push @subtype, ("B-0703","B-0713","B-0715","B-0736");

my %b8;
$b8{"B-0801"} = "HLA00146";	#B*08:01:01:01, Bw6
my %b802;
$b802{"B-0802"} = "HLA00147";	#B*08:02, Bw4
$bw{"B-0801"} = "Bw6"; $bw{"B-0802"} = "Bw4";
$group{"B-0801"} = "B8"; $group{"B-0802"} = "B802";
$base{"B-0801"} = "B8"; $base{"B-0802"} = "B8";
push @subtype, ("B-0802");

#my %b12;
my %b44;
$b44{"B-4402"} = "HLA00318";	#B*44:02:01:01, Bw4
$b44{"B-4404"} = "HLA00321";	#B*44:04
$b44{"B-4408"} = "HLA00325";	# B*44:08
$bw{"B-4402"} = "Bw4"; $bw{"B-4404"} = "Bw4"; $bw{"B-4408"} = "Bw4";
$group{"B-4402"} = "B44"; $group{"B-4404"} = "B44"; $group{"B-4408"} = "B44";
$base{"B-4402"} = "B44"; $base{"B-4404"} = "B44"; $base{"B-4408"} = "B44";
push @subtype, ("B-4404", "B-4408");

my %b45;
$b45{"B-4501"} = "HLA00329";	#B*45:01:01:01, Bw6
$bw{"B-4501"} = "Bw6";
$group{"B-4501"} = "B45";
$base{"B-4501"} = "B45";

my %b13;
$b13{"B-1302"} = "HLA00153";	#B*13:02:01:01, Bw4
$bw{"B-1302"} = "Bw4";# $bw{"B1309"} = "Bw6";
$group{"B-1302"} = "B13";# $group{"B1309"} = "B13";
$base{"B-1302"} = "B13";# $base{"B1309"} = "B13";
#push @subtype, ("B1309");

my %b14;
$b14{"B-1401"} = "HLA00157";	#B*14:01:01:01, Bw6
$b14{"B-1402"} = "HLA00158";	#B*14:02:01:01, Bw6
$bw{"B-1401"} = "Bw6"; $bw{"B-1402"} = "Bw6";
$group{"B-1401"} = "B14"; $group{"B-1402"} = "B14";
$base{"B-1401"} = "B64"; $base{"B-1402"} = "B65";
push @subtype, ("B-1402");

my %b15;	# this name does not reflect allele name, but convenience to separate Bw6
$b15{"B-1501"} = "HLA00162";	#15:01:01:01, Bw6
$b15{"B-1502"} = "HLA00165";	#15:02:01:01, Bw6
$b15{"B-1512"} = "HLA00175";	#15:12:01, Bw6
$b15{"B-1514"} = "HLA00177";	#15:14, Bw6
$b15{"B-1510"} = "HLA00173";	#15:10:01, Bw6
$b15{"B-1503"} = "HLA00166";	#15:03:01:01, Bw6
$b15{"B-1540"} = "HLA00203";	# B*15:40:01
$b15{"B-1552"} = "HLA01053";	# B*15:52
$b15{"B-1542"} = "HLA00205";	# B*15:42
$b15{"B-1524"} = "HLA00187";	#B*15:24:01, Bw4, B62Bw4
$b15{"B-1513"} = "HLA00176";	#15:13:01, Bw4
$b15{"B-1537"} = "HLA00200";	# B*15:37 362 bp
$b15{"B-1538"} = "HLA00201";	# B*15:38:01 362 bp

$bw{"B-1501"} = "Bw6"; $bw{"B-1502"} = "Bw6"; $bw{"B-1512"} = "Bw6"; $bw{"B-1514"} = "Bw6"; $bw{"B-1510"} = "Bw6"; $bw{"B-1503"} = "Bw6";
$bw{"B-1540"} = "Bw6"; $bw{"B-1552"} = "Bw6"; $bw{"B-1542"} = "Bw6"; $bw{"B-1524"} = "Bw4";  $bw{"B-1513"} = "Bw4";
$bw{"B-1537"} = "Bw6"; $bw{"B-1538"} = "Bw6";# $bw{"B-1520"} = "Bw6"; $bw{"B-4802"} = "Bw6";
$group{"B-1501"} = "B15"; $group{"B-1502"} = "B15"; $group{"B-1512"} = "B15"; $group{"B-1514"} = "B15"; $group{"B-1510"} = "B15"; $group{"B-1503"} = "B15"; 
$group{"B-1540"} = "B15"; $group{"B-1552"} = "B15"; $group{"B-1542"} = "B15"; $group{"B-1524"} = "B15"; $group{"B-1513"} = "B15";
$group{"B-1537"} = "B15"; $group{"B-1538"} = "B15";# $group{"B-1520"} = "B15"; $bw{"B-4802"} = "B15";
$base{"B-1501"} = "B62"; $base{"B-1502"} = "B75"; $base{"B-1512"} = "B76"; $base{"B-1514"} = "B76"; $base{"B-1510"} = "B71"; $base{"B-1503"} = "B72"; 
$base{"B-1540"} = "B62"; $base{"B-1552"} = "B71"; $base{"B-1542"} = "B62"; $base{"B-1524"} = "B62"; $base{"B-1513"} = "B77";
$base{"B-1537"} = "B71"; $base{"B-1538"} = "B62";

my%b63;
$b63{"B-1517"} = "HLA00180";	#15:17:01:01, Bw4  
$b63{"B-1516"} = "HLA00179";	# B*15:16:01:01
$bw{"B-1517"} = "Bw4"; $bw{"B-1516"} = "Bw4";
$group{"B-1517"} = "B63"; $group{"B-1516"} = "B63";
$base{"B-1517"} = "B63"; $base{"B-1516"} = "B63";

my %b1523;
$b1523{"B-1523"} = "HLA00186";	# B*15:23
$bw{"B-1523"} = "Bw4"; $group{"B-1523"} = "B1523"; $base{"B-1523"} = "B71";
my %b17;	# this name does not reflect allele name, but convenience to separate Bw4
$b17{"B-5701"} = "HLA00381";	#57:01:01:01, Bw4
$b17{"B-5801"} = "HLA00386";	#58:01:01:01, Bw4
$bw{"B-5701"} = "Bw4"; $bw{"B-5801"} = "Bw4";
$group{"B-5701"} = "B17"; $group{"B-5801"} = "B17";
$base{"B-5701"} = "B57"; $base{"B-5801"} = "B58";

push @subtype, ("B-1502","B-1512","B-1514","B-1510","B-1503","B-1540","B-1552","B-1542","B-1523","B-1524","B-1517","B-1513","B-1537","B-1538", "B-1516");

my %b38;	#B16
$b38{"B-3801"} = "HLA00267";	#B*38:01:01:01, Bw4
$b38{"B-3803"} = "HLA00270";	# B*38:03
$b38{"B-3806"} = "HLA01265";	# B*38:06
$bw{"B-3801"} = "Bw4"; $bw{"B-3803"} = "Bw4"; $bw{"B-3806"} = "Bw4";
$group{"B-3801"} = "B38"; $group{"B-3803"} = "B38"; $group{"B-3806"} = "B38";
$base{"B-3801"} = "B38"; $base{"B-3803"} = "B38"; $base{"B-3806"} = "B38";
push @subtype, ("B-3803","B-3806");
my %b39;	#B16
$b39{"B-3901"} = "HLA00271";	#B*39:01:01:01, Bw6
$b39{"B-3902"} = "HLA00274";	#B*39:02:01, Bw6
$b39{"B-3910"} = "HLA00284";	# B*39:10:01, Bw6
$bw{"B-3901"} = "Bw6"; $bw{"B-3902"} = "Bw6"; $bw{"B-3910"} = "Bw6";
$group{"B-3901"} = "B39"; $group{"B-3902"} = "B39"; $group{"B-3910"} = "B39";
$base{"B-3901"} = "B39"; $base{"B-3902"} = "B3902"; $base{"B-3910"} = "B39";
push @subtype, ("B-3902","B-3910");

my %b18;
$b18{"B-1801"} = "HLA00213";	#B*18:01:01:01, Bw6
$b18{"B-1806"} = "HLA00218";	#B*18:06, Negative
my %b1809;
$b1809{"B-1809"} = "HLA01131";	#B*18:09, Bw4
$bw{"B-1801"} = "Bw6"; $bw{"B-1806"} = "Negative"; $bw{"B-1809"} = "Bw4";
$group{"B-1801"} = "B18"; $group{"B-1806"} = "B18"; $group{"B-1809"} = "B1809";
$base{"B-1801"} = "B18"; $base{"B-1806"} = "B18"; $base{"B-1809"} = "B18";
push @subtype, ("B-1806","B-1809");


my %b49;	# B21, Bw4
$b49{"B-4901"} = "HLA00340";	#B*49:01:01:01, Bw4
my %b50;	#B21, Bw6
$b50{"B-4005"} = "HLA00296";	#B*40:05:01:01, Bw6
$b50{"B-5001"} = "HLA00341";	#B*50:01:01:01, Bw6
$bw{"B-4901"} = "Bw4"; $bw{"B-4005"} = "Bw6"; $bw{"B-5001"} = "Bw6";
$group{"B-4901"} = "B49"; $group{"B-4005"} = "B50"; $group{"B-5001"} = "B50";
$base{"B-4901"} = "B49"; $base{"B-4005"} = "B4005"; $base{"B-5001"} = "B50";
push @subtype, ("B-4005");

my %b22;	# Bw6
$b22{"B-5401"} = "HLA00367";
$b22{"B-5501"} = "HLA00368";
$b22{"B-5601"} = "HLA00376";	#B*56:01:01:01
$bw{"B-5401"} = "Bw6"; $bw{"B-5501"} = "Bw6"; $bw{"B-5601"} = "Bw6";
$group{"B-5401"} = "B22";# $group{"B5412"} = "B22";
$group{"B-5501"} = "B22";# $group{"B5503"} = "B22";
$group{"B-5601"} = "B22";
$base{"B-5401"} = "B54";# $base{"B5412"} = "B54";
$base{"B-5501"} = "B55";# $base{"B5503"} = "B55";
$base{"B-5601"} = "B56";

my %b27;
$b27{"B-2705"} = "HLA00225";	# B*27:05:02:01: used more common allele than HLA00220";	#Bw4
$bw{"B-2705"} = "Bw4";

my %b2708;
$b2708{"B-2708"} = "HLA00229";	#Bw6
$b2708{"B-2712"} = "HLA00233";	# B*27:12:01:01, Bw6
$bw{"B-2708"} = "Bw6"; $bw{"B-2712"} = "Bw6";
$group{"B-2705"} = "B27"; $group{"B-2708"} = "B2708";
$group{"B-2712"} = "B2708";	# ATTENTION, use B2708 residues, problem is residue 67
$base{"B-2705"} = "B27"; $base{"B-2708"} = "B2708";
$base{"B-2712"} = "B2708";	# B2712 is closer to B4804 than B2708
push @subtype, ("B-2708","B-2712"); 
my %b35;
$b35{"B-3501"} = "HLA00237";	#B*35:01:01:01, Bw6
$b35{"B-3510"} = "HLA00247";	#B*35:10, Bw6, residue 63
$b35{"B-3515"} = "HLA00252";	#B*35:15:01, Bw6, residue 63
$b35{"B-3519"} = "HLA00256";	#B*35:19, Bw6, residue 45
$bw{"B-3501"} = "Bw6"; $bw{"B-3510"} = "Bw6"; $bw{"B-3515"} = "Bw6"; $bw{"B-3519"} = "Bw6";# $bw{"B3574"} = "Negative";
$group{"B-3501"} = "B35"; $group{"B-3510"} = "B35"; $group{"B-3515"} = "B35"; $group{"B-3519"} = "B35";# $group{"B3574"} = "B35";
$base{"B-3501"} = "B35"; $base{"B-3510"} = "B35"; $base{"B-3515"} = "B35"; $base{"B-3519"} = "B35";# $base{"B3574"} = "B35";
push @subtype, ("B-3510","B-3515","B-3519");
my %b37;
$b37{"B-3701"} = "HLA00265";	#B*37:01:01:01, Bw4
$b37{"B-3702"} = "HLA00266";	# B*37:02, Bw4
my %b3705;
$b3705{"B-3705"} = "HLA01359";	#B*37:05, Bw6
$bw{"B-3701"} = "Bw4"; $bw{"B-3702"} = "Bw4"; $bw{"B-3705"} = "Bw6";
$group{"B-3701"} = "B37"; $group{"B-3702"} = "B37"; $group{"B-3705"} = "B3705";
$base{"B-3701"} = "B37"; $base{"B-3702"} = "B37"; $base{"B-3705"} = "B37";
push @subtype, ("B-3702","B-3705"); 
my %b40;
$b40{"B-4001"} = "HLA00291";	#B*40:01:01, Bw6
$b40{"B-4002"} = "HLA00293";	#B*40:02:01:01, Bw6
$b40{"B-4016"} = "HLA00307";	#B*40:16:01:01
$b40{"B-4021"} = "HLA00983";	# B*40:21
$bw{"B-4001"} = "Bw6"; $bw{"B-4002"} = "Bw6"; $bw{"B-4016"} = "Bw6"; $bw{"B-4021"} = "Bw6";
$group{"B-4001"} = "B40"; $group{"B-4002"} = "B40"; $group{"B-4016"} = "B40"; $group{"B-4021"} = "B40";
$base{"B-4001"} = "B60"; $base{"B-4002"} = "B61"; $base{"B-4016"} = "B60"; $base{"B-4021"} = "B60";
push @subtype, ("B-4002","B-4016","B-4021");	#,"B4050"
my %b48;
$b48{"B-4801"} = "HLA00335";	#B*48:01:01:01, Bw6
$b48{"B-4804"} = "HLA00338";	#B*48:04:01:01, Bw6
$b48{"B-4805"} = "HLA00339";	# B*48:05
$bw{"B-4801"} = "Bw6"; $bw{"B-4804"} = "Bw6"; $bw{"B-4805"} = "Bw6";
$group{"B-4801"} = "B48"; $group{"B-4804"} = "B48"; $group{"B-4805"} = "B48";
$base{"B-4801"} = "B48"; $base{"B-4804"} = "B48"; $base{"B-4805"} = "B48";
push @subtype, ("B-4804","B-4805");

my %b41;
$b41{"B-4101"} = "HLA00312";	#B*41:01:01:01, Bw6
$group{"B-4101"} = "B41";
$base{"B-4101"} = "B41";
$bw{"B-4101"} = "Bw6";

my %b42;
$b42{"B-4201"} = "HLA00315";	#B*42:01:01:01, Bw6
$bw{"B-4201"} = "Bw6";
$group{"B-4201"} = "B42";
$base{"B-4201"} = "B42";

my %b46;
$b46{"B-4601"} = "HLA00331";	#B*46:01:01:01
$bw{"B-4601"} = "Negative";# $bw{"B4640"} = "Negative";
$group{"B-4601"} = "B46";# $group{"B4640"} = "B46";
$base{"B-4601"} = "B46";# $base{"B4640"} = "B46";

my %b47;
$b47{"B-4701"} = "HLA01437";	#B*47:01:01:02,  B61Bw4 belongs to B47
$b47{"B-4047"} = "HLA:HLA01748";	# B*40:47, B60Bw4, partial
$bw{"B-4701"} = "Bw4"; $bw{"B-4047"} = "Bw4";
$group{"B-4701"} = "B47"; $group{"B-4047"} = "B47";
$base{"B-4701"} = "B47"; $base{"B-4047"} = "B60";
push @subtype, ("B-4047");

my %b53;
$b53{"B-5301"} = "HLA00364";	#B*53:01:01:01
$bw{"B-5301"} = "Bw4";
$group{"B-5301"} = "B53";
$base{"B-5301"} = "B53";

my %b59;
$b59{"B-5901"} = "HLA00389";	#B*59:01:01:01
$bw{"B-5901"} = "Bw4";
$group{"B-5901"} = "B59";
$base{"B-5901"} = "B59";

my %b67;
$b67{"B-6701"} = "HLA00390";	# B*67:01:01
$b67{"B-6702"} = "HLA01374";	# B*67:02:01:01
$bw{"B-6701"} = "Bw6";
$bw{"B-6702"} = "Negative";
$group{"B-6701"} = "B67";
$group{"B-6702"} = "B67";
$base{"B-6701"} = "B67";
$base{"B-6702"} = "None";
push @subtype, ("B-6702");
my %b73;
$b73{"B-7301"} = "HLA00392";
$bw{"B-7301"} = "Negative";
$group{"B-7301"} = "B73";
$base{"B-7301"} = "B73";
my %b78;
$b78{"B-7801"} = "HLA16309";	# B*78:01:01:02 362 bp
$bw{"B-7801"} = "Bw6";
$group{"B-7801"} = "B78";
$base{"B-7801"} = "B78";
my %b81;
$b81{"B-8101"} = "HLA00398";
$bw{"B-8101"} = "Bw6";
$group{"B-8101"} = "B81";
$base{"B-8101"} = "B81";
my %b82b83;
$b82b83{"B-8201"} = "HLA00399";	#B*82:01:01:01
$bw{"B-8201"} = "Bw6";# $bw{"B-8301"} = "Bw6";
$group{"B-8201"} = "B82B83";# $group{"B-8301"} = "B82B83";
$base{"B-8201"} = "B82";# $base{"B-8301"} = "B83";


sub HLAB {
	my $gene = "B";
	return $gene;
}

sub HLAB_LEADER {
	my $leader = 23;		# B specific
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

sub BASETYPE {		# known serotype
	my %unique;
	my @basetype;
	foreach my $value (sort values %base ) {
		unless ( exists $unique{ $value } ) {
			push @basetype, $value;
			$unique{ $value } = 0;
		}
	}
	my $basetype_ref = \@basetype;
}

sub BROAD {
	my %broad;
	foreach my $base ( keys %base ) {
		if (( $base{ $base } eq "B51" ) || ( $base{ $base } eq "B5102" ) ||
			( $base{ $base } eq "B5103" ) || ( $base{ $base } eq "B52" )) {
			$broad{ $base } = "B5";
		}
		elsif ( $base{ $base } eq "B703" ) {
			$broad{ $base } = "B7";
		}
		elsif (( $base{ $base } eq "B62" ) || ( $base{ $base } eq "B63" ) ||
			( $base{ $base } eq "B75" ) || ( $base{ $base } eq "B76" ) ||
			( $base{ $base } eq "B77" )) {
			$broad{ $base } = "B15";
		}
		elsif (( $base{ $base } eq "B71" ) || ( $base{ $base } eq "B72" )) {
			$broad{ $base } = "B70";
		}
		elsif (( $base{ $base } eq "B64" ) || ( $base{ $base } eq "B65" )) {
			$broad{ $base } = "B14";
		}
		elsif (( $base{ $base } eq "B38" ) || ( $base{ $base } eq "B39" ) || ( $base{ $base } eq "B3902" )) {
			$broad{ $base } = "B16";
		}
		elsif (( $base{ $base } eq "B54" ) || ( $base{ $base } eq "B55" ) ||
			( $base { $base } eq "B56" )) {
			$broad{ $base } = "B22";
		}
		elsif (( $base{ $base } eq "B60" ) || ( $base{ $base } eq "B61" )) {
			$broad{ $base } = "B40";
		}
		elsif (( $base{ $base } eq "B44" ) || ( $base{ $base } eq "B45" )) {
			$broad{ $base } = "B12";
		}
		elsif (( $base{ $base } eq "B49" ) || ( $base{ $base } eq "B50" ) || ( $base{ $base } eq "B4005")) {
			$broad{ $base } = "B21";
		}
		elsif (( $base{ $base } eq "B57" ) || ( $base{ $base } eq "B58" )) {
			$broad{ $base } = "B17";
		}
		elsif ( $base{ $base } eq "B2708" ) {
			$broad{ $base } = "B27";
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
	push @combined, @b5; 
	push @combined, @b7; 
	push @combined, @b713; 
	push @combined, @b8; 
	push @combined, @b802; 
	push @combined, @b44; 	#B12
	push @combined, @b45; 	#B12
	push @combined, @b13;
	push @combined, @b14;
	push @combined, @b15;
	push @combined, @b63;
	push @combined, @b1523;
	push @combined, @b17;
	push @combined, @b38;	#B16
	push @combined, @b39;	#B16
	push @combined, @b18;
	push @combined, @b1809;
	push @combined, @b49;
	push @combined, @b50;
	push @combined, @b22;
	push @combined, @b27;
	push @combined, @b2708;
	push @combined, @b35;
	push @combined, @b37;
	push @combined, @b3705;
	push @combined, @b40;
	push @combined, @b48;
	push @combined, @b41;
	push @combined, @b42;
	push @combined, @b46;
	push @combined, @b47;
	push @combined, @b53;
	push @combined, @b59;
	push @combined, @b67;
	push @combined, @b73;
	push @combined, @b78;
	push @combined, @b81;
	push @combined, @b82b83;
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
	if ( $serotype eq "B5" ) {
		@residues = @b5; 
	}
	elsif ( $serotype eq "B7" ) {
		@residues = @b7; 
	}
	elsif ( $serotype eq "B713" ) {
		@residues = @b713; 
	}
	elsif ( $serotype eq "B8" ) {
		@residues = @b8;
	}
	elsif ( $serotype eq "B802" ) {
		@residues = @b802;
	}
	elsif ( $serotype eq "B44" ) {
		@residues = @b44;
	}
	elsif ( $serotype eq "B45" ) {
		@residues = @b45;
	}
	elsif ( $serotype eq "B13" ) {
		@residues = @b13;
	}
	elsif ( $serotype eq "B14" ) {
		@residues = @b14;
	}
	elsif ( $serotype eq "B15" ) {
		@residues = @b15;	#Bw6
	}
	elsif ( $serotype eq "B63" ) {
		@residues = @b63;	#Bw4
	}
	elsif ( $serotype eq "B1523" ) {
		@residues = @b1523;	#Bw4
	}
	elsif ( $serotype eq "B17" ) {
		@residues = @b17;	#Bw4
	}
	elsif ( $serotype eq "B38" ) {
		@residues = @b38;
	}
	elsif ( $serotype eq "B39" ) {
		@residues = @b39;
	}
	elsif ( $serotype eq "B18" ) {
		@residues = @b18;
	}
	elsif ( $serotype eq "B1809" ) {	#Bw4
		@residues = @b1809;
	}
	elsif ( $serotype eq "B49" ) {
		@residues = @b49;
	}
	elsif ( $serotype eq "B50" ) {
		@residues = @b50;
	}
	elsif ( $serotype eq "B22" ) {
		@residues = @b22;
	}
	elsif ( $serotype eq "B27" ) {
		@residues = @b27;
	}
	elsif ( $serotype eq "B2708" ) {
		@residues = @b2708;
	}
	elsif ( $serotype eq "B35" ) {
		@residues = @b35;
	}
	elsif ( $serotype eq "B37" ) {
		@residues = @b37;
	}
	elsif ( $serotype eq "B3705" ) {
		@residues = @b3705;
	}
	elsif ( $serotype eq "B40" ) {
		@residues = @b40;
	}
	elsif ( $serotype eq "B48" ) {
		@residues = @b48;
	}
	elsif ( $serotype eq "B41" ) {
		@residues = @b41;
	}
	elsif ( $serotype eq "B42" ) {
		@residues = @b42;
	}
	elsif ( $serotype eq "B46" ) {
		@residues = @b46;
	}
	elsif ( $serotype eq "B47" ) {
		@residues = @b47;
	}
	elsif ( $serotype eq "B53" ) {
		@residues = @b53;
	}
	elsif ( $serotype eq "B59" ) {
		@residues = @b59;
	}
	elsif ( $serotype eq "B67" ) {
		@residues = @b67;
	}
	elsif ( $serotype eq "B73" ) {
		@residues = @b73;
	}
	elsif ( $serotype eq "B78" ) {
		@residues = @b78;
	}
	elsif ( $serotype eq "B81" ) {
		@residues = @b81;
	}
	elsif ( $serotype eq "B82B83" ) {
		@residues = @b82b83;
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

	if ( $serotype eq "B5" ) {
		%ref = %b5; 
	}
	elsif ( $serotype eq "B7" ) {
		%ref = %b7; 
	}
	elsif ( $serotype eq "B713" ) {
		%ref = %b713; 
	}
	elsif ( $serotype eq "B8" ) {
		%ref = %b8;
	}
	elsif ( $serotype eq "B802" ) {
		%ref = %b802;
	}
	elsif ( $serotype eq "B44" ) {
		%ref = %b44;
	}
	elsif ( $serotype eq "B45" ) {
		%ref = %b45;
	}
	elsif ( $serotype eq "B13" ) {
		%ref = %b13;
	}
	elsif ( $serotype eq "B14" ) {
		%ref = %b14;
	}
	elsif ( $serotype eq "B15" ) {
		%ref = %b15;
	}
	elsif ( $serotype eq "B63" ) {
		%ref = %b63;
	}
	elsif ( $serotype eq "B1523" ) {
		%ref = %b1523;
	}
	elsif ( $serotype eq "B17" ) {
		%ref = %b17;
	}
	elsif ( $serotype eq "B38" ) {
		%ref = %b38;
	}
	elsif ( $serotype eq "B39" ) {
		%ref = %b39;
	}
	elsif ( $serotype eq "B18" ) {
		%ref = %b18;
	}
	elsif ( $serotype eq "B1809" ) {
		%ref = %b1809;
	}
	elsif ( $serotype eq "B49" ) {
		%ref = %b49;
	}
	elsif ( $serotype eq "B50" ) {
		%ref = %b50;
	}
	elsif ( $serotype eq "B22" ) {
		%ref = %b22;
	}
	elsif ( $serotype eq "B27" ) {
		%ref = %b27;
	}
	elsif ( $serotype eq "B2708" ) {
		%ref = %b2708;
	}
	elsif ( $serotype eq "B35" ) {
		%ref = %b35;
	}
	elsif ( $serotype eq "B37" ) {
		%ref = %b37;
	}
	elsif ( $serotype eq "B3705" ) {
		%ref = %b3705;
	}
	elsif ( $serotype eq "B40" ) {
		%ref = %b40;
	}
	elsif ( $serotype eq "B48" ) {
		%ref = %b48;
	}
	elsif ( $serotype eq "B41" ) {
		%ref = %b41;
	}
	elsif ( $serotype eq "B42" ) {
		%ref = %b42;
	}
	elsif ( $serotype eq "B46" ) {
		%ref = %b46;
	}
	elsif ( $serotype eq "B47" ) {
		%ref = %b47;
	}
	elsif ( $serotype eq "B53" ) {
		%ref = %b53;
	}
	elsif ( $serotype eq "B59" ) {
		%ref = %b59;
	}
	elsif ( $serotype eq "B67" ) {
		%ref = %b67;
	}
	elsif ( $serotype eq "B73" ) {
		%ref = %b73;
	}
	elsif ( $serotype eq "B78" ) {
		%ref = %b78;
	}
	elsif ( $serotype eq "B81" ) {
		%ref = %b81;
	}
	elsif ( $serotype eq "B82B83" ) {
		%ref = %b82b83;
	}
	else {		# "ALL", combine all
		%ref = (%b5,%b7,%b713,%b8,%b802,%b44,%b45,%b13,%b14,%b15,%b63,%b1523,%b17,%b38,%b39,%b18,%b1809,%b49,%b50,%b22,
		%b27,%b2708,%b35,%b37,%b3705,%b40,%b48,%b41,%b42,%b46,%b47,%b53,%b59,%b67,%b73,%b78,%b81,%b82b83);
	}
	
	return $ref_ref;
}

sub SERO {	# serotype and subtype information
	my @sero;
	my %ref = (%b5,%b7,%b713,%b8,%b802,%b44,%b45,%b13,%b14,%b15,%b63,%b1523,%b17,%b38,%b39,%b18,%b1809,%b49,%b50,%b22,
	%b27,%b2708,%b35,%b37,%b3705,%b40,%b48,%b41,%b42,%b46,%b47,%b53,%b59,%b67,%b73,%b78,%b81,%b82b83);
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
	my %tmp = (%b5,%b7,%b713,%b8,%b802,%b44,%b45,%b13,%b14,%b15,%b63,%b1523,%b17,%b38,%b39,%b18,%b1809,%b49,%b50,%b22,
	%b27,%b2708,%b35,%b37,%b3705,%b40,%b48,%b41,%b42,%b46,%b47,%b53,%b59,%b67,%b73,%b78,%b81,%b82b83);
	my %ref;
	my $key_ref = \%ref;
	for my $key ( sort keys %tmp ) {
		if ( $key =~ /B-510/ ) {
			$ref{$key} = "B\\*51";
		}				# B5
		elsif ( $key =~ /B-52/ ) {
			$ref{$key} = "B\\*52";
		}
		elsif ( $key =~ /B-07/ ) { 
			$ref{$key} = "B\\*07";
		}
		elsif ( $key =~ /B-08/ ) { 
			$ref{$key} = "B\\*08";
		}
		elsif ( $key =~ /B-44/ ) {
			$ref{$key} = "B\\*44";
		}
		elsif ( $key eq "B-4501" ) {
			$ref{$key} = "B\\*45";
		}
		elsif ( $key =~ /B-13/ ) {
			$ref{$key} = "B\\*13";
		}
		elsif ( $key =~ /B-14/ ) { 
			$ref{$key} = "B\\*14";
		}
		elsif ( $key eq "B-5701" ) {
			$ref{$key} = "B\\*57";
		}
		elsif ( $key eq "B-5801" ) {
			$ref{$key} = "B\\*58";
		}
		elsif ( $key =~ /B-15/ ) { 
			$ref{$key} = "B\\*15";
		}
		elsif ( $key =~ "B-38" ) {
			$ref{$key} = "B\\*38";
		}
		elsif ( $key =~ /B-39/ ) { 
			$ref{$key} = "B\\*39";
		}
		elsif ( $key =~ /B-18/ ) { 
			$ref{$key} = "B\\*18";
		}
		elsif ( $key eq "B-4901" ) {
			$ref{$key} = "B\\*49";
		}
		elsif ( $key eq "B-5001" ) {
			$ref{$key} = "B\\*50";
		}
		elsif ( $key eq "B-5401" ) {
			$ref{$key} = "B\\*54";
		}
		elsif ( $key eq "B-5501" ) {
			$ref{$key} = "B\\*55";
		}
		elsif ( $key =~ /B-56/ ) {
			$ref{$key} = "B\\*56";
		}
		elsif ( $key =~ /B-27/ ) { 
			$ref{$key} = "B\\*27";
		}
		elsif ( $key =~ /B-35/ ) {
			$ref{$key} = "B\\*35";
		}
		elsif ( $key =~ /B-37/ ) { 
			$ref{$key} = "B\\*37";
		}
		elsif ( $key =~ /B-40/ ) { 
			$ref{$key} = "B\\*40";
		}
		elsif ( $key =~ /B-41/ ) { 
			$ref{$key} = "B\\*41";
		}
		elsif ( $key =~ /B-48/ ) { 
			$ref{$key} = "B\\*48";
		}
		elsif ( $key =~ /B-42/ ) { 
			$ref{$key} = "B\\*42";
		}
		elsif ( $key =~ /B-46/ ) {
			$ref{$key} = "B\\*46";
		}
		elsif ( $key =~ /B-47/ ) {
			$ref{$key} = "B\\*47";
		}
		elsif ( $key =~ /B-53/ ) {
			$ref{$key} = "B\\*53";
		}
		elsif ( $key =~ /B-59/ ) {
			$ref{$key} = "B\\*59";
		}
		elsif ( $key =~ /B-67/ ) { 
			$ref{$key} = "B\\*67";
		}
		elsif ( $key =~ /B-73/ ) {
			$ref{$key} = "B\\*73";
		}
		elsif ( $key =~ /B-78/ ) {
			$ref{$key} = "B\\*78";
		}
		elsif ( $key =~ /B-81/ ) {
			$ref{$key} = "B\\*81";
		}
		elsif ( $key =~ /B-82/ ) {
			$ref{$key} = "B\\*82";
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
	my $seq = "N" x 25;

	$partial{ "B-3803" } = $seq;	# still partial
	$partial{ "B-4047" } = $seq;	# still partial
	$partial{ "B-4021" } = $seq;	# still partial

		
	return $partial_ref;
}

# modified on October 18 2022
sub WHO {
	my %who;
	my $whotype_ref = \%who;
	$who{"B-5101"} = "B51"; $who{"B-5102"} = "B5102"; $who{"B-5103"} = "B5103"; $who{"B-5201"} = "B52"; $who{"B-0702"} = "B7"; $who{"B-0703"} = "B703"; #B5, B7
	$who{"B-0801"} = "B8"; $who{"B-4402"} = "B44"; $who{"B-4501"} = "B45"; $who{"B-1302"} = "B13"; $who{"B-1401"} = "B64"; $who{"B-1402"} = "B65"; # B8, B12, B13, B14
	$who{"B-1501"} = "B62"; $who{"B-1516"} = "B63"; $who{"B-1517"} = "B63"; $who{"B-1502"} = "B75"; $who{"B-1512"} = "B76"; $who{"B-1514"} = "B76"; $who{"B-1513"} = "B77"; #B15
	$who{"B-3801"} = "B38"; $who{"B-3901"} = "B3901"; $who{"B-3902"} = "B3902"; $who{"B-5701"} = "B57"; $who{"B-5801"} = "B58"; $who{"B-1801"} = "B18"; # B16, B17, B18
	$who{"B-4005"} = "B4005"; $who{"B-4901"} = "B49"; $who{"B-5001"} = "B50"; $who{"B-5401"} = "B54"; $who{"B-5501"} = "B55";$who{"B-5601"} = "B56"; # B21, B22
	$who{"B-2705"} = "B27"; $who{"B-2708"} = "B2708"; $who{"B-3501"} = "B35"; $who{"B-3701"} = "B37"; $who{"B-4001"} = "B60";$who{"B-4002"} = "B61"; # B27, B35, B37, B40 
	$who{"B-4101"} = "B41"; $who{"B-4201"} = "B42"; $who{"B-4601"} = "B46"; $who{"B-4701"} = "B47"; $who{"B-4801"} = "B48"; $who{"B-5301"} = "B53"; # B41, B42, B46, B47, B48, B53
	$who{"B-5901"} = "B59"; $who{"B-6701"} = "B67"; $who{"B-1510"} = "B71"; $who{"B-1503"} = "B72"; $who{"B-7301"} = "B73"; $who{"B-7801"} = "B78"; #B59, B67, B70, B73, B78
	$who{"B-8101"} = "B81"; $who{"B-8201"} = "B82";# $who{"B-8301"} = "B83";
	return $whotype_ref;
}


1;
