#!/usr/bin/perl -w
# Openfile.pm

# Author: Kazutoyo Osoegawa, Ph.D.
# email: kazutoyo@stanford.edu
# phone: 650-724-0169

# line 18: modified from "\n" to "\r\n"
# Last modified on February 2 2022

use strict;

package Openfile;

sub open_file {
	print "Enter a file name to open:\n";
	my $enter = <STDIN>;
	open(FILENAME, $enter) || die "Error: $!";
	# The open() function opens a file indicated and attaches it to a filehandle called FILENAME.
	
	my $line;
	my @list = ();		
	for ($line = 0; <FILENAME>; $line++) {
	# Calling reading operator (<>) returns the next line from the given filehandle.
		if ($_ eq "\r\n") {		## To remove new line characters for PC
			next;
		}
		
	chomp ($_);
	$list[$line] = $_;
	}
	close FILENAME;
	return @list;
}


# Open a file
# Enter each line in an array format
sub open_file_from_list {

	my ($enter) = @_;

	open(FILENAME, $enter) || die "Error: $!";
	# The open() function opens a file indicated and attaches it to a filehandle called FILENAME.
	
	my @list = ();		
	for (my $line = 0; <FILENAME>; $line++) {
	# Calling reading operator (<>) returns the next line from the given filehandle.
		if ($_ eq "\r\n") {		## To remove new line characters for PC
			next;
		}
		
	chomp ($_);
	$list[$line] = $_;
	}
	close FILENAME;
	return @list;
}


1;
# It must return the value of 1 at the end of the Perl module script
# this is to indicate to the script using the module that it loaded successfully.
