#!/usr/bin/perl
use strict;
no strict 'refs';
use warnings;
use Cwd qw(cwd);

# This script changes the letters to numbers for a given keyfile. The format of the keyfile should be
# PLATE <tab> ROW <tab> COLUMN <tab> ORF

my @NEWKEYFILE;
my @info;
my $file = "keyfile.txt";								# Open INPUT file "keyfile.txt" or report an error				
open (INPUT, $file ) or die "$file cannot be opened. $!";																															
my @KEYFILE= <INPUT>;									# put the data into an array @KEYFILE
close INPUT or die "Cannot close the file: $!";			# close the file

my $keysize = scalar@KEYFILE;
if ($keysize == 1) {
	die "Your keyfile only has one entry, this is likely because keyfile.txt does not have standard line endings, therefore is recognised as one big line of data! Please save it with different line endings!\n";
}

foreach (@KEYFILE) {									# iterate over the KEYFILE
	@info=split/\t/;									# split the data by tabs
	
	chomp($info[0]);									# remove any trailing invisible characters (e.g. tabs or newlines) from the data
	chomp($info[1]);
	chomp($info[2]);
	chomp($info[3]);
	my $orf_pattern='^[Y|y][A-Pa-p][L|R|l|r][0-9]{3}[W|w|C|c](?:-[A-Za-z])?$';# define what a normal ORF would look like
	my $orf = ($info[3]);									# get the orf from the data
	if (($orf !~ /$orf_pattern/)&&($orf !~ /BLANK/)&&($orf !~ /CONTROL/)&&($orf !~ /ORF/)) { #check it looks OK, otherwise....
		print "Warning: $orf doesn't look like a real ORF name in your keyfile.\n";# print an error
	}							
	$info[1] = letterToNumber($info[1]);				# change row letter to a number
	if (($info[1] !~ /^\d+$/)&&($info[1] !~ /Row/)) {	# check that the row is really a number, otherwise...
		print "Warning: A row number $info[1] is in your keyfile!\n";# print a warning
	}
	if (($info[2] !~ /^\d+$/)&&($info[2] !~ /Column/)) {# check that the column is a number, otherwise...
		print "Warning: A column number $info[2] is in your keyfile!\n";# print a warning
	}
	if (($info[0] !~ /^\d+$/)&&($info[2] !~ /Plate/)) {
		print "Warning: A plate number $info[0] is in your keyfile!\n"; # print a warning
	}
	
	my $entry = "$info[0]\t$info[1]\t$info[2]\t$info[3]\n";	# create a modified entry for the keyfile
	push (@NEWKEYFILE, $entry);								# add the entry to a newkeyfile array
}

open(my $fh, ">newkeyfile.txt") or die "results file cannot be opened. $!";	# write output file called "results.tab"
print $fh "@NEWKEYFILE";

################################## NUMBERTOLETTER ################################
sub letterToNumber {							# Changes row letters to numbers

$_ = $info[1];

	if ( $_ eq "A"){ $_ = 1; }				# there must be a better way to do this!
	elsif( $_ eq "B" ) { $_ = 2; }
	elsif( $_ eq "C" ) { $_ = 3; }
	elsif( $_ eq "D" ) { $_ = 4; }
	elsif( $_ eq "E" ) { $_ = 5; }
	elsif( $_ eq "F" ) { $_ = 6; }
	elsif( $_ eq "G" ) { $_ = 7; }
	elsif( $_ eq "H" ) { $_ = 8; }
	elsif( $_ eq "I" ) { $_ = 9; }
	elsif( $_ eq "J" ) { $_ = 10; }
	elsif( $_ eq "K" ) { $_ = 11; }
	elsif( $_ eq "L" ) { $_ = 12; }
	elsif( $_ eq "M" ) { $_= 13; }
	elsif( $_ eq "N" ) { $_ = 14; }
	elsif( $_ eq "O" ) { $_ = 15; }
	elsif( $_ eq "P" ) { $_ = 16; }
	
	$info[1] = $_;
}