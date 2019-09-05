#!/usr/bin/perl
use strict;
no strict 'refs';
use warnings;
use Cwd qw(cwd);

### This script uses two files; 1. a keyfile and 2. a "colonyAreas.txt" file (the latter usually from cm-engine)
# The keyfile must be called "keyfile.txt" and have this format....
# Plate <tab> Row <tab> Column <tab> ORF 
# The keyfile can have more data but the first four columns of data must be as listed above and in the same order.

# The output from cm_engine should be called "colonyAreas.txt" and be in the usual format.
# It is important that the image names used are in the standard format, which is ...
# plasmidName,PlateNumber,etc

#### ******* SET THE NUMBER OF COLONIES & STRAINS ON YOUR PLATES *******
# The numbers here should be 96, 384 or 1536
my $platesize = 1536  ;						# how many colonies were on each plate

my $density = 384;							# how many STRAINS on each plate (this should be the same as the keyfile)


# **************************************************************************

our $repeats=$platesize/$density;						# repeats is the number of replicates of each strain in your data, typically 4 or 16
our ($maxcol, $maxrow);
if (($repeats == 1)||($repeats == 4)||($repeats == 16)) {} else {	# check if the number of repeats makes sense, if not report an error
	die "You specified a keyfile density of $density strains per plate and a platesize of $platesize colonies per plate, this gives $repeats repeats of each strain (this number should be 1, 4 or 16)!\n"
}

if ($platesize == 1536) {								# if the data are in 1536 format...
	$maxcol = 47;										# 47 should be the penultimate column
	$maxrow = 31;										# 31 should be the penultimate row
} elsif ($platesize == 384) {							# if the data are in 384 format...
	$maxcol = 23;										# 23 should be the penultimate column
	$maxrow = 15;										# 15 should be the penultimate row
} elsif ($platesize == 96) {							# if the data are in 96 format...
	$maxcol = 11;										# 11 should be the penultimate column
	$maxrow = 7;										# 7 should be the penultimate row
} else {												# otherwise print a warning
	print "WARNING: Have you set the platesize and density density correctly? Your settings are currently $platesize colonies per plate, with $density strains per plate. Meaning that you have $repeats repeats on each plate.\n ";
}

#chdir("/Users/thorpep/Desktop/Peter/Programs/PERL/colArea_extractor/newColanalyser");
my $file = "keyfile.txt";								# Open INPUT file "keyfile.txt" or report an error				
open (INPUT, $file ) or die "$file cannot be opened. $!";																															
my @KEYFILE= <INPUT>;									# put the data into an array @KEYFILE
close INPUT or die "Cannot close the file: $!";			# close the file

my $keysize = scalar@KEYFILE;
if ($keysize == 1) {
	die "Your keyfile only has one entry, this is likely because keyfile.txt does not have standard line endings, therefore is recognised as one big line of data! Please save it with different line endings!\n";
}

our %keyHash;											# define variables for the script
our @info;
our @data;
our @RESULTS;
our ($i,$j,$row,$col,$orf,$platename,$plasmidname,$platenum,$rownum,$colnum,$keylook);


# This section of code uses your keyfile to create a new keyfile (in the form of a hash) that is in the same spatial organisation as 
# your screen, so potentially it could take a 96 strains per plate keyfile and create a new keyfile in 1536 format.
# The new keyfile has the 'keys' defined as plateNumber_rowNumber_columnNumber, e.g. 1_1_4 (plate 1, row 1, column 4)
# The values of the keyfile hash are just the ORF names

foreach (@KEYFILE) {									# iterate over the KEYFILE
	@info=split/\t/;									# split the data by tabs
	
	chomp($info[0]);									# remove any trailing invisible characters (e.g. tabs or newlines) from the data
	chomp($info[1]);
	chomp($info[2]);
	chomp($info[3]);
	my $orf_pattern='^[Y|y][A-Pa-p][L|R|l|r][0-9]{3}[W|w|C|c](?:-[A-Za-z])?$';# define what a normal ORF would look like
	$orf = ($info[3]);	
	chomp($orf);
	trim ($orf);								# get the orf from the data
	if (($orf !~ /$orf_pattern/)&&($orf !~ /BLANK/)&&($orf !~ /ORF/)) { #check it looks OK, otherwise....
		print "Warning: $orf doesn't look like a real ORF name in your keyfile.\n";# print an error
	}							
	$info[1] = letterToNumber($info[1]);				# change row letter to a number
	if (($info[1] !~ /^\d+$/)&&($info[1] !~ /Row/)) {	# check that the row is really a number, otherwise...
		print "Warning: A row number $info[1] is in your keyfile!\n";# print a warning
	}
	if (($info[2] !~ /^\d+$/)&&($info[2] !~ /Column/)) {# check that the column is a number, otherwise...
		print "Warning: A column number $info[2] is in your keyfile!\n";# print a warning
	}
	
	if (($repeats == 1)&&($info[0] !~ m/Plate #/)) {	# if there is one repeat of each strain on the plates e.g. your keyfile and data are both 385 format
		my $entry = "$info[0]_$info[1]_$info[2]";		# then do not change the row and column numbers, make an entry with the plate_row_column
		$keyHash{ $entry } = $orf;						# put the data into a hash with keys being a number from 1 to the end of the plate.
	} elsif (($repeats == 4)&&($info[0] !~ m/Plate #/))	{ # otherwise, if the data are repeated four times...
		foreach  $i (0..1){								# iterate twice...
			foreach  $j (0..1) {						# iterate twice...
				$row = ($info[1]*2)-$i;					# new row number is 2x old row, minus either 0 or 1
				$col = ($info[2]*2)-$j;					# new column number is 2x old column, minus either 0 or 1
				my $entry = "$info[0]_".$row."_"."$col";# make an entry with the plate_row_column
				$keyHash{ $entry } = $orf;				# put the data into a hash with keys being a number from 1 to the end of the plate.
			}
		}
	} elsif (($repeats == 16)&&($info[0] !~ m/Plate #/))	{# otherwise, if the data are repeated four times...
		foreach  $i (0..3){								# iterate four times...
			foreach  $j (0..3) {						# iterate four times...
				$row = ($info[1]*4)-$i;					# new row number is 4x old row, minus 0, 1, 2 or 3 
				$col = ($info[2]*4)-$j;					# new column number is 4x old row, minus 0, 1, 2 or 3 
				my $entry = "$info[0]_".$row."_"."$col";# make an entry with the plate_row_column
				$keyHash{ $entry } = $orf;				# put the data into a hash with keys being a number from 1 to the end of the plate.
			}
		}	
	}									
}

# This section of code goes through your colonyAreas data and looks up each ORF to create a new file with the ORFs listed.

$file = "colonyAreas.txt";								# Open INPUT file "keyfile.txt" or report an error				
open (INPUT, $file ) or die "$file cannot be opened. $!";																															
my @COLONYSIZES= <INPUT>;									# put the data into an array @KEYFILE
close INPUT or die "Cannot close the file: $!";			# close the file

$keysize = scalar@COLONYSIZES;
if ($keysize == 1) {
	die "Your dataset only has one entry, this is likely because colonyAreas.txt does not have standard line endings, therefore is recognised as one big line of data! Please save it with different line endings!\n";
}

print "\nThis is a list of each plasmid and plate in your dataset\n";

foreach (@COLONYSIZES) {								# Iterate over each line of the ColonyAreas file
	@info=split/\t/;									# split the data by tabs
	chomp($info[0]);									# remove any trailing invisible characters (e.g. tabs or newlines) from the data
	my $text='^[A-Za-z]';								# define text
	if ($info[0] =~ /$text/ ){							# check if this is the plate definition, if so...
		@data=split/,/;									# split the data by commas
		chomp($data[0]);								# remove any trailing invisible characters (e.g. tabs or newlines) from the data
		$data[1] =~ s/.tif//;							# if the second data line includes the file suffix, stip it out.
		chomp($data[1]);								# remove any trailing invisible characters (e.g. tabs or newlines) from the data
		trim ($data[1]);								# trim any whitespace around the platename
		if (defined $data[2]) {							# if there is more information in your imagenames
			$data[2] =~ s/.tif//;						# if it includes the file suffix, stip it out.
			chomp($data[2]);							# remove any trailing invisible characters (e.g. tabs or newlines) from the data
			trim ($data[2]);							# trim any whitespace around the extra_data
			if ($data[2] =~ /^[\w\W\d]/) {				#  if the extra data looks like number or text
				#print "I found some extra information! It is >$data[2]<\n";
				$data[0] = $data[0] . "_" . $data[2];	# add the extra info onto the plasmid name
			}
		}
		$plasmidname = $data[0];						# define the plasmid name as the initial data
		if ($data[1] !~ /^\d+$/) {						# if the data is not numeric...
			print "****************\nWarning: One of your plate numbers is $data[1], this is not numeric and so was not included in the output!\n*****************\n";	#print an error
			$platenum = "error";						# define platenum as an error
		} else {
			$platenum = $data[1];						# define the plate number as the second item 
		}					
		print "Plate $platenum , plasmid $plasmidname.\n"; # print out the plasmid and plate
		$rownum = 1;									# reset the row to 1
		$colnum = 1;									# reset the column to 1
	} elsif (($info[0]<0)||($info[0]>1000)) {			# Otherwise, check if the colony sizes look normal
		print "One of the colony sizes is $info[0]! This has been excluded.";	# report an error if the colony sizes are out of the range
	} else {											# Otherwise,...
		$keylook = $platenum . "_" . $rownum . "_" . $colnum; # look in the keyfile hash, for a key that matches plate_row_col (e.g. 1_2_5)
		if ((exists $keyHash{$keylook})&&($platenum ne "error")) {				# check if the current position is in the keyfile data hash, if so...
			$orf = $keyHash{$keylook};					# Find the ORF name that corresponds to this position
			chomp ($orf);
			my $entry = "$plasmidname\t$platenum\t$rownum\t$colnum\t$orf\t$info[0]\n"; # create an entry which includes all the data including colony size
			push (@RESULTS, $entry);					#put this entry into the results
		} elsif ($platenum ne "error") {										# otherwise print a warning
			print "Warning: No ORF was found in your keyfile data for Plate $platenum, Row $rownum, Column $colnum.\n"
		}
		if ($rownum > $maxrow) {						# if we have reached the last row in the data...
			$rownum = 1;								# reset the column number to 1
			$colnum = $colnum + 1;						# and add one to the column number
		} else {										# otherwise...
			$rownum =$rownum + 1;						# add one to the row number
		}
	}
}

my $entry = "Plasmid\tPlate\tRow\tColumn\tORF\tColony size\n";				# create a row of headings
unshift (@RESULTS, $entry);													# put this at the top of the results array
open(my $fh, ">results.tab") or die "results file cannot be opened. $!";	# write output file called "results.tab"
print $fh "@RESULTS";


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
########################      trim      ##############################
#subroutine to trim off the white space from both ends of each string
sub trim {
	$_[0]=~ s/^\s+//;
	$_[0]=~ s/\s+$//;
	return;
}

