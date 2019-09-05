#!/usr/bin/perl
use strict;
no strict 'refs';
use warnings;

############################## Z score smoothing script  #####################################
# This script takes an input table with the following format, usually copied from a spreadsheet.
# Please ensure that you have saved the input table using a text editor that uses standard line breaks
# Note that saving the file from excel often does not have standard line breaks.
#
# plate <tab> row <tab> column <tab> ORF <tab> score <tab> score etc... as many scores as you like
# please include column headings in the first row, particularly if you have many scores.
# The first column heading must be the word "plate" (case not important), otherwise the heading won't be recognised
#
# There can be any number of plates, and each plate can have any number of rows and columns
# Columns are numerical, where as rows are letters (A, B, C etc)
#
# The script adjusts the scores (assumed to be Z scores) based upon column and row averages relative to the median plate value
# The scores should be numerical not text, although the script makes some attempt to parse out text from the scores
# 
#
# Written in 2013 by Peter Thorpe

# ******** define variables here ***********#
my @sorted=();								#
my $row;									#
my $col;									#
#*******************************************#

opendata();									# run the opendata subroutine (puts the data into an array and works out how many plates there are)

for (my $s=1; $s<(1+$main::size); $s++) {	# iterate over each screen in the input file
	
	%main::newscore = ();					# clear the newscore hash
	%main::oldscore = ();					# clear the oldscore hash

	my $num=$s-1;							
	$main::xname = $main::NAMES[$s-1];		# get the name of the current screen from the NAMES array
	#print "Screen $s has the heading $main::xname\n";
			
	for ($main::n=1; $main::n<($main::max+1); $main::n++) {	# iterate over each plate (plate number is defined by $n), up to the max number of plates
		
		platemedian($main::medcheck="firstpass",%main::HoH);# run the platemedian subroutine (get the median score for plate $n)
		smoothcol();						# run the smoothcol subroutine (adjust the z scores based on column averages)
		platemedian($main::medcheck="secondpass",%main::newscore);# run the platemedian subroutine (get the median score for plate $n)
		smoothrow(%main::newscore);			# run the smoothrow subroutine (adjust the z scores based on row averages)
		levelplates();						# run the levelplates subroutine
	}
	collateresults(@main::RESULTS);			# run the collateresults subroutine
}

sortlabelsave(@main::RESULTS);				# run the sortlabelsave subroutine

####################################### opendata ##########################################
# puts the data from results.txt into an array called @ARRAY, lines without a numeric z score are omitted
# also calculates the maximum number of plates +1 = $max
sub opendata {
	
	@main::NAMES=();
	%main::HoH = ();										# clear the HoH hash
	
	my $file = "results.txt";								# Open INPUT file "results.txt" or report an error	
	open(INPUT, $file) or die "results.txt file cannot be opened. $!";																															
	my @TOTARRAY= <INPUT>;									# put the data into an array @TOTARRAY
	close INPUT or die "Cannot close the file: $!";			# close the file
	
	my $alen=scalar(@TOTARRAY);								# get the length of the input array
	if ($alen < 2) {										# check that there are at least 2 lines of data in the input file or terminate with an error message
		die "Your input file was not in the correct format, only $alen of data were found. Did you save the 'results.txt' file in the correct format using a text editor?";
	}
	
	my @firstline = split(/\t/,$TOTARRAY[0]);				# get the first line of the array
	$main::size=scalar(@firstline);							# get the size of each line
	$main::size = $main::size - 4;							# subtract the four columns that are not scores
	
	$main::max = 0;											# set $max to zero
	$main::maxr = 0;										# set $maxr to zero
	$main::maxc = 0;										# set $maxc to zero
	
	foreach (@TOTARRAY) {
		my @info=split/\t/;									# splite the data by tabs
		foreach (@info) {									# iterate over the data on this line and
			chomp($_);										# and remove newline characters
			trim($_);										# and remove any flanking whitespaces
		}
		if ($info[0] =~ m/Plate/i) {						# if this is the headings row
			for (my $q=0; $q<$main::size; $q++) {			# for each screen...
				push(@main::NAMES, "$info[$q+4]");			# put the heading name into the @NAMES array
			}
		} elsif (($info[0] !~ m/Plate/i)&&(defined ($info[3])&&($info[3] ne ""))) {	# if this is not the headings row and has an orf
			$info[0] =~ s/\D//g;							# remove any non-digits from the plate number
			$info[2] =~ s/\D//g;							# remove any non-digits from the column number
			# check if the plate, row, column and orf are all in the correct format, if not report an error
			if (($info[0] !~ /[\d|\d\d]/)||($info[1] !~ /[\d|\d\d]/)||($info[2] !~ /[\d|\d\d]/)||($info[3] !~ /[Y|y][A-Pa-p][R|L|r|l][0-9]{3}[W|w|C|c](?:-[A-Za-z])?/)) {
				print "There is a format error with ORF $info[3] on plate $info[0], row $info[1], column $info[2], it has not been included in the analysis.\n";
			} else {
				my $entry="$info[0]\t$info[1]\t$info[2]\t$info[3]"; # define the entry based on plate, row, column and orf
				push (@main::RESULTS, "$entry\n");				# put the identifiers into a RESULTS array
				if ($info[0] > $main::max) {					# check if this is the largest plate number
					$main::max = $info[0];						# if so, set this as the maximum plate number
				}
			}
		}
	}
	
	my $limit=$main::size-1;								# limit is 1 less than the number of screens
	
	for my $index(0..$limit) {								# iterate based over the number of screens, from zero to limit.
		my $pos=$index+4;									# set the position of the score for this screen
		
		foreach (@TOTARRAY) {
			my @info=split/\t/;								# splite the data by tabs
			foreach (@info) {								# iterate over the data on this line and
				chomp($_);									# remove newline characters
				trim($_);									# and remove spaces and tabs 
			}
		 	if (($info[0] !~ m/Plate/i)&&(defined $info[$pos])) {	# if this is not the headings row
								
				$info[0] =~ s/\D//g;						# remove any non-digits from the plate number
				$info[2] =~ s/\D//g;						# remove any non-digits from the column number
							
				$main::row = $info[1];						# define the row number
				letterToNumber($main::row);					# Changes row letters to numbers using the letterToNumber subroutine
				$info[1] = $main::row;						# reinsert the row, now as a number into the @info array
				
				my $codec="$info[0],$info[1],$info[2]";		# make a codec identifier for each plate,row,
				
				if ($info[0] > $main::max) {				# check if this is the largest plate number
					$main::max = $info[0];					# if so, set this as the maximum plate number
				}
				if ($info[1] > $main::maxr) {				# check if this is the largest row number
					$main::maxr = $info[1];					# if so, set this as the maximum row number
				}
				if ($info[2] > $main::maxc) {				# check if this is the largest column number
					$main::maxc = $info[2];					# if so, set this as the maximum column number
				}
				# check if the data exists and is numeric...
				if ($info[$pos] =~ m/^-?(0|([1-9][0-9]*))(\.[0-9]*)?([eE][-+]?[0-9]*)?$/) {	
					# regex matches any number of intergers followed by . followed by anynumber of intergers followed by E or e followed by anynumber of integers
					my $score = ($info[$pos]);				# get the score from the appropriate position in the data
					my $sname = $main::NAMES[$index];		# get the screen name from the NAMES array
					my $codec = "$info[0],$info[1],$info[2]";# define a universal identifier plate,row,column			
					$main::HoH {$sname}->{$codec} = $score;	# put the screen name, orf and score into a hash of hashes 
				}
			} 
		}
	}
	return(@main::RESULTS,@main::NAMES, $main::size, $main::max, %main::HoH);																				
}

#################################### platemedian ##########################################
# This subroutine takes 
# the outputs from this routine are two arrays, @plate is just a subset of the array for this ($n) specific plate
# @z is simply a list of all the scores from this plate
sub platemedian {
	my @plate=();											# clear the @plate array
		
	if ($main::medcheck eq "firstpass") {	
		for (my $col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...
			for (my $row = 1; $row<($main::maxr+1); $row++) { # for each row 1 through max..
				my $codec = "$main::n,$row,$col";			# define a plate code for 
				if (exists $main::HoH{$main::xname}{$codec}) {
					my $score = $main::HoH{$main::xname}{$codec};# if so, define the score
					chomp($score);
					push (@plate,$score)
				}
			}
		}
	}
	
	if ($main::medcheck eq "secondpass") {	
		for ($col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...
			for (my $row = 1; $row<($main::maxr+1); $row++) { 	# for each row 1 through max..
				my $codec = "$main::n,$row,$col";			# define the codec for each position on the plate
				if (exists $main::newscore{$codec}) {		# check if that position exists in the newscore hash
					my $score = $main::newscore{$codec};	# if so, define the score
					chomp($score);							# remove any newline characters
					push (@plate,$score)					# add the score to the plate array
				}
			}
		}
	}
	median(@plate);											# run the median subroutine
}

######################################### median ##########################################
# extracts the median number from the input array
sub median {
	@sorted = sort { $a <=> $b } @_;						# sort the input array by number
	my $len = scalar(@sorted);								# get the size of the sorted list ($len)
	my $median=0;											# set the median value to zero
	if ($len > 0) {											# if there are values in the array
		if (@sorted % 2) {									# if this length is an odd number								
	    	$median = $sorted[ int($len/2) ];				# the median value is the mid value in the list
		} else {											# otheriwise...
	    	$median = ($sorted[int($len/2)-1] + $sorted[int($len/2)])/2;	# the median is the average of the two middle values
		}
	}
	$main::median=$median;									# put the resulting median into 
	return ($main::median);
}

###################################### smoothcol ##########################################
# This subrountine to adjust scores based upon column averages. The mean of each column is 
# compared to the plate median. The difference between these two is used to correct the 
# individual scores in the column. Hence, in the resulting data the mean of the data in each
# column should equal the plate median.
sub smoothcol {

	for ($col=1; $col<($main::maxc+1); $col++) {			# for each column 1 through max...
		@main::c=();										# clear the @c array
		for ($row = 1; $row<($main::maxr+1); $row++) { 		# for each row
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::HoH{$main::xname}{$codec}) {	# check if this identifier is in the hash of hases
				my $score = $main::HoH{$main::xname}{$codec};# if so, define the score
				chomp ($score);								# remove any newline characters
				if (($score==$score*1)&&($score ne "")) {	# check if the score is numeric
						push(@main::c,  $score);			# put the score into the @c array
				}
			}
		}
				
		average(@main::c);									# run the average subroutine (get the average scores in each column)
		
				# check that there are sufficient values to average and that the correction is not a crazy value
		if (($main::count<4)||($main::correction>4)||($main::correction<-4)) {						# otherwise...
				#print "Error with plate $main::n, column $main::c, correction is $main::correction, ";	# report an error...
			if ($main::count<4) {
				#print "there are only $main::count samples to average, so no column correction was applied.\n"; # report an error...
			} else {
				#	print "the correction value is $main::correction, which is outside the acceptable range - no column correction was applied.\n";
			}
			$main::correction=0;							# set $correction to zero, since there was an error
		}
		
		for ($row = 1; $row<($main::maxr+1); $row++) { 		# for each row
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::HoH{$main::xname}{$codec}) {	# check if this identifier is in the hash of hases
				my $score = $main::HoH{$main::xname}{$codec};# if so, define the score
				chomp ($score);								# remove any newline characters
				if (($score==$score*1)&&($score ne "")) {	# check if the score is numeric
						my $nscore = $score - $main::correction;# calculate the newscore
						$main::oldscore{$codec} = $score;	# put the oldscore in the oldscore hash
						$main::newscore{$codec} = $nscore;	# put the newscore in the newscore hash
				}
			}
		}
	}
}

####################################### smoothrow ##########################################
# This subrountine to adjust scores based upon column averages. The mean of each column is 
# compared to the plate median. The difference between these two is used to correct the 
# individual scores in the column. Hence, in the resulting data the mean of the data in each
# column should equal the plate median.
sub smoothrow {
	
	for ($row=1; $row<($main::maxr+1); $row++) {			# for each row 1 through max...
		@main::r=();										# clear the @r array
		
		for ($col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...	
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::newscore{$codec}) {			# check if this identifier is in the newscore hash
				my $score = $main::newscore{$codec};				# get the score from the newscore hash
				chomp ($score);								# remove any newline characters
				push(@main::r,$score);						# put the score into the @c array
			}
		}
					
		average(@main::r);									# run the average subroutine (get the average scores in each row)
		
		# check that there are sufficient values to average and that the correction is not a crazy value
		if (($main::count<5)||($main::correction>4)||($main::correction<-4)) {# otherwise...
				#print "Error with plate $main::n, row $row, correction is $main::correction, ";	# report an error
			if ($main::count<5) {
				#print "there are only $main::count samples to average, so no row correction was applied.\n";
			} else {
				print "the correction value is $main::correction, which is outside the acceptable range - no row correction was applied.\n";
			}
			$main::correction=0;							# set $correction to zero, since there was an error
		}
		$col=1;												# reset $col to 1
		for ($col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...	
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::newscore{$codec}) {			# check if this identifier is in the newscore hash
				my $score = $main::newscore{$codec};		# get the score from the newscore hash
				chomp ($score);								# remove any newline characters
				my $nscore = $score - $main::correction;	# calculate the newscore
				delete $main::newscore{$codec};				# delete the column adjusted score
				$main::newscore{$codec} = $nscore;			# replace the cell with a new cell with an extra plasmid
			}
		}
	}
}

######################################## average ##########################################
# This subroutine calculates the median value from an input list (@_)
# This value is confusingly termed "mean" for historical reasons - it is not the mean! 
# This median is then compared with the median plate value (the median score of the plate)
# and returns a correction value, which is the difference between the median of the plate
# and the medain of the current input. 
sub average {
	
	@sorted = sort { $a <=> $b } @_;						# sort the input array by number
	$main::count = scalar(@sorted);							# get the size of the array
	my $mean=0;											# set the median value to zero
	if ($main::count > 0) {											# if there are values in the array
		if (@sorted % 2) {									# if this length is an odd number								
	    	$mean = $sorted[ int($main::count/2) ];		# the median value is the mid value in the list
		} else {											# otheriwise...
	    	$mean = ($sorted[int($main::count/2)-1] + $sorted[int($main::count/2)])/2;	# the median is the average of the two middle values
		}
		$main::correction = $mean-$main::median;			# calculate the correction value relative to the median plate value
	} else {												# if there were no values...
		$main::correction = 0;								# set the ratio to zero
	}
	return ($main::count, $main::correction);	
}

#################################### levelplates ##########################################
# This subroutine determines the median scores for a given plate and
# adjusts all the scores on the plate such that the median is zero
sub levelplates {
	
	@main::levels=();										# clear the levels array
	for (my $row=1; $row<($main::maxr+1); $row++) {			# for each row 1 through max...
		for ($col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...	
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::newscore{$codec}) {			# check if this identifier is in the newscore hash
				my $score = $main::newscore{$codec};		# get the score from the newscore hash
				chomp ($score);								# remove any newline characters
				push(@main::levels,$score);					# put the score into the @levels array
			}
		}
	}
	
	median(@main::levels);									# run the median subroutine
	
	for ($row=1; $row<($main::maxr+1); $row++) {			# for each row 1 through max...
		for ($col=1; $col<($main::maxc+1); $col++) {		# for each column 1 through max...	
			my $codec = "$main::n,$row,$col";				# make an identifier
			if (exists $main::newscore{$codec}) {			# check if this identifier is in the newscore hash
				my $score = $main::newscore{$codec};		# get the score from the newscore hash
				chomp ($score);								# remove any newline characters
				$score = $score - $main::median;			# correct the scores based upon the plate median values
				delete $main::newscore{$codec};				# delete the column adjusted score
				$main::newscore{$codec} = $score;			# replace the cell with a new cell with an extra plasmid
			}
		}
	}
}

################################## collateresults #########################################
# add the latest results into the results array
sub collateresults {
	my @TEMPRESULTS=();										# blank the tempresults array
	foreach(@main::RESULTS) {								# iterate over the results array
		my @info=split/\t/;									# split the data by tabs, into an array @info
		my $isize=scalar(@info);							# get the size of the line
		my $entry="";										# blank the $entry variable
		foreach(@info) {									# iterate over the data on this line
			chomp($_);										# remove any newline characters
			trim($_);
			$entry=$entry . $_ . "\t";						# add each piece of data to the entry
		}
		
		$main::row = $info[1];								# define the row number
		letterToNumber($main::row);							# Changes row letters to numbers using the letterToNumber subroutine
		$info[1] = $main::row;								# reinsert the row, now as a number into the @info array
		my $codec = "$info[0],$info[1],$info[2]";			# define a universal identifier plate,row,column
			
		if (exists $main::oldscore{$codec}) {				# if the codec name exists in the oldscore hash
			my $oscore=$main::oldscore{$codec};				# get the original score for this orf
			if (($oscore==$oscore*1)&&($oscore ne "")) {
				trim($oscore);								# remove any flanking whitespace
				chomp($oscore);
				my $score=$main::newscore{$codec};			# get the new score for this orf
				chomp($score);
				$entry = $entry . "$score\t$oscore\n";		# add both to the results array
			}
		} else {
			$entry = $entry . "\t\n";						# add two tabs to the results array
		}
		push (@TEMPRESULTS, $entry);	
	}
	@main::RESULTS=@TEMPRESULTS;
}

#################################### sortlabelsave ########################################
# sort the results, add labels and save
sub sortlabelsave {
	
	@main::RESULTS = sort {									# sort the results array based on plate, column, row
		(split '\t', $a)[0] <=> (split '\t', $b)[0]||(split '\t', $a)[1] cmp (split '\t', $b)[1]||(split '\t', $a)[2] <=> (split '\t', $b)[2];		
	} @main::RESULTS;

	my $ndata="Plate\tRow\tColumn\tORF\t";					# define the initial column headings
	foreach(@main::NAMES) {									# for each of the names in the original data
		trim($_);
		chomp($_);
		my $l1=$_ . " corrected scores";					# make a heading for the new data
		my $l2=$_ . " original scores";						# make a heading for the old data
		$ndata=$ndata . "$l1\t$l2\t";						# add the new data to the headings entry
	}
	$ndata="$ndata\n";										# add a newline to the end of the headings
	unshift (@main::RESULTS,$ndata);						# add the column headings to the start of the results
	
	my $handle = "smoothed";
	my $filename="smoothed.tab";
	open($handle, ">$filename") or die "results file cannot be opened. $!";	# write output file
	print $handle "@main::RESULTS";							# "smoothed.tab"
}

#################################### letterToNumber #######################################
# converts row letters to numbers
sub letterToNumber {

	my $row=$main::row;									
	$row =~ tr/[ABCDEFGHI]/[123456789]/;					# use a simple translate for single digits

	if ( $row eq "J"){ $row = 10; }							# individually change the other letters
	elsif( $row eq "K" ) { $row = 11; }
	elsif( $row eq "L" ) { $row = 12; }
	elsif( $row eq "M" ) { $row = 13; }
	elsif( $row eq "N" ) { $row = 14; }
	elsif( $row eq "O" ) { $row = 15; }
	elsif( $row eq "P" ) { $row = 16; }
	$main::row=$row;
	return($row);
}

####################################	TRIM      #########################################
#subroutine to trim off the white space and tabs from both ends of each string
sub trim {
	$_[0]=~ s{^(.*)\t(.*?)$}{$1$2}gm;
	$_[0]=~ s/^\s+//;				
	$_[0]=~ s/\s+$//;
	return;
}