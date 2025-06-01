use strict;
use warnings;
use Getopt::Long qw(GetOptions); # for commandline parameters
my$file; # variable for the filename of the BOLD snapshot file
my%done; # hash that is used to check the status of a species name
my%europ; # hash that is used to check the location of a species name
my%true; # hash to assign the valid name to synonyms
#-------------------------------------------
#-------------------Usage-------------------
#-------------------------------------------
# This script extracts data from BOLD snapshot datasets 
# Usage: >read_snapshot_data.pl -f <BOLD snapshot file>
#-------------------------------------------
GetOptions 	(
	'file=s' => \$file, # define parameter handle 
	) 
	or die "Usage: $0 -f filename\n"; # report an input arror due to wrong usage
#-------------------------------------------
# read species names and synonyms from file
open my $dat,'<','all_specs_and_syn.csv' or die $!;
print "\ncollect species names\n";
while(<$dat>){
		my$line=$_; # define the variable $line for the current line 
		chomp $line; # remove line ending 
		my@names=split(/;/,$line); # create an array with all the names and synonyms
		foreach(@names) {
				my$spec=$_; # set the name as value for the variable $spec
				chomp$spec; # remove possible line endings
				$true{$spec}=$names[0]; # set first name in the array as valid name for all following names in the array
				$europ{$spec}=1; # define the has %europe as present for each name in the list
			}
	}
close $dat; # reset the filehandle
open $dat,'<',"$file" or die $!; # open BOLD spanshot *.tsv file and define the filehandle DAT
print "read BOLD data\n"; # show that the script has reached the second part of the process
unlink "filtered_informations.csv" if -e "filtered_informations.csv"; # delete the output file if it already exists 
open my $out,'>>','filtered_informations.csv' or die $!; # open the output file
print $out "sampleid	species	real_species	ident_rank	voucher_type	seqlength	ambiguities	museum_id	institut	identifier	country	BIN\n"; # print output header in the output file
my$count=0; # set the value of $count to 0
my$count_samples=0;
my$line_num=0;
while(<$dat>) {
		if($line_num==0) {
				$line_num=1;
				next;
			}
		my$line=$_; # the whole line that is currently being examined
		my@line=split(/\t/,$line);
		next unless defined $europ{$line[16]};
		my@data=read_bold_public($line); # using the subroutine read_bold_public on the current data 
		unless(defined $europ{$data[1]}){next} # skip line if the species is not present in the all_specs_and_syn.csv file
		my$progress=0; # setting the value for $progress to 0 
		foreach(@data) {
				my$val=$_; # current value (subroutine output)
				chomp$val if defined $val; # remove line ending if $val is defined
				print  $out "$val" if defined $val; # print data from subroutine output array
				if($progress==1) # indication for $array[1] which contains the species name
					{
						print $out "\t$true{$val}" # print valid name of species ($val) as seperate collumn in the output file
					}
				print  $out "\t"; # add separator to the output to end the collumn
				$progress++; # increase the value $progress
			}
		$count_samples++;	
		print  $out "\n"; # add line ending to start in a new line in the next loop
		$count++ unless defined $done{$true{$data[1]}}; # count species only once
		$done{$true{$data[1]}}=1; # define the hash %done for the current species name
	}
print "species:	$count\n"; # print out number of species (includes only valid species)
print "samples:	$count_samples\n"; # print out number of samples in the filtered data
#-------------------------------------------
#subroutine that codes all the wanted data in the BOLD snapshot dataset
sub read_bold_public {
		my$line=$_[0];

		my@code=split(/\t/,$line); # converting the table contents into an array

		my$species=$code[16]; # species name
		chomp$species;
		my$processid=$code[0]; # process ID
		my$sampleid=$code[1]; # sample ID
		my$gb_acs=$code[49]; # GenBank Accession nr
		my$specid=$code[2]; # species ID
		my$museumid=$code[3]; # museum ID
		my$inst=$code[5]; # institution were the sample/tissue/holotype etc is stored
		my$voucher_type=$code[19]; # type of sample (Holotype/Paratype/DNA/frozen etc.)
		my$identification=$code[7]; # Taxonomic identification of the specimen
		my$identifier=$code[18]; # who identified the sample
		my$ident_method=$code[47]; # identification method (morphology/molecular etc)
		my$col_date=$code[21]; # collection date
		my$BIN=$code[6];
		my$country=$code[35]; # country where the sample was collected
		chomp $country; # remove possible lineendings
		my$site=$code[40]; # collection site
		my$coord=$code[28]; # coordinates of the collection site
		my$ident_rank=$code[55]; # identification rank (species/family/genus etc)		
		my$nucl=$code[51]; # the sequence
		chomp $nucl; # remove possible lineendings
		my@nucl=split(//,$nucl);
		AGAIN:
		if($nucl[0]=~/^[^AGTC]$/) {
				shift @nucl; # remove first element in the list
				$nucl=join('',@nucl); # sequence string is equal to modified array		
				goto AGAIN # go back to check last array element
			}
		if($nucl[-1]=~/^[^AGTC]$/) {	
				pop @nucl; # remove last element in array
				$nucl=join('',@nucl); # redefine sequence string as joined array
				goto AGAIN # go back to check again
			}	
		my$seqlen=length $nucl; #length of the whole sequence
		my$counter=$nucl=~tr/A,G,T,C,-//; # number of non ambigous sites
		my$ambigu=$seqlen-$counter; # sequence length minus the non ambiguous sites is equal to the number of ambiguities
		return($sampleid,$species,$ident_rank,$voucher_type,$seqlen,$ambigu,$museumid,$inst,$identifier,$country,$BIN); #output value list
	}



