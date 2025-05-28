#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my$filename;
my$subfamily;
my$db;
my%BINs;
GetOptions(
    'in=s' => \$filename, # input filename
    'sub=s' => \$subfamily, # comma separated list of families that need further splitting	
    'db=s' => \$db # defines if db files will be created at the end
) or die "define input\n";

#USAGE: perl createxml_family.pl -in <input file> -sub<subfamily, subfamily,...> 

# Check if required args are present
die "Input file is not defined. Use --in=FILENAME\n" unless defined $filename;
#die "Subfamily is not defined. Use --sub=SUBFAMILY\n" unless defined $subfamily; 

goto SKIP unless defined $subfamily;
my@subfamily=split(",",$subfamily);
my%subfam;
foreach(@subfamily){
	$subfam{$_}=1;
}
SKIP:

#my$filename=$ARGV[0];
my@blacklist=("SPECIES_ID","TYPE_SPECIMEN","SEQ_QUALITY","HAS_IMAGE","COLLECTORS","COLLECTION_DATE","COUNTRY","REGION","SECTOR","SITE","COORD","IDENTIFIER","ID_METHOD","INSTITUTION","PUBLIC_VOUCHER","MUSEUM_ID","nuc","elev_accuracy","primers_forward","primers_reverse");
my%notthisone;
my$name;
my%avail;
my%start;
my%family;
open(DAT,"<all_specs_and_syn.csv");
while(<DAT>){
	my$line=$_;
	chomp $line;
	my@array=split(";",$line);
	foreach(@array){
		$avail{$_}=1;
		#print "$_\n";
	}
}

print $filename;
print "next step\n";
if($filename=~/^(.*)\.\w+/){
	$name=$1;
}
foreach(@blacklist){
	$notthisone{$_}=1;
}
print $name;
open my $fh, '<', $filename or die "$filename: $!";

system "rm -r family_output" if -e "family_output";
system "mkdir family_output";

my $sub;
my %done;
my @header;
my $counter=0;
my $speccount=0;
my%tax;
my@BINlist;
while(my $line = <$fh> ) {
	chomp $line;
	my$rank=20;
	$speccount++;
	if($speccount==5000){
		#goto SKIP;
	}
	my @collumns= split ("\t",$line);

	next if $collumns[24]=~/^None$/;
	next if $collumns[9]=~/^None$/;	
	next unless defined $collumns[94];
	next unless $collumns[9]=~/.+/;
	next unless $collumns[24]=~/.+/;

	if(defined $subfam{$collumns[20]}){
		system "mkdir family_output/$collumns[$rank]" unless -e "family_output/$collumns[$rank]";
		open(OUT,">>family_output/$collumns[$rank]/$collumns[21].xml") unless $collumns[$rank] eq "family";	
		$rank=21;
		#print "$collumns[9] $collumns[21]\n";

		if(defined $BINs{$collumns[9]}){
			$BINs{$collumns[9]}="$BINs{$collumns[9]};$collumns[21]" unless defined $done{"$collumns[9]$collumns[21]"};
			$done{"$collumns[9]$collumns[21]"}=1;
			push(@BINlist,$collumns[9]);
			$tax{$collumns[9]}=$collumns[20];
		}
		else{
			$BINs{$collumns[9]}=$collumns[21];
			$done{"$collumns[9]$collumns[21]"}=1;
		}
	}
	else{
	open(OUT,">>family_output/$collumns[$rank].xml") unless $collumns[$rank] eq "family";
	}
	#unless(defined $avail{$collumns[25]}){
	#	next;
	#}
	unless ($collumns[$rank] eq "family"){
	print OUT "<records>\n" unless defined $start{$collumns[$rank]};
	$start{$collumns[$rank]}=1;
	}
	if($counter == 0)
		{
			#print OUT "<records>\n";
			@header= split("\t",$line);
			$counter++;
			next;
		}		
	if(defined $avail{$collumns[24]}){

	}	
	else
	{
		next;
	}
	$header[96]="nothing";	
	#print $line;
	print OUT "	<record>\n";
	#print OUT "		<Keep>0</Keep>\n";

	my$length=@collumns;
	
	my $count_header=0;
	foreach(@collumns) {

		my$var=$header[$count_header];
		my$value=$_;


		# Rename COLLECTORS to HAS_COLLECTOR
		if ($var eq 'bold_recordset_code_arr') {
	        $var = 'recordset_code_arr';
	    }
	    if($header[$count_header] eq "country\/ocean") {
			$var ="country_ocean";
		}
		if($header[$count_header] eq "province\/state") {
			$var ="province_state";
		}

		if($value=~/^(.*)&(.*)$/) {
			$value="$1 et $2";
			$value =~ s/&/and\;/g;
		}
		if($value=~/^(.*)\<(.*)$/) {

			$value =~ s/</ \;/g;
		}		
		#next if defined $notthisone{$var};
		unless(defined $notthisone{$var}){
		print OUT "		<$var>$value<\/$var>\n" if defined $var;
		}
		$count_header++;
	}

	#print OUT "		<additionalStatus><\/additionalStatus>\n";
	print OUT "	<\/record>\n";
	close OUT;

}
SKIP:
my@files=glob("family_output/*.xml");
foreach(@files){
open(OUT,">>$_");

print OUT "<\/records>\n";
close OUT;	

}

foreach(@subfamily){
print "\n$_\n";
my@files=glob("family_output/$_/*.xml");
foreach(@files){
	my$path=$_;
	print "\n$path\n";
	open(OUT,">>$path");

	print OUT "<\/records>\n";
	close OUT;	

print $_;
}
}

foreach(@BINlist){
	if ($BINs{$_}=~/;/){
		my$id=$_;
		open(OUT,">>family_output/$tax{$id}/subfamily_shared_BINs.csv");
		print OUT "$id;$BINs{$id};$tax{$id}\n" unless defined $done{"$id;$BINs{$id};$tax{$id}"};
		$done{"$id;$BINs{$id};$tax{$id}"}=1;
		close OUT;
	}
	
}
