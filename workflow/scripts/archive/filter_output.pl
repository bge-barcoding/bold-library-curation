use strict;
use warnings;
open my $dat, '<', 'draft_filtered_informations.csv' or die $!;
my$count=0;
while(<$dat>)
	{
		my$line=$_;
		my@array=split(/	/,$line);

		my$ID=$array[0];
		my$species=$array[1];
		my$real_species=$array[2];
		my$ident_rank=$array[3];
		my$voucher_type=$array[4];
		my$seqlenght=$array[5];
		my$ambiguities=$array[6];
		my$museum_id=$array[7];
		my$institut=$array[8];
		my$identifier=$array[9];
		my$country=$array[10];
		my$county=$array[11];
		next unless $seqlenght =~/^\d+/;
		if($seqlenght>500)
			{
				print "$seqlenght\n" unless $seqlenght=~/\d+/;
				$count++;
			}	
	}
print $count;	