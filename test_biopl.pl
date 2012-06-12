require "subfuncs.pl";
use Bio::SeqIO;
use constant CENTER_X => 600; # X coordinate of circle center500
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object1000
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

my $usage  = "test_biopl.pl \n";

my $gb_file = shift or die $usage;

my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
my $seq_object = $seqio_object->next_seq;

while ($seq_object) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			my $name = main_name_for_gb_feature($feat_object);
			my @locations = $feat_object->location->each_Location;
			my $loc_str = "";
			foreach $loc (@locations) {
				my $start = $loc->start;
				my $end = $loc->end;
				$loc_str .= "($start..$end)";
			}
			print "$name\t$loc_str\n";
		}
	}
	$seq_object = $seqio_object->next_seq;
}

