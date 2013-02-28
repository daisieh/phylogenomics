require "subfuncs.pl";

my $fastafile = shift;
my $labelfile = shift;

my $fa_aln = make_aln_from_fasta_file ($fastafile, 0);

my $labels;
$labels = make_label_lookup ($labelfile);

foreach my $seqio ($fa_aln->each_seq) {
	my $key = $seqio->id();
	my $label = $key;
	if (exists $labels->{$key}) {
		$label = $labels->{$key};
	}
	print ">$label\n";
	print $seqio->seq() . "\n";
}
