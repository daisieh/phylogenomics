require "subfuncs.pl";

my $fastafile = shift;

my $fa_aln = make_aln_from_fasta_file ($fastafile, 0);

my $str = $fa_aln->consensus_iupac();

# print ">$fastafile\n$str\n";

my $aln_length = $fa_aln->length();
for (my $i=1; $i<=$aln_length; $i++) {
	$str =~ s/^(.)//;
	my $char = $1;
	if ($char !~ m/[agctAGCT-]/) {
		$char =~ s/M/A\/C/i;
		$char =~ s/R/A\/G/i;
		$char =~ s/W/A\/T/i;
		$char =~ s/S/C\/G/i;
		$char =~ s/Y/C\/T/i;
		$char =~ s/K/G\/T/i;
		$char =~ s/V/A\/C\/G/i;
		$char =~ s/H/A\/C\/T/i;
		$char =~ s/D/A\/G\/T/i;
		$char =~ s/B/C\/G\/T/i;
		$char =~ s/N/A\/C\/G\/T/i;
	} else {
		$char = "";
	}
	print "$i\t$char\n";
}
