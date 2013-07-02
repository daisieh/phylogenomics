use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Pod::Usage;
use File::Basename;
use Getopt::Long;

require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($fastafile, $out_file, $help, $consensus, $iupac) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$out_file,
            'help|?' => \$help,
            'iupac' => \$iupac,
            'consensus!' => \$consensus) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $out_fh = STDOUT;
if ($out_file) {
	open $out_fh, ">", $out_file;
}

my $fa_aln = make_aln_from_fasta_file ($fastafile, 0);

my $str = $fa_aln->consensus_iupac();

if ($consensus) {
	print $out_fh ">$fastafile\n$str\n";
} else {
	my $aln_length = $fa_aln->length();
	for (my $i=1; $i<=$aln_length; $i++) {
		$str =~ s/^(.)//;
		my $char = $1;
		if ($char !~ m/[agctAGCT-]/) {
			if ($iupac) {

			} else {
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
			}
		} else {
			$char = "";
		}
		print $out_fh "$i\t$char\n";
	}
}

if ($out_file) {
	close $out_fh;
}
