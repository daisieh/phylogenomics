require "subfuncs.pl";

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my ($inputfile, $labelfile, $outfile) = 0;
GetOptions ('input|samples=s' => \$inputfile,
            'labels|names=s' => \$labelfile,
			'outfile=s' => \$outfile) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($inputfile && $labelfile) {
    pod2usage(-msg => "Error: an option was mis-specified\n", -exitval => 2);
}


open FH, "<", $inputfile;
my $out_fh;
if ($outfile) {
	open $out_fh, ">", $outfile;
} else {
	$out_fh = STDOUT;
}
my @input_lines = <FH>;

my $labels;
$labels = make_label_lookup ($labelfile);

foreach my $line (@input_lines) {
	if ($line =~ /([\(\,]+)(.*?)([\:\,\)])(.*)/) {
		# we're in a newick tree
		my $templine = $line;
		$line = "";
		while ($templine =~ /^(.*?[\(\,]+)(.*?)([\:\,\)])(.*)/) {
			my $key = $2;
			my $label = $key;
			if (exists $labels->{$key}) {
				$label = $labels->{$key};
			}
			$line .= $1 . "$label" . $3;
			$templine = $4;
		}
		$line .= "$templine\n";
	} elsif ($line =~ /^(.+?)(\s+.+)$/) {
		# it is a NEXUS/phylip-type line
		my $key = $1;
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		$line = $label . $2 . "\n";
	} elsif ($line =~ /^>(.+?)(\s*.*)$/) {
		# it is a fasta-type line
		my $key = $1;
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		$line = $label . $2 . "\n";
	}

 	print $out_fh "$line";
}

close $out_fh;
