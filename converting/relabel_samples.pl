require "subfuncs.pl";

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my ($inputfile, $labelfile, $outfile) = 0;
GetOptions ('input|samples=s' => \$inputfile,
            'labels|names=s' => \$labelfile,
			'outfile=s' => \$outfile) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($inputfile && $labelfile) {
    pod2usage(-msg => "Error: an option was misspecified\n", -exitval => 2);
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
my $i=0;
foreach my $line (@input_lines) {
	if ($line =~ /([\(\,]+)(.*?)([\:\,\)])(.*)/) {
		# we're in a newick tree
		my $templine = $line;
		$line = "";
		while ($templine =~ /^(.*?[\(\,]+)(.*?)([\:\,\)])(.*)/) {
			my ($head, $key, $tail, $remainder) = ($1, $2, $3, $4);
			my $label = $key;
			if (exists $labels->{$key}) {
				$label = $labels->{$key};
			}
			$line .= "$head$label$tail";
			$templine = $remainder;
		}
		$line .= "$templine\n";
		$templine = $line;
		$line = "";
		while ($templine =~ /^(.+?[\(\,]+.*?\,)(.*?)([\:\)])(.*)/) {
			my ($head, $key, $tail, $remainder) = ($1, $2, $3, $4);
			if ($key =~ /(\(.+\,)(.+)/) {
				$head = $head . $1;
				$key = $2;
			}
			my $label = $key;
			if (exists $labels->{$key}) {
				$label = $labels->{$key};
			}
			$line .= "$head$label$tail";
			$templine = $remainder;
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
	} elsif ($line =~ /^>(.+?)(\s+.*)$/) {
		# it is a fasta-type line
		my $key = $1;
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		$line = ">$label" . $2 . "";
	}
 	print $out_fh "$line";
 	$i++;
}

close $out_fh;

__END__

=head1 NAME

slice_fasta_file

=head1 SYNOPSIS

relabel_samples -input input_file -labels label_file [-outfile output_file]

=head1 OPTIONS

  -input|samples:    file containing names to replace
  -labels|names:     tab-delimited lookup file of names and their replacements
  -outfile:          output file name

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut
