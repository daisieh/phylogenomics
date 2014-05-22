use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(make_label_lookup);

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my ($inputfile, $labelfile, $outfile, $simplename, $dnacode) = 0;
GetOptions ('input|samples=s' => \$inputfile,
            'labels|names=s' => \$labelfile,
            'simple' => \$simplename,
            'dna|code' => \$dnacode,
			'outfile=s' => \$outfile) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($inputfile && $labelfile) {
    pod2usage(-msg => "Error: an option was misspecified\n", -exitval => 2);
}


open FH, "<", $inputfile;
my @input_lines = <FH>;
close FH;

my $out_fh;
if ($outfile) {
	print "opening $outfile...\n";
	open $out_fh, ">", $outfile or die "couldn't open $outfile\n";
} else {
	$out_fh = STDOUT;
}

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
			$key = get_key($key, $dnacode);
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
			$key = get_key($key, $dnacode);
			my $label = $key;
			if (exists $labels->{$key}) {
				$label = $labels->{$key};
			}
			$line .= "$head$label$tail";
			$templine = $remainder;
		}
		$line .= "$templine\n";
	} elsif ($line =~ /^>(.+?)(\s.*|$)/) {
		# it is a fasta-type line
		my ($key, $remainder) = ($1,$2);
		$key = get_key($key, $dnacode);
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		chomp $remainder;
		if ($simplename) {
			$remainder = "";
		}
		$line = ">$label$remainder\n";
	} elsif ($line =~ /^(\s*\d+\s+)(.+?)(,*)$/) {
		# it is a TRANSLATE-type line (numbered list of taxa)
		print "translate $line";
		my ($beginning, $key, $remainder) = ($1,$2,$3);
		$key = get_key($key, $dnacode);
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		chomp $remainder;
		$line = "$beginning$label$remainder\n";
	} elsif ($line =~ /^(.+?)(\s+.+)$/) {
		# it is a generic NEXUS/phylip-type line
		my ($key, $remainder) = ($1,$2);
		$key = get_key($key, $dnacode);
		my $label = $key;
		if (exists $labels->{$key}) {
			$label = $labels->{$key};
		}
		chomp $remainder;
		$line = "$label$remainder\n";
	}
 	print $out_fh "$line";
 	$i++;
}

close $out_fh;

sub get_key {
	my $key = shift;
	my $dnacode = shift;

	if ($dnacode) {
		$key =~ s/.*DNA(\d+).*/$1/;
	}

	return $key;
}
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
