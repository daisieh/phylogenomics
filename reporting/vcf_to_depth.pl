use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(sample_list vcf_to_depth);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $samplefile = 0;
my $help = 0;
my $outfile = "";

GetOptions ('samples|input|vcf=s' => \$samplefile,
            'outputfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my @samples = @{sample_list ($samplefile)};
print "samples ". join (",",@samples). "\n";
if ($samplefile =~ /(.*?)\.vcf$/) {
	@samples = ($1);
}

foreach my $sample (@samples) {
	$name = basename($sample);
	print "processing $name...\n";
	my $vcf_file = $sample . ".vcf";
	if ($outfile ne "") {
		vcf_to_depth ($vcf_file,$outfile);
	} else {
		vcf_to_depth ($vcf_file);
	}
}


__END__

=head1 NAME

vcf_to_depth

=head1 SYNOPSIS

vcf_to_depth -samplefile -output [-threshold]

=head1 OPTIONS

  -samples|input|vcf:   name of sample or list of samples to convert
  -outputfile:      optional: prefix of output fasta file

=head1 DESCRIPTION

Generates a summary file with the read depths of each position in the inputted vcf file(s).

=cut

