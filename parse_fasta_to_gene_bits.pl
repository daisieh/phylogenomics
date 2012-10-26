use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
require "subfuncs.pl";

my ($fa_file, $gb_file, $out_file, $help) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'genbank|gb_file=s' => \$gb_file,
            'outputfile=s' => \$out_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ((@ARGV == 0) || ($help)) {
    pod2usage(-verbose => 1);
}

print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "parse_fasta_to_genes.pl genbank_file fasta_file out_file\n";
my $gb_file   = shift or die $usage;
my $fa_file = shift or die $usage;
my $out_file = shift or die $usage;

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns;

my $result_str = get_locations_from_genbank_file ($gb_file, "CDS");
my @gene_bits = split (/\n/,$result_str);

# we don't need the last line; it just says the size of the gb source sequence
while ({pop @gene_bits} =~ /source/) {
}

foreach my $gene_bit (@gene_bits) {
    $gene_bit =~ /(.+?)\t(.+?)\t(.+)\s*.*/;
    my $name = $1;
    my $start = $2;
    my $end = $3;
    my $curr_slice = $whole_aln->slice($start, $end);
    $curr_slice->description($name);
    push @gene_alns, $curr_slice;
}

open my $gene_file, ">$out_file.fasta";

foreach my $aln (@gene_alns) {
	my $gene_name = $aln->description();
	foreach my $seq ($aln->each_seq()) {
		my $name = $seq->id() . "_$gene_name: " . $seq->length();
		print $gene_file ">$name\n";
		print $gene_file $seq->seq() . "\n";
	}
}
close $gene_file;

__END__

=head1 NAME

parse_fasta_to_gene_bits

=head1 SYNOPSIS

parse_fasta_to_gene_bits [-fasta fa_file] [-genbank gb_file] [-outputfile output_file]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -outputfile:      prefix of output files

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut

