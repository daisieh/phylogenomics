use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Pod::Usage;
use File::Basename;
use Getopt::Long;

require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($fa_file, $gb_file, $out_file, $help, $multiple) = 0;
GetOptions ('fasta=s' => \$fa_file,
            'genbank|gb_file=s' => \$gb_file,
            'outputfile=s' => \$out_file,
            'help|?' => \$help,
            'multiple!' => \$multiple) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns = @{parse_aln_into_genes($whole_aln, $gb_file,"CDS")};

if ($multiple) {
	foreach my $aln (@gene_alns) {
		my $gene_name = $aln->description();
		open my $gene_file, ">", "$out_file/$gene_name.fasta";
		foreach my $seq ($aln->each_seq()) {
			my $name = $seq->id();
			print $gene_file ">$name\n";
			print $gene_file $seq->seq() . "\n";
		}
		close $gene_file;
	}
} else {
	open my $gene_file, ">", "$out_file.fasta";

	foreach my $aln (@gene_alns) {
		my $gene_name = $aln->description();
		foreach my $seq ($aln->each_seq()) {
			my $name = $seq->id() . "_$gene_name: " . $seq->length();
			print $gene_file ">$name\n";
			print $gene_file $seq->seq() . "\n";
		}
	}
	close $gene_file;
}

__END__

=head1 NAME

parse_fasta_to_genes

=head1 SYNOPSIS

parse_fasta_to_genes [-fasta fa_file] [-genbank gb_file] [-outputfile output_file] [-multiple]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -outputfile:      name of either output file or directory name, if -multiple is specified
  -multiple:		output separate fasta files for each CDS (-outputfile is directory name)

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut

