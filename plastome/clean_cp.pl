#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../parsing/";
use Blast qw (parse_xml revcomp_hsp);
use Genbank qw (parse_genbank get_sequence);
use lib "$FindBin::Bin/../";
use Subfunctions qw (parse_fasta reverse_complement split_seq find_sequences consensus_str);
use File::Temp qw (tempfile);
use Data::Dumper;
use YAML::Tiny;

my $help = 0;
my $outfile = "";
my $reffile = "";
my $contigfile = "";


GetOptions ('reffile=s' => \$reffile,
			'contigfile=s' => \$contigfile,
			'outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $refseq = "";
if ($reffile =~ /\.gb$/) {
	my $gb = parse_genbank($reffile);
	$refseq = get_sequence($gb);
} else {
	my ($ref_hash, $ref_array) = parse_fasta($reffile);
	$refseq = $ref_hash->{@$ref_array[0]};
}
my $reflen = length ($refseq);

my ($fh, $to_align) = tempfile();
# my $to_align = "$outfile.toaln.fasta";
# open $fh, ">", $to_align;

print $fh ">reference\n$refseq\n";

# confirm ordering of outputs

# (undef, my $aligned_file) = tempfile(OPEN => 0);
my $aligned_file = "$outfile.aln.fasta";

my ($contig_hash, $contig_array) = parse_fasta($contigfile);

# output all of the final_contigs as one sequence for alignment.
print $fh ">$contigfile\n";

foreach my $c (@$contig_array) {
	print $fh $contig_hash->{$c} . "\n";
}
close $fh;

system ("mafft  --retree 2 --inputorder $to_align > $aligned_file");

my ($aln_hash, undef) = parse_fasta($aligned_file);

my $aligned_ref = $aln_hash->{"reference"};
my $aligned_seq = $aln_hash->{"$contigfile"};

# trim the beginning of the sequence, if necessary
if ($aligned_ref =~ /^(-+)/) {
	my $offset = length $1;
	(undef, $aligned_seq, undef) = split_seq ($aligned_seq, $offset + 1, length $aligned_seq);
}

# trim the end of the sequence, if necessary
if ($aligned_ref =~ /(-+)$/) {
	my $offset = length $1;
	(undef, $aligned_seq, undef) = split_seq ($aligned_seq, 1, (length $aligned_seq) - $offset);
}

# remove remaining gaps
$aligned_seq =~ s/-//g;

open OUTFH, ">", "$outfile.cleaned.fasta";
print OUTFH ">$contigfile\n$aligned_seq\n";
close OUTFH;

__END__

=head1 NAME

clean_cp.pl

=head1 SYNOPSIS

clean_cp.pl [-reffile reffile] [-contigfile contigfile] [-outputfile outputfile]

=head1 OPTIONS

  -reffile:         fasta file of reference plastome
  -contigfile:      fasta file of putative cp contigs
  -outputfile:      name of output file

=head1 DESCRIPTION

Cleans a draft plastome, such as the draft.fasta output from contigs_to_cp.pl

=cut
