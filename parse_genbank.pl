use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Bio::AlignIO;

require "subfuncs.pl";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "test_biopl.pl in_file out_file result_dir\n";
my $in_file   = shift or die $usage;
my $out_file = shift or die $usage;

$in_file =~ /.*\.(.+)/;
my $in_format = $1;

$out_file =~ /.*\.(.+)/;
my $out_format = $1;

if (($out_format eq "fa") || ($out_format eq "fas") || ($out_format eq "fasta")) {
	$out_format = "fasta";
} elsif (($out_format eq "nex") || ($out_format eq "nexus")) {
	$out_format = "nexus";
} elsif (($out_format eq "gb") || ($out_format eq "genbank")) {
	$out_format = "genbank";
}
#my $seqio_object = Bio::SeqIO->new(-file => $in_file);
my $in = Bio::AlignIO->newFh(-file => $in_file);
#my $seq_object = $seqio_object->next_seq;

my $seq_out = Bio::AlignIO->newFh(
                              -file   => ">$out_file",
                              -format => $out_format,
                              );

print $seq_out $_ while <$in>;
# write each entry in the input file to the output file
# while (my $inseq = $seqio_object->next_seq) {
#     $seq_out->write_seq($inseq);
# }
