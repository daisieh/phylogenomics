#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Blast qw (parse_xml);
use Genbank qw (parse_genbank parse_gene_array_to_features set_sequence get_sequence);
use Subfunctions qw (parse_fasta blast_to_genbank align_regions_to_reference align_hits_to_ref);
use Data::Dumper;

my $help = 0;
my $outfile = "";
my $gbfile = "";
my $fastafile = "";
my $orgname = "";
my $samplename = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('reference|gb|genbank=s' => \$gbfile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'organism=s' => \$orgname,
			'sample=s' => \$samplename,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($gbfile !~ /\.gb$/) {
	print "reference file needs to be a fully annotated Genbank file.\n";
	exit;
}

if ($outfile eq "") {
	$outfile = "output";
}

if ($gbfile eq "") {
	print "need to supply a Genbank reference file (-reference).\n";
	exit;
}

if ($fastafile eq "") {
	print "need to supply a fasta-formatted plastome sequence to annotate (-fasta).\n";
	exit;
}

my ($fastahash, $fastaarray) = parse_fasta($fastafile);

# there should be only one key, so just one name.
if ($samplename eq "") {
	$samplename = @$fastaarray[0];
}

my $queryseq = $fastahash->{@$fastaarray[0]};

# write out $outfile.regions
my ($result_hash, $result_array) = blast_to_genbank ($gbfile, $fastafile);

my @finished_array = ();
foreach my $subj (@$result_array) {
	if (exists $result_hash->{$subj}->{'hsps'}) {
		print "$subj\n" . Dumper (align_hits_to_ref ($result_hash->{$subj}, $queryseq));

	} else {
		push @finished_array, $subj;
	}
}

my $genbank_header = "$samplename [organism=$orgname][moltype=Genomic DNA][location=chloroplast][topology=Circular][gcode=11]";
open FASTA_FH, ">", "$outfile.fsa";
print FASTA_FH ">$genbank_header\n$queryseq\n";
close FASTA_FH;

my $gene_array = align_regions_to_reference ($result_hash, \@finished_array, $gbfile);

open TBL_FH, ">", "$outfile.tbl";
print TBL_FH Genbank::write_sequin_tbl ($gene_array, $genbank_header);
close TBL_FH;

# need to annotate inverted repeats

print "finding inverted repeats\n";
# my ($fh, $refblast) = tempfile();
my $refblast = "temp.xml";
system("blastn -query $fastafile -subject $fastafile -outfmt 5 -out $refblast -evalue 1e-200");

my $self_array = Blast::parse_xml ("$refblast");
my @irs = ();
foreach my $hit (@$self_array) {
	my @hsps = sort Blast::sort_regions_by_start @{$hit->{"hsps"}};
	foreach my $hsp (@hsps) {
		# only look at identical pieces that are smaller than the entire reference
		my $querylen = $hsp->{"query-to"} - $hsp->{"query-from"};
		if (($querylen < 50000) && ($querylen > 10000)) {
			push @irs, $hsp;
		}
	}
}

if (@irs > 2) {
	die "Error! There seem to be more than two inverted repeats. Are you sure this is a plastome sequence?";
}
@irs = sort Blast::sort_hsps_by_query_start @irs;

# all of the IR information is in one of the IRs, so shift.
my $ir = shift @irs;
my $irb = $ir->{'query-from'} . ".." . $ir->{'query-to'};
my $ira = $ir->{'hit-to'} . ".." . $ir->{'hit-from'};

#      repeat_region   84922..112749
#                      /note="inverted repeat B"
#                      /rpt_type=inverted

__END__

=head1 NAME

plann.pl

=head1 SYNOPSIS

plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile [-organism "Genus species"] [-sample samplename]

=head1 OPTIONS

  -reference:       a well-annotated plastome reference sequence, in genbank file format
  -fastafile:       the plastome sequence to be annotated, in fasta format
  -outfile:         the output name (default is "output")
  -organism:        [optional: scientific name for Genbank annotation]
  -sample:          the name of the plastome sample (default is the name in the fasta file)

=head1 DESCRIPTION


=cut
