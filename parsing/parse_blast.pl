#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Blast qw (parse_xml compare_hsps compare_regions);
use Genbank qw (parse_genbank write_features_as_fasta);
use Subfunctions qw (parse_fasta);
use Data::Dumper;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $help = 0;
my $outfile = "";
my $gbfile = "";
my $fastafile = "";


GetOptions ('reference=s' => \$gbfile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}


if ($gbfile !~ /\.gb$/) {
	print "reference file needs to be a fully annotated Genbank file.\n";
	exit;
}

my $gene_array = parse_genbank($gbfile);

open FAS_FH, ">", "$gbfile.fasta";
print FAS_FH write_features_as_fasta ($gene_array);
close FAS_FH;
my ($ref_hash, $ref_array) = parse_fasta ("$gbfile.fasta", 1);

print "running blastn -query $fastafile -subject $gbfile.fasta -outfmt 5 -out $outfile.xml -word_size 10\n";
system("blastn -query $fastafile -subject $gbfile.fasta -outfmt 5 -out $outfile.xml -word_size 10");

print "parsing results\n";

my $hit_array = parse_xml ("$outfile.xml");
my $hits = {};
foreach my $hit (@$hit_array) {
	my $subject = $hit->{"subject"}->{"name"};
	my $query = $hit->{"query"}->{"name"};
	my @hsps = sort compare_hsps @{$hit->{"hsps"}};
	my $best_hit = shift @hsps;
	$hits->{$subject}->{"hsp"} = $best_hit;
	$hits->{$subject}->{"slen"} = $hit->{"subject"}->{"length"};
	if ($best_hit->{"hit-from"} < $best_hit->{"hit-to"}) {
		$hits->{$subject}->{"orientation"} = 1;
	} else {
		$hits->{$subject}->{"orientation"} = -1;
	}
}
open OUTFH, ">", "$outfile.regions" or die "couldn't create $outfile";
my $prev_loc = 0;
my @sorted_ref_array = sort compare_regions @$ref_array;
foreach my $subj (@sorted_ref_array) {
	$subj =~ s/\t.*$//;
	if (!(exists $hits->{$subj}->{"hsp"})) {
		next;
	}
	my $adjusted_hseq = $hits->{$subj}->{"hsp"}->{"hseq"};
	$adjusted_hseq =~ s/-//g;
	my $hlen = length $adjusted_hseq;
	print OUTFH "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\n";
}

close OUTFH;



__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

parse_blast -reference genbank.gb -fasta fastafile [-outputfile output_file]

=head1 OPTIONS

  -fastafile:       fasta sequence to blast
  -reference:       genbank file to blast against
  -outputfile:      name of output file

=head1 DESCRIPTION


=cut
