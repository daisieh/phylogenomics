#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Blast qw (parse_xml);
use lib "$FindBin::Bin/../";
use Subfunctions qw (parse_fasta);
use Data::Dumper;

my $help = 0;
my $outfile = "";
my $reffile = "";
my $fastafile = "";


GetOptions ('reffile=s' => \$reffile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print "running blastn\n";
system("blastn -query $fastafile -subject $reffile -outfmt 5 -out $outfile.xml -word_size 10");

print "parsing results\n";

my ($ref_hash, $ref_array) = parse_fasta($reffile);
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

open OUTFH, ">", $outfile or die "couldn't create $outfile";

foreach my $subj (@$ref_array) {
	$subj =~ s/\t.*$//;
	if (!(exists $hits->{$subj}->{"hsp"})) {
		print "$subj: no hits\n";
		next;
	}
	my $adjusted_hseq = $hits->{$subj}->{"hsp"}->{"hseq"};
	$adjusted_hseq =~ s/-//g;
	my $hlen = length $adjusted_hseq;
	print OUTFH "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\n";
	if ($hlen == $hits->{$subj}->{"slen"}) {
		if ($hits->{$subj}->{"orientation"} > 0) {
			print "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\n";
		} else {
			print "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\tinverted match\n";
		}
	} else {
		print "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\tpartial match: $hlen ne $hits->{$subj}->{slen}\n";
	}
}

close OUTFH;

sub compare_hsps {
	my $score = $b->{"bit-score"} - $a->{"bit-score"};
	if ($score == 0) {
		my $b_direction = ($b->{"query-to"} - $b->{"query-from"})/($b->{"hit-to"} - $b->{"hit-from"});
		my $a_direction = ($a->{"query-to"} - $a->{"query-from"})/($a->{"hit-to"} - $a->{"hit-from"});
		if ($b_direction > $a_direction) {
			$score = 1;
		} elsif ($a_direction < $b_direction) {
			$score = -1;
		} else {
			$score = 0;
		}
	}
	return $score;
}


__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

First run: blastn -query comparison.fasta -subject reference.fasta -outfmt 3 -out blast_file
Then run: parse_blast [-blast blast_file] [-outputfile output_file]

=head1 OPTIONS

  -blast:           "outfmt 3" formatted blastn output
  -outputfile:      name of output file

=head1 DESCRIPTION

Parses an "outfmt 3" formatted blastn file to generate a list of regions to be used in
Genbank annotations.

=cut
