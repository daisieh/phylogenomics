#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use FindBin;
use File::Spec;
use File::Path qw (make_path);
use lib "$FindBin::Bin";
use Blast qw (parse_xml sort_hsps_by_match);
use lib "$FindBin::Bin/..";
use Subfunctions qw (parse_fasta);
use Data::Dumper;

my $help = 0;
my $genefile = "";
my $outdir = "";
my $refdir = "";
my $querydir = "";
my $xmldir = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('reference=s' => \$refdir,
			'genes=s' => \$genefile,
			'query=s' => \$querydir,
			'output=s' => \$outdir,
			'xml=s' => \$xmldir,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

$outdir = File::Spec->rel2abs($outdir);
unless (-d $outdir) {
	make_path($outdir);
}

my @genes = ();

open FH, "<", $genefile;
while (my $line = readline FH) {
	chomp $line;
	push @genes, $line;
}

foreach my $gene (@genes) {
	my $queryfile = File::Spec->catfile ($querydir, "$gene.fasta");
	my $reffile = File::Spec->catfile ($refdir, "$gene.fasta");
	my $outfile = File::Spec->catfile ($outdir, "$gene");

	my $blastfile = "$outfile.xml";
	if ($xmldir eq "") {
		my $cmd = "blastn -query $queryfile -subject $reffile -outfmt 5 -out $blastfile -word_size 10";
		print "running $cmd\n";
		system($cmd);
	} else {
		$blastfile = File::Spec->catfile ($xmldir, "$gene.xml");
	}

	print "parsing results for $blastfile\n";
	my $hit_array = parse_xml ($blastfile);

	my ($ref_hash, $ref_array) = parse_fasta($reffile);

	my $hits = {};
	foreach my $hit (@$hit_array) {
		my $subject = $hit->{"subject"}->{"name"};
		my $query = $hit->{"query"}->{"name"};
		my @hsps = sort { sort_hsps_by_match ($a, $b) } @{$hit->{"hsps"}};
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
		print OUTFH "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\t$hits->{$subj}->{'hsp'}->{'identity'}\t$hits->{$subj}->{'slen'}\n";
		my $pident = ($hits->{$subj}->{'hsp'}->{'identity'} / $hits->{$subj}->{"slen"}) * 100;
		print "\t percent identity is $hits->{$subj}->{hsp}->{'identity'} / $hits->{$subj}->{'slen'} = $pident\n";
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
}

print "finished blast_list\n";
__END__

=head1 NAME

blast_list -gene genelist.txt -ref refdir/ -query querydir/ -out outputdir/

=head1 SYNOPSIS

GetOptions ('reference=s' => \$refdir,
			'genes=s' => \$genefile,
			'query=s' => \$querydir,
			'output=s' => \$outdir,
			'xml=s' => \$xmldir,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

=head1 DESCRIPTION

For a list of genes, takes the corresponding ref sequence (in fasta form, parsed out from a GFF), the corresponding query sequence that you want to match to the ref (in fasta form), and outputs the blastn result and a list of best matching regions for each of the reference pieces.

=cut
