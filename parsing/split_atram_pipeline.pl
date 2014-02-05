#!/usr/bin/perl
use File::Spec;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions;

my $fastafile = shift;
my $outdir = shift;

if (! defined $outdir) { print "\n\tUsage: split_atram_pipeline.pl fastafile outdir\n\n"; exit; }

my ($seqs, $seqnames) = parse_fasta($fastafile);
$outdir = File::Spec->rel2abs($outdir);

open LISTFH, ">", File::Spec->catfile( $outdir, "fastalist.txt" );

foreach my $t (@$seqnames) {
	my $clean_name = $t;
	$clean_name =~ s/\|.*$//;
	my $seqfile = File::Spec->catfile( $outdir, "$clean_name.fasta" );
	print LISTFH "$clean_name\t$seqfile\n";
	open FH, ">", $seqfile;
	print FH ">$t\n$seqs->{$t}\n";
	close FH;
}

close LISTFH;
