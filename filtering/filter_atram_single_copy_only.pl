#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec qw (rel2abs);
use File::Temp qw (tempfile);
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(parse_fasta);

my $validatefile = 0;
my $contigdir = 0;
my $outdir = 0;

GetOptions ('validatefile=s' => \$validatefile,
            'contigdir=s' => \$contigdir,
            'outdir=s' => \$outdir,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}


$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);


my $singlecopy_contigs = {};
open VALFH, "<", $validatefile or die "couldn't open validatefile $validatefile";
foreach my $line (<VALFH>) {
#Potri.001G166200.1	1e-90	457	810	81.18	self	4.22_len_1868_cov_11.1
	if ($line =~ /(.+?)\t.+?\tself\t(.*)$/) {
		$singlecopy_contigs->{$1} = $2;
	}
}
close VALFH;
print "found " . (keys %$singlecopy_contigs) . " single copy genes\n";

open OUTFH, ">", File::Spec->catfile($outdir, "genelist.txt");
foreach my $key (keys %$singlecopy_contigs) {
	print "finding $singlecopy_contigs->{$key} in $key.best.fasta\n";
	my $name = $key;
	if ($key =~ /(.+)\.\d$/) {
		$name = $1;
	}
	print OUTFH "$name\n";
	my $bestfile = File::Spec->catfile($contigdir, "$key.best.fasta");
	my ($contigs, $contigarray) = parse_fasta($bestfile);
	my $fh_name = File::Spec->catfile($outdir, "$name.fasta");
	open FH, ">", $fh_name or die "couldn't create file $fh_name";
	print FH ">$name\n$contigs->{$singlecopy_contigs->{$key}}\n";
	close FH;
	print "writing $name to $fh_name\n";
}

close OUTFH;
