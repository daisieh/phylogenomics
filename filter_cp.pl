#!/usr/bin/perl
use strict;
use File::Basename;


my $usage = "perl " . basename($0);
$usage .=	" <fastafile> <cp.fasta> <resultfile>\n\n";

my $fastafile = shift or die "$usage";
my $cpfile = shift or die "$usage";
my $resultfile = shift or die "$usage";
my $filterfile = "$fastafile.filtered";

print "blasting...";

my @blastargs = ("blastn", "-evalue", "20", "-subject", "$cpfile", "-query", "$fastafile", "-outfmt", "6 qseqid bitscore", "-out", "$filterfile");
system (@blastargs);

print "done (results in $filterfile)\n";

open my $F, "<$fastafile" or die "couldn't open fasta file";
my $fs = readline $F;

open my $this_file, "<$filterfile" or die "error opening filtered seqs file";

my $filtering = readline $this_file;
if (eof($this_file)) {
	print "no sequences blasted to e. coli\n";
	$filtering = "%%%";
}

while ($filtering =~ m/^\s+$/) {
	$filtering = readline $this_file;
	if (eof($this_file)) {
		print "%%%";
		$filtering = "%%%";
		last;
	}
}

print "filtering...";

open my $outfile, ">$resultfile.fasta" or die "couldn't create result file";
truncate $outfile, 0;

while ($fs ne "") {
	#first line:
	my $first_header = "$fs";
	$fs = readline $F;

	#second line:
	my $sequence = "$fs";
	$fs = readline $F;
	while (($fs !~ m/>.*/)) {
		chomp $sequence;
		$sequence .= "$fs";
		$fs = readline $F;
		if (eof($F)) {
			last;
		}
	}
	if ($filtering ne "%%%") {
# 		print "checking $filtering...";
		$filtering =~ /(.+?)\t/;
		my $curr_seq = $1;

		if ($first_header =~ m/$curr_seq/) {
			$filtering = readline $this_file;
			print $outfile "$first_header$sequence";
			while ($filtering =~ m/$curr_seq/) {
				$filtering = readline $this_file;
				if ($filtering eq "") {
					$filtering = "%%%";
					last;
				}
			}
		}
	}
}

print "done\n";

close $this_file;
close $F;
close $outfile;
