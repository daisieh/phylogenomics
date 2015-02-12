#!/usr/bin/env perl
use strict;
use File::Temp qw (tempfile tempdir);
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(split_seq reverse_complement meld_matrices);
use Blast qw (blast_to_ref);


if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $ref_file = 0;
my $align_file = 0;
my $out_file = 0;
my $separate = 0;
my $help = 0;
my $blast_file = "";
my $evalue = 10;
my $ref_out = 0;
my $no_blast = 0;
my $debug = 0;

GetOptions ('fasta|input=s' => \$align_file,
            'outputfile=s' => \$out_file,
            'reference=s' => \$ref_file,
            'separate' => \$separate,
            'blastfile=s' => \$blast_file,
            'evalue=f' => \$evalue,
            'include_ref' => \$ref_out,
            'no_blast|noblast' => \$no_blast,
            'debug' => \$debug,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($ref_file and $align_file and $out_file) {
    pod2usage(-msg => "Must specify reference file, alignment file, and output file.");
}

print $runline;

my @references = ();
my @refids = ();
open refFH, "<:crlf", $ref_file;
my $refseq = "";
my $refid = "";

while (my $line = readline refFH) {
	if ($line =~ />(.+)$/) {
		$refid = $1;
		$refid =~ s/\s/_/g;
		push @refids, $refid;
		if ($refseq ne "") {
			push @references, "$refseq";
			$refseq = "";
		}
	} else {
 		chomp $line;
		$refseq .= $line;
	}
}

push @references, "$refseq";

close refFH;

my $result_matrices = ();

my (undef, $tempreffile) = tempfile(UNLINK => 1);
if ($blast_file eq "") {
	(undef, $blast_file) = tempfile(UNLINK => 1);
}

if ($no_blast != 1) {
	system ("makeblastdb -in $ref_file -dbtype nucl -out $tempreffile.db");
	system ("blastn -task blastn -evalue $evalue -query $align_file -db $tempreffile.db -outfmt 5 -out $blast_file");
} else { debug ("no blast\n"); }
$result_matrices = blast_to_ref("$blast_file");

for (my $i=0;$i<@refids;$i++) {
	$result_matrices->{$refids[$i]}->{'reference'} = $references[$i];
}

if ($separate == 0) {
	$refid = join ("|", @refids);
} else {
	$refid = "reference";
}

for (my $i=0;$i<@refids;$i++) {
	my $key = $refids[$i];
	my $reflen = length($references[$i]);
	$result_matrices->{$key}->{$refid} = delete ($result_matrices->{$key}->{'reference'});
	foreach my $k (keys (%{$result_matrices->{$key}})) {
		# pad out the sequence at the end so that they're aligned for matrixmeld.
		${$result_matrices->{$key}}{$k} .= "-" x ($reflen - length(${$result_matrices->{$key}}{$k}));
	}
}

if ($separate == 0) {
	my ($res1, $res2) = meld_matrices(\@refids, $result_matrices);
	my %mastertaxa = %{$res1};
	my %regiontable = %{$res2};

	my $value = $mastertaxa{$refid};
	open (fileOUT, ">", "$out_file.fasta");
	if ($ref_out == 1) {
		print fileOUT ">$refid\n$value\n";
	}
	delete $mastertaxa{$refid};
	delete $mastertaxa{"length"};

	# currently printing in fasta format: perhaps add a flag to alter this?
	foreach my $key ( keys %mastertaxa ) {
		my $value = $mastertaxa{$key};
		print fileOUT ">$key\n$value\n";
	}
	close fileOUT;
} else {
	foreach my $refname (@refids) {
		open (fileOUT, ">", "$out_file.$refname.fasta");
		my $value = ${$result_matrices->{$refname}}{$refid};
		if ($ref_out == 1) {
			print fileOUT ">$refid\n$value\n";
		}
		delete ${$result_matrices->{$refname}}{$refid};
		delete ${$result_matrices->{$refname}}{"length"};
		foreach my $key ( keys (%{$result_matrices->{$refname}}) ) {
			my $value = ${$result_matrices->{$refname}}{$key};
			print fileOUT ">$key\n$value\n";
		}
		close fileOUT;
	}
}

__END__

=head1 NAME

blast_to_ref

=head1 SYNOPSIS

blast_to_ref -fasta fastafile -reference reffile -output outputfile

=head1 OPTIONS

  -fasta|input:     fasta file of aligned sequences.
  -reference:       fasta file with sequences of interest.
  -outputfile:      output file name.
  -blastfile:	    optional: if specified, put BLASTN output in this file.
  -evalue:	        optional: if specified, sets evalue threshold for BLASTN hits.
  -include_ref:     optional: if specified, includes the reference sequence in the fasta output.
  -no_blast:        optional: (use with -blastfile flag) uses the specified blast xml file as input (for debugging)
  -debug:           optional: spews lots of information to STDERR.

=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

