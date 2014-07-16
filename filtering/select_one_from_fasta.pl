#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use File::Path qw (make_path);
use File::Spec qw (catfile rel2abs);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $seq_name, $help, $split_all) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'split_all' => \$split_all,
            'sequence=s' => \$seq_name,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

unless (($fastafile && $resultfile)) {
    my $msg = qq{Error: an option was mis-specified:
    fasta=$fastafile
    outputfile=$resultfile
    };

    pod2usage(-msg => $msg, -exitval => 2);
}

# make the output file an absolute path, just to be safe.
my $output_path = File::Spec->rel2abs($resultfile);

if ($split_all) {
	print "$output_path\n";
	# check to see if the path for the output_file exists; if not, create it.
	unless (-d $output_path) {
		make_path ($output_path);
	}
}

open fileIN, "<:crlf", "$fastafile" or die "couldn't open $fastafile";
my $input = readline fileIN;
my $taxonlabel = "";
my $sequence = "";
while (defined $input) {
	if ($input =~ /^>(.+)\s*$/) {
		print "$input";
		if ($sequence ne "") { #are we done with a sequence?
			if ($split_all) {
				open FH, ">", File::Spec->catfile($output_path, "$taxonlabel.fasta");
				print FH ">$taxonlabel\n$sequence\n";
				close FH;
			} elsif ($taxonlabel eq $seq_name) {
				last;
			}
		}
		$sequence = "";
		$taxonlabel = $1;
		$taxonlabel =~ s/\s+/_/g;
		$taxonlabel =~ s/_$//;
	} else {
		$input =~ /^\s*(.+)\s*$/;
		$sequence .= "$1\n";
	}
	$input = readline fileIN;
}

close fileIN;

if ($split_all) {
	open FH, ">", File::Spec->catfile($output_path, "$taxonlabel.fasta");
	print FH ">$taxonlabel\n$sequence\n";
	close FH;
} elsif ($taxonlabel eq $seq_name) {
	open FH, ">", "$resultfile.fasta";
	print FH ">$taxonlabel\n$sequence\n";
	close FH;
}


__END__

=head1 NAME

select_one_from_fasta

=head1 SYNOPSIS

select_one_from_fasta [-fasta fa_file] [-outputfile output_file] [-sequence seq_name]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -outputfile:      output file name
  -sequence:        name of sequence to be selected.

=head1 DESCRIPTION

Finds a single sequence from a fasta file and outputs to a separate file.

=cut
