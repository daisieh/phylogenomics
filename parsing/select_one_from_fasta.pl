use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "bioperl_subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $seq_name, $help) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
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

open FH, "<", "$fastafile" or die "couldn't open $fastafile";
my $line = readline FH;

while ($line !~ /$seq_name/) {
	$line = readline FH;
}

open RESULT_FH, ">", "$resultfile"."\.fasta";
print RESULT_FH $line;
$line = readline FH;

while ($line !~ />/) {
	print RESULT_FH $line;
	$line = readline FH;
}

close RESULT_FH;
close FH;

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
