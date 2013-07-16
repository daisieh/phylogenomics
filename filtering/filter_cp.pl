#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $usage = "perl " . basename($0);
$usage .=	" <fastafile> <cp.fasta> <resultfile>\n\n";
print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $reffile, $resultfile, $task) = 0;
my $evalue = 10;
my $run_blast = 1; # run blast by default
my $task = "keep-hits";
my $task_switch;
GetOptions ('fasta=s' => \$fastafile,
            'reference=s' => \$reffile,
            'outputfile=s' => \$resultfile,
            'task:s' => \$task,
            'discard-hits!' => \$task_switch,
            'evalue:i' => \$evalue,
            'blast!' => \$run_blast) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);


if ($task_switch) {
    $task = "discard-hits";
}

unless (($fastafile && $reffile && $resultfile) || (!$run_blast && $fastafile && $resultfile)) {
    my $msg = qq{Error: an option was mis-specified:
    fasta=$fastafile
    reference=$reffile
    outputfile=$resultfile
    task=$task
    };

    pod2usage(-msg => $msg, -exitval => 2);
}

my $filterfile = "$resultfile.filtered";

if ($run_blast) {
    print "blasting...\n";

    my @blastargs = ("blastn", "-evalue", "$evalue", "-subject", "$reffile", "-query", "$fastafile", "-outfmt", "6 qseqid bitscore evalue", "-out", "$filterfile");
    system (@blastargs);

    print "done (results in $filterfile)\n";
} else {
    print "not running blast, using results from $filterfile\n";
}

print "filtering...\n";

open fastaFH, "<$fastafile" or die "couldn't open fasta file";
open my $filterFH, "<$filterfile" or die "error opening filtered seqs file";
open my $outfile, ">$resultfile.fasta" or die "couldn't create result file";
truncate $outfile, 0;

if (eof($filterFH)) {
    print "no sequences blasted to reference.\n";
    exit;
}

my $filter_seq = find_next_seq (0, $filterFH);
my $fasta_line = readline fastaFH;
while (($fasta_line ne "") && $filter_seq) {
    #first line:
    my $first_header = "$fasta_line";
    $fasta_line = readline fastaFH;
    #second line:
    my $sequence = "$fasta_line";
    $fasta_line = readline fastaFH;
    while (($fasta_line !~ m/>.*/)) {
        $sequence .= "$fasta_line";
        $fasta_line = readline fastaFH;
        if (eof(fastaFH)) {
            last;
        }
    }
    if ($first_header =~ m/$filter_seq/) {
        if ($task ne "discard-hits") {
            print $outfile "$first_header$sequence";
        }
        $filter_seq = find_next_seq ($filter_seq, $filterFH);
    } else {
        if ($task eq "discard-hits") {
            print $outfile "$first_header$sequence";
        }
    }
}

print "done\n";

close $filterFH;
close fastaFH;
close $outfile;

sub find_next_seq {
    my $curr_seq = shift;
    my $FH = shift;

    my ($filtering, $seq);
    while ($filtering = readline $FH) {
        $filtering =~ /(.+?)\s/;
        $seq = $1;
        if ($seq !~ m/$curr_seq/) {
            return $seq;
        }
    }
    return undef;
}

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

filter_cp [options]

=head1 OPTIONS

    -fasta:     fasta file to filter
    -reference: fasta file of the reference to be filtered against
    -task:      "discard-hits" means to discard hits to the reference
                "keep-hits" means to keep hits
    -evalue:    sets the evalue for blastn (default is 10)
    --blast/noblast:    runs blastn or uses previously generated filtered list (default is to run blastn)

=head1 DESCRIPTION

=cut
