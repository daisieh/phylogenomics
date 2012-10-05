#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $usage = "perl " . basename($0);
$usage .=	" <fastafile> <cp.fasta> <resultfile>\n\n";

my ($fastafile, $reffile, $resultfile, $task) = 0;
my $evalue = 10;
my $task = "keep-hits";
GetOptions ('fasta=s' => \$fastafile,
            'reference=s' => \$reffile,
            'task=s' => \$task,
            'evalue:i' => \$evalue) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($fastafile && $reffile && $resultfile) {
    my $msg = qq{Error: an option was mis-specified:
    fasta=$fastafile
    reference=$reffile
    task=$task
    };
    pod2usage(-msg => $msg, -exitval => 2);
}

my $filterfile = "$resultfile.filtered";

print "blasting...";

my @blastargs = ("blastn", "-evalue", "$evalue", "-subject", "$reffile", "-query", "$fastafile", "-outfmt", "6 qseqid bitscore evalue", "-out", "$filterfile");
system (@blastargs);

print "done (results in $filterfile)\n";

open my $F, "<$fastafile" or die "couldn't open fasta file";
my $fs = readline $F;

open my $this_file, "<$filterfile" or die "error opening filtered seqs file";

my $filtering = readline $this_file;
if (eof($this_file)) {
	print "no sequences blasted to reference.\n";
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
    if ($task eq "discard-hits") {
        last;
    } else {
        if ($filtering ne "%%%") {
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
}

print "done\n";

close $this_file;
close $F;
close $outfile;

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

=head1 DESCRIPTION

=cut
