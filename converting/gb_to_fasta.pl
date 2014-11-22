#!/usr/bin/env perl

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Genbank qw(parse_genbank get_sequence get_name);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $gbfile = "";
my $outfile = "";
my $help = 0;

GetOptions ('input=s' => \$gbfile,
            'outputfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $gbdata = parse_genbank($gbfile);

open FH, ">", $outfile or die "couldn't create output file $outfile";
print FH ">" . get_name() . "\n" . get_sequence() . "\n";
close FH;
