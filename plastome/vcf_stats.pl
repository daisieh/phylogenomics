use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($vcf_file, $alt_thresh, $outfile) = 0;
my $keepfiles = 0;
my $help = 0;
my $cov_thresh = 1000;

GetOptions ('samples|input|vcf=s' => \$vcf_file,
            'outputfile:s' => \$outfile,
            'keepfiles!' => \$keepfiles,
            'threshold:i' => \$alt_thresh,
            'minimum|coverage:i' => \$cov_thresh,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

open VCF_FH, "<", $vcf_file or die "couldn't open input file $vcf_file.";
$vcf_file =~ /(.*)\.vcf/;
my $seqname = $1;
unless ($outfile) {
	$outfile = $seqname;
}

# eat header:
my $line = readline VCF_FH;
while ($line =~ m/^#/) {
	$line = readline VCF_FH;
}

my @positions;
my $last_pos = 0;
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	/home/charles/Analysis/ChloroplastGenomes/GRA10_DNA341.bam
#Populus_trichocarpa_chloroplast	1	.	A	.	84	.	DP=18;AF1=7.924e-12;CI95=1.5,0;DP4=18,0,0,0;MQ=59	PL	0
my $i=1;
while ($line !~ m/^\s*$/) {
	my ($chrom,$pos,undef,$ref,$alt,$qual,$filter,$info,undef) = split (/\t/,$line);
	if ($pos != $i) { # there is a gap in the vcf file: add Ns until we get to the next listed position.
		my @this_pos = (0,"N",".",0,0);
		push @positions, \@this_pos;
	} else {
		$info =~ /DP=(.*?);/;
		my $depth = $1;
		$info =~ /.*?DP4=(\d+?),(\d+?),(\d+?),(\d+?);/;
		my $ref_depth = $1 + $2;
		my $alt_depth = $3 + $4;
		if ($pos != $last_pos) {
			my @this_pos = ($depth,$ref,$alt,$ref_depth,$alt_depth);
			push @positions, \@this_pos;
		}
		$last_pos = $pos;
		$line = readline VCF_FH;
	}
	$i++;
}

close VCF_FH;
open FH, ">", $outfile.".stats";

print "$seqname " . @positions . "\n";
print FH "pos\tdepth\tref_freqs\talt_freqs\tref_ratio\talt_ratio\tkeep\n";
for (my $i=0;$i<@positions; $i++) {
	my $pos = $i+1;
	my $ptr = @positions[$i];
	my @this_pos = @$ptr;
	my $depth = @this_pos[0];
	my $ref = @this_pos[1];
	my $alt = @this_pos[2];
	my $ref_depth = @this_pos[3];
	my $alt_depth = @this_pos[4];
	my $keep = 0;
	my $alt_ratio = 0;
	my $ref_ratio = 0;
	if (($alt_depth+$ref_depth) > 0 ) {
		$alt_ratio = $alt_depth/($alt_depth+$ref_depth);
		$ref_ratio = $ref_depth/($alt_depth+$ref_depth);
	}

	if ($alt_ratio > 0.8) {
		$keep = 1;
	} elsif ($ref_ratio > 0.8) {
		$keep = 1;
	}
	print FH "$pos\t$depth\t$ref_depth\t$alt_depth\t$ref_ratio\t$alt_ratio\t$keep\n";
}

close FH;

__END__

=head1 NAME

vcf2fasta

=head1 SYNOPSIS

vcf_stats -vcf_file -output [--keepfiles]

=head1 OPTIONS

  -samples|input|vcf:   vcf file to convert
  -outputfile:      prefix of output fasta file
  -keepfiles:       optional: keep data files that are created (default is no)
=head1 DESCRIPTION

vcf_stats outputs read stats for a vcf file

=cut

