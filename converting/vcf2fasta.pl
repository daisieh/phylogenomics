use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $samplefile = 0;
my $keepfiles = 0;
my $help = 0;
my $cov_thresh = 10;
my $alt_thresh = 0.8;
my $keep_ambig = 0;
my $outfile = "result";

GetOptions ('samples|input=s' => \$samplefile,
            'outputfile:s' => \$outfile,
            'keepfiles!' => \$keepfiles,
            'threshold:i' => \$alt_thresh,
            'minimum|coverage:i' => \$cov_thresh,
            'iupac|ambiguity' => \$keep_ambig,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my @samples = @{sample_list ($samplefile)};
if ($samplefile =~ /(.*?)\.vcf$/) {
	@samples = ($1);
}

if (@samples == 1) {
	unless ($outfile) {
		$outfile = $samples[0];
	}
}
my $result_str = "";

foreach my $sample (@samples) {
	my $vcf_file = $sample . ".vcf";
	open VCF_FH, "<", $vcf_file or die "couldn't open input file $vcf_file.";

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

	my $seq = "";
	print "processed sample $sample with length " . @positions . "\n";
	for (my $i=0;$i<@positions; $i++) {
		my $ptr = @positions[$i];
		my @this_pos = @$ptr;
		my $depth = @this_pos[0];
		my $ref = @this_pos[1];
		my $alt = @this_pos[2];
		my $ref_depth = @this_pos[3];
		my $alt_depth = @this_pos[4];
		if ($depth > $cov_thresh) {
			if ($alt_depth/($alt_depth+$ref_depth) > $alt_thresh) {
				if ($keep_ambig) {
					$alt = get_iupac_code([$ref, $alt]);
				}
				$seq .= $alt;
			} elsif ($ref_depth/($alt_depth+$ref_depth) > $alt_thresh) {
				$seq .= $ref;
			} else {
				$seq .= "-";
			}
		} else {
			$seq .= "-";
		}
	}

	$result_str .= ">$sample\n$seq\n";
}

open FH, ">", $outfile.".fasta";
print FH "$result_str";
close FH;

__END__

=head1 NAME

vcf2fasta

=head1 SYNOPSIS

vcf2fasta -samplefile -output [-threshold] [--keepfiles]

=head1 OPTIONS

  -samples|input|vcf:   name of sample or list of samples to convert
  -outputfile:      optional: prefix of output fasta file
  -keepfiles:       optional: keep data files that are created (default is no)
  -threshold:		optional: threshold ratio to accept alternate allele (default is 0.8)
  -min|coverage:    optional: coverage required to call base (default is 1000)

=head1 DESCRIPTION

Graphs coverage depths using depth data generated from samtools mpileup.

=cut

