use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $samplefile = 0;
my $help = 0;
my $cov_thresh = 5;
my $pl_thresh = 0;
my $outfile = "";

GetOptions ('samples|input|vcf=s' => \$samplefile,
            'outputfile=s' => \$outfile,
            'threshold|pl=i' => \$pl_thresh,
            'minimum|coverage|reads=i' => \$cov_thresh,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my @samples = @{sample_list ($samplefile)};
if ($samplefile =~ /(.*?)\.vcf$/) {
	@samples = ($1);
}

my $result_str = "";
my $name;

my @AMBIGUOUS_POS = (0,"N",".",0,0,0);

foreach my $sample (@samples) {
	$name = basename($sample);
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
		my @this_pos = ();
		my ($chrom,$pos,undef,$ref,$alt,$qual,$filter,$info,undef,$pl) = split (/\t/,$line);
		if ($pos != $i) { # there is a gap in the vcf file: add Ns until we get to the next listed position.
			push @this_pos, @AMBIGUOUS_POS;
			push @positions, \@this_pos;
		} else {
			if ($line =~ /INDEL/) { # we don't handle indels; drop these lines.
				push @this_pos, @AMBIGUOUS_POS;
				push @positions, \@this_pos;
			} else {
				$info =~ /DP=(.*?);/;
				my $depth = $1;
				$info =~ /.*?DP4=(\d+?),(\d+?),(\d+?),(\d+?);/;
				my $ref_depth = $1 + $2;
				my $alt_depth = $3 + $4;
				chomp $pl;
				if ($pos != $last_pos) {
					push @this_pos, ($depth,$ref,$alt,$ref_depth,$alt_depth,$pl);
					push @positions, \@this_pos;
				}
				$last_pos = $pos;
			}
			$line = readline VCF_FH;
		}
		$i++;
	}

	close VCF_FH;

	my $seq = "";
	print "processed sample $name with length " . @positions . "\n";
	for (my $i=0;$i<@positions; $i++) {
		my $ptr = @positions[$i];
		my @this_pos = @$ptr;
		my $pos = $i+1;
		my $depth = @this_pos[0];
		my $ref = @this_pos[1];
		my $alt = @this_pos[2];
		my $ref_depth = @this_pos[3];
		my $alt_depth = @this_pos[4];
		my $pl = @this_pos[5];

		if ($depth > $cov_thresh) {
			if ($pl =~ /,/) {
				my @pls = split(/,/, $pl);
				my $alleles = get_allele_str([$ref, $alt]);
				my @genotypes = @{get_ordered_genotypes($alleles)};
				$alt = "N";
				foreach (my $g=0;$g<@pls;$g++) {
					if ($pls[$g] <= $pl_thresh) {
						$alt = get_iupac_code($genotypes[$g]);
					}
				}
				$seq .= $alt;
			} elsif ($pl <= $pl_thresh) {
				$seq .= $ref;
			} else {
				$seq .= "N";
			}
		} else {
			$seq .= "n";
		}
	}

	$result_str .= ">$name\n$seq\n";
}

if (@samples == 1) {
	unless ($outfile ne "") {
		$outfile = $name;
	}
}


open FH, ">", $outfile.".fasta";
print FH "$result_str";
close FH;



__END__

=head1 NAME

vcf2fasta

=head1 SYNOPSIS

vcf2fasta -samplefile -output [-threshold]

=head1 OPTIONS

  -samples|input|vcf:   name of sample or list of samples to convert
  -outputfile:      optional: prefix of output fasta file
  -threshold|pl:	optional: maximum threshold for Phred-scaled likelihood (lower is best, default is 0)
  -min|coverage:    optional: read depth required to call base (default is 5)

=head1 DESCRIPTION

Graphs coverage depths using depth data generated from samtools mpileup.

=cut

