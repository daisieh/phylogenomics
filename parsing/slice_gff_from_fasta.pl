use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw (make_path);
use lib "$FindBin::Bin/../";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code);
use lib "$FindBin::Bin/";
use GFF qw (feature_to_seq subseq_from_fasta parse_gff_block parse_attributes);

my $gff_file = "";
my $gene = "";
my $fastafile = "";
my $outfile = "";
my $genefile = "";

GetOptions ('gfffile=s' => \$gff_file,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
			'single=s' =>\$gene,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}


if ($outfile eq "") {
	$outfile = "$gene.fasta";
}

my @genes = ();
my $outdir = "";
if ($genefile ne "") {
	if ($gene ne "") {
		pod2usage(-msg => "Either name a single gene or input a file with a list of genes.", -exitval => 2);
	} else {
		$outdir = File::Spec->rel2abs($outfile);
		unless (-d $outdir) {
			make_path($outdir);
		}
		open GENE_FH, "<", $genefile or die "couldn't open gene file.";
		while (my $line = readline GENE_FH) {
			chomp $line;
			push @genes, $line;
		}
	}
} else {
	if ($gene eq "") {
		pod2usage(-msg => "Either name a single gene or input a file with a list of genes.", -exitval => 2);
	}
	push @genes, $gene;
}

@sorted_genes = sort @genes;

open FH, "<", $gff_file;
my $gff_block = "";
my $lastline = "";
foreach my $gene (@sorted_genes) {
	my $in_gene = 0;
	while (my $line = readline FH) {
		# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
		if (($lastline . $line) =~ /gene.*$gene/) {
			if ($lastline =~ /gene.*$gene/) {
				$gff_block = $lastline . $line;
			} else {
				$gff_block = $line;
			}
			$lastline = "";
			$in_gene = 1;
		} elsif (($line =~ /gene/) && ($in_gene == 1)) {
			$lastline = $line;
			last;
		} elsif ($in_gene == 1) {
			$gff_block .= $line;
		}
	}

	if ($gff_block eq "") {
		print "No gene named $gene found.\n";
		exit;
	}
	my $gff_hash = parse_gff_block ($gff_block);
	$gff_hash->{"sequence"} = subseq_from_fasta ($fastafile, $gff_hash->{"start"}, $gff_hash->{"end"});

	if (@sorted_genes > 1) {
		$outfile = File::Spec->catfile($outdir, "$gene.fasta");
	}

	open OUT_FH, ">", $outfile or die "couldn't create $outfile";
	print OUT_FH ">$gff_hash->{Name}.gene\n$gff_hash->{sequence}\n";

	my $params = { 'padded' => 0, 'separate' => 1 };
	for (my $mRNA_num = 1; $mRNA_num <= (keys $gff_hash->{"mRNA"}); $mRNA_num++) {
		print "mRNA $mRNA_num:\n";
		my @feature_types = ("five_prime_UTR","exon","three_prime_UTR","CDS");
		foreach my $type (@feature_types) {
			my $seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{$type}, $params);
			print @$seq . " $type\n";
			if ((ref $seq) =~ /ARRAY/ ) {
				for (my $i=1; $i<=@$seq; $i++) {
					print OUT_FH ">$gff_hash->{Name}.$mRNA_num.$type.$i\n@$seq[$i-1]\n";
				}
			}
		}
	}
	close OUT_FH;
}

close FH;


