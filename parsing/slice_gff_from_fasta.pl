#!/usr/bin/env perl
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use lib "$FindBin::Bin/..";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code subseq_from_fasta);
use lib "$FindBin::Bin";
use GFF qw (feature_to_seq parse_gff_block parse_attributes export_gff_block read_gff_block write_gff_file set_gff_sequence);

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

open my $fh, "<", $gff_file;

foreach my $gene (@sorted_genes) {
	$gff_block = read_gff_block($fh, $gene);

	if ($gff_block eq "") {
		print "No gene named $gene found.\n";
		next;
	}
	my $gff_hash = parse_gff_block ($gff_block);
	my $sequence = $gff_hash->{"sequence"};
	if (! (exists $gff_hash->{"sequence"})) {
		set_gff_sequence ($gff_hash, subseq_from_fasta ($fastafile, $gff_hash->{"start"}, $gff_hash->{"end"}));
	}

	if (@sorted_genes > 1) {
		$outfile = File::Spec->catfile($outdir, "$gene.fasta");
	}

	write_gff_file ($gff_hash, "$outfile.gff");

	open OUT_FH, ">", $outfile or die "couldn't create $outfile";
	print OUT_FH ">$gff_hash->{Name}.gene\n$gff_hash->{sequence}\n";

	my $params = { 'padded' => 0, 'separate' => 1 };
	for (my $mRNA_num = 1; $mRNA_num <= (keys $gff_hash->{"mRNA"}); $mRNA_num++) {
		(undef, my $seq, undef) = split_seq($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"start"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"end"});
		print OUT_FH ">$gff_hash->{Name}.$mRNA_num\n$seq\n";
		my @feature_types = ("five_prime_UTR","exon","three_prime_UTR","CDS");
		foreach my $type (@feature_types) {
			$seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{$type}, $params);
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


__END__

=head1 NAME

select_one_from_fasta

=head1 SYNOPSIS

GetOptions ('gfffile=s' => \$gff_file,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
			'single=s' =>\$gene,
            'help|?' => \$help)

=head1 DESCRIPTION

Finds a single sequence from a fasta file and outputs to a separate file.

=cut
