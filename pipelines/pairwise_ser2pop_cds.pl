#!/usr/bin/env perl
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code subseq_from_fasta parse_fasta);
use lib "$FindBin::Bin/../parsing";
use GFF qw (feature_to_seq parse_gff_block parse_attributes export_gff_block read_gff_block write_gff_file set_gff_sequence);

my $gff_file = "";
my $refdir = "";
my $outdir = "";

GetOptions ('gfffile=s' => \$gff_file,
			'refdir=s' => \$refdir,
			'output=s' => \$outdir,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $outdir = File::Spec->rel2abs($outdir);
unless (-d $outdir) {
	make_path($outdir);
}

my $pop_name = "";
my $sal_name = "";
if ($gff_file =~ /(Ser\.ph\.)(.+)\.gff/) {
	$pop_name = "$2";
	$sal_name = "$1$2";
}

my $fastafile = File::Spec->catfile($refdir, "$pop_name.fasta");

my ($pophash, $poparr) = parse_fasta ($fastafile);
my @genes = ();
my $sequence = "";
open my $fh, "<", $gff_file or die "couldn't open $gff_file";
my $in_gff = 1;
my $in_fasta = 0;
# my $line;
while (my $line = readline $fh) {
	if ($line =~ /##gff/) {
		$in_gff = 1;
		next;
	}
	if ($line =~ /##FASTA/) {
		$in_gff = 0;
		$in_fasta = 1;
		$line = readline $fh;
		next;
	}
	if ($in_gff == 1) {
	#Chr01	phytozome9_0	gene	1	819	.	+	.	ID=Potri.001G001000;Name=Potri.001G001000;
		seek($fh, -length($line), 1);
		my $gff_block = read_gff_block ($fh);
		push @genes, $gff_block;
		next;
	}
	if ($in_fasta == 1) {
		chomp $line;
		$sequence .= $line;
		next;
	}
}

@sorted_genes = sort @genes;

foreach my $gff_block (@genes) {
	my $gff_hash = parse_gff_block ($gff_block);
	set_gff_sequence ($gff_hash, $sequence);
	my $gene = $gff_hash->{"ID"};

	$outfile = File::Spec->catfile($outdir, "$gene.fasta");

	open OUT_FH, ">", $outfile or die "couldn't create $outfile";

	my $params = { 'padded' => 0, 'separate' => 1 };
	my $mRNA_num = 1;
	(undef, my $seq, undef) = split_seq($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"start"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"end"});
	$seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"CDS"}, $params);
	if ((ref $seq) =~ /ARRAY/ ) {
		print OUT_FH ">$sal_name.$mRNA_num.CDS\n";
		for (my $i=1; $i<=@$seq; $i++) {
			print OUT_FH "@$seq[$i-1]";
		}
		print OUT_FH "\n>$pop_name.$mRNA_num.CDS\n";
		for (my $i=1; $i<=@$seq; $i++) {
			my $cds = "$pop_name.$mRNA_num.CDS.$i";
			print OUT_FH "$pophash->{$cds}";
		}
		print OUT_FH "\n";
	}

	close OUT_FH;
}

close FH;


__END__

=head1 NAME

select_one_from_fasta

=head1 SYNOPSIS

GetOptions ('gfffile=s' => \$gff_file,
			'refdir=s' => \$refdir,
			'outdir=s' => \$outdir,
            'help|?' => \$help)

=head1 DESCRIPTION

Finds a single sequence from a fasta file and outputs to a separate file.

=cut
