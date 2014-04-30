use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw (make_path);
use lib "$FindBin::Bin/../";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code parse_fasta);
use lib "$FindBin::Bin/";
use GFF qw (feature_to_seq subseq_from_fasta parse_gff_block parse_attributes export_gff_block read_gff_block);

my $gff_file = "";
my $gene = "";
my $fastafile = "";
my $outfile = "";
my $genefile = "";
my $blastfile = "";

GetOptions ('gfffile=s' => \$gff_file,
			'blastdir=s' => \$blastfile,
			'fastadir=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my @genes = ();
my $outdir = File::Spec->rel2abs($outfile);
unless (-d $outdir) {
	make_path($outdir);
}
open GENE_FH, "<", $genefile or die "couldn't open gene file.";
while (my $line = readline GENE_FH) {
	chomp $line;
	push @genes, $line;
}

@sorted_genes = sort @genes;

open my $fh, "<", $gff_file;

foreach my $gene (@sorted_genes) {
	$gff_block = read_gff_block($fh, $gene);

	if ($gff_block eq "") {
		print "No gene named $gene found.\n";
		exit;
	}

	my $gff_hash = parse_gff_block ($gff_block);
# 	print Dumper($gff_hash);
	(my $seqhash, undef) = parse_fasta (File::Spec->catfile ($fastafile, "$gene.fasta"));
	$gff_hash->{"sequence"} = $seqhash->{$gene};

	$gff_hash->{"start"} = 1;
	$gff_hash->{"end"} = length ($gff_hash->{"sequence"});
	$gff_hash->{"seqid"} = $gene;
	$gff_hash->{"source"} = "Ser_aTRAM";

	open BLASTFH, "<", File::Spec->catfile ($blastfile, $gene);
	foreach my $line (<BLASTFH>) {
		if ($line =~ /$gene\.gene/) {
			next;
		}

# 		Potri.001G000200.1.exon.1	525	608
		my ($longname,$start,$end) = split (/\t/, $line);
		chomp $end;
		if ($longname =~ /$gene\.(\d+)\.(.+?)\.(\d+)/) {
			my $mRNA = $1;
			my $type = $2;
			my $num = $3;
# 			print "replacing $longname end " . $gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"end"};
			$gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"start"} = $start;
			$gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"end"} = $end;
# 			print " with " . $gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"end"} . "\n";
		}
	}

	open OUTFH, ">", File::Spec->catfile ($outdir, "$gene.gff");
	print OUTFH "##gff-version 3\n";
	print OUTFH export_gff_block ($gff_hash);
	print OUTFH "##FASTA\n";
	print OUTFH ">$gene\n$seqhash->{$gene}\n";
	close OUTFH;
}

close $fh;


