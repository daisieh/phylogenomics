require "bioperl_subfuncs.pl";
use Bio::Tools::Run::Phylo::PAML::Codeml;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($gb_file, $fa_file, $output_name, $keepfiles) = 0;
GetOptions ('genbank|gb=s' => \$gb_file,
            'fasta=s' => \$fa_file,
			'keepfiles' => \$keepfiles,
            'output=s' => \$output_name) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns = ();
if ($gb_file) {
	@gene_alns = @{parse_aln_into_genes($whole_aln, $gb_file, 1)};
} else {
	$whole_aln->description("gene");
	push @gene_alns, $whole_aln;
}

open my $total_file, ">$output_name.total";
truncate $total_file, 0;

print $total_file "running codeml on $gb_file and $fa_file\n";

open my $formatted_file, ">$output_name.formatted";
truncate $formatted_file, 0;


my $paml_exec = Bio::Tools::Run::Phylo::PAML::Codeml->new
			   ( -params => { 'runmode' => -2, 'seqtype' => 1, 'model' => 1} );

#$paml_exec->save_tempfiles(1);
if ($keepfiles) {
	$paml_exec->save_tempfiles(1);
	print "temporary files are located in " . $paml_exec->tempdir() . "\n";
}

print $total_file scalar @gene_alns . " genes in the analysis\n";
print $total_file $whole_aln->num_sequences() . " sequences in the analysis\n\n";
print $total_file "gene\tseq1\tseq2\tdN/dS\tdN\tdS\n";


foreach my $aln (@gene_alns) {
	my $name = $aln->description();
	print "$name\n";
	my $resultstr = $name;
	$paml_exec->alignment($aln);
	$paml_exec->outfile_name("$output_name"."_$name.mlc");
	my ($rc,$parser) = $paml_exec->run();
	if ($rc == 0) {
		my $t = $paml_exec->error_string();
  		print $total_file "problem in $name: $t\n";
	} else {
		while( my $result = $parser->next_result ) {
			my @otus = $result->get_seqs();
			my $MLmatrix = $result->get_MLmatrix();
			#	this loop set is excessively complicated because I am trying to get the output to correspond to yn00's output block.
			my $j = 1;
			for (my $i=1;$i<=$j+1;$i++) {
				if ($i == scalar @otus) { last; }
				for ($j=0;$j<$i;$j++) {
					my $dN = $MLmatrix->[$j]->[$i]->{dN};
					my $dS = $MLmatrix->[$j]->[$i]->{dS};
					my $kaks =$MLmatrix->[$j]->[$i]->{omega};
					my $seq1 = @otus[$i]->display_name();
					my $seq2 = @otus[$j]->display_name();
					$resultstr .= "\t$seq1#$seq2#" . $kaks;
					print $total_file "$name\t$seq1\t$seq2\t$kaks\t$dN\t$dS\n";
				}
			}
		}
	}
	print $formatted_file $resultstr, "\n";
}


__END__

=head1 NAME

pairwise_omega

=head1 SYNOPSIS

pairwise_omega -fasta fastafile [-genbank genbankfile] -output outfile_prefix

Performs pairwise codeml on the sequences given in the fasta file. If an optional genbank
file is specified, will parse the fasta alignment into separate genes for analysis.

=head1 OPTIONS

    -fasta          fasta file of whole aligned sequences
    -genbank|gb     (optional) genbank file listing the genes to analyze
    -output         output file name

=head1 DESCRIPTION

=cut
