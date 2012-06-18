use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Yn00;
use Bio::Tools::Run::Phylo::PAML::Codeml;

my $aln = Bio::SimpleAlign->new();
my @seqs = (
		["Pdel181_DNA2", "ATGGATATAGTAAGTCTCGCTTCGGCTGCTTTGATGGTAGTCTTTACATTTTCCCTTTCACTCGTAGTATGGGGAAGAAGTGGACTC"],
		["Pfre186_DNA2", "ATGGATATAGTAAGTCTCGCTTCGGCTGCTTTGATGGTAGTCTTTACATTTTCCCTTTCACTCGTAGTATGGGGAAGAAGTGGACTC"],
		["Pgra187_DNA2", "ATGGATATAGTAAGTCTCGCTTGGGCTGCTTTGATGGTAGTCTTTACATTTTCCCTTTCACTCGTAGTATGGGGAAGAAGTGGACTC"],
		["Phet26_DNA21", "ATGGATATAGTAAGTCTCGCTTGGGCTGCTTTGATGGTAGTCTTTACATTTTCCCTTTCACTCGTAGTATGGGGAAGAAGTGGACTC"]
		);

for (my $i=0; $i<= $#seqs; $i++) {
	my $outseq = Bio::LocatableSeq->new(-seq => $seqs[$i][1], -id => $seqs[$i][0]);
	$aln->add_seq ($outseq);
}

my $paml_exec = Bio::Tools::Run::Phylo::PAML::Codeml->new
			   ( -params => { 'runmode' => -2,
							  'seqtype' => 1,
							} );
# my $paml_exec = Bio::Tools::Run::Phylo::PAML::Yn00->new();
print "There were " . scalar @seqs . " seqs input into the analysis.\n";
$paml_exec->alignment($aln);
my ($rc,$parser) = $paml_exec->run();
if ($rc == 0) {
	my $t = $paml_exec->error_string();
	print "problem in $name: $t\n";
}

while( my $result = $parser->next_result ) {
	my @otus = $result->get_seqs();
	print "There were " . scalar @otus . " seqs output from the analysis.\n";
	print "seq1\tseq2\tdN/dS\tdN\tdS\n";

	my $MLmatrix = $result->get_MLmatrix();
	my $j = 1;
	for (my $i=1;$i<=$j+1;$i++) {
		if ($i == scalar @otus) { last; }
		for ($j=0;$j<$i;$j++) {
			my $dN = $MLmatrix->[$j]->[$i]->{dN};
			my $dS = $MLmatrix->[$j]->[$i]->{dS};
			my $kaks =$MLmatrix->[$j]->[$i]->{omega};
			my $seq1 = @otus[$i]->display_name();
			my $seq2 = @otus[$j]->display_name();
			print "$seq1\t$seq2\t$kaks\t$dN\t$dS\n";
		}
	}
}
