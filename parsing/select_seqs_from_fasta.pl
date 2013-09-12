use strict;

my $fastafile_1 = shift;
my $seq_blast = shift;

my $seq_names = "$seq_blast.sorted";
my $outfile = "$seq_names.fasta";
system ("cat $seq_blast | sort | uniq -u > $seq_names");

system ("gawk '{if (NF == 0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile_1 | sort | gawk '{ print \">\"\$1\"\\n\"\$2}' FS=\",\" > $fastafile_1.sorted");


open LIST_FH, "<", "$seq_names";
open FA1_FH, "<", "$fastafile_1.sorted";
open OUT_FH, ">", "$outfile";

my $curr_name = readline LIST_FH;

while ($curr_name) {
	chomp $curr_name;
	my $fa_seq1 = (readline FA1_FH) . (readline FA1_FH);

	if ($fa_seq1 eq "") { last; }

	if ($fa_seq1 =~ /$curr_name/) {
		print OUT_FH "$fa_seq1";
		$curr_name = readline LIST_FH;
	}
}

close LIST_FH;
close FA1_FH;
close OUT_FH;
