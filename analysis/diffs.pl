use strict;
use File::Temp qw/ tempfile /;

my $fastafile = shift;

my (undef, $tempfile) = tempfile(UNLINK => 1);
flattenfasta ($fastafile, $tempfile, "\t");

open FH, "<", $tempfile;
my @seqs = ();
my @names = ();
foreach my $line (<FH>) {
	$line =~ />(.*?)\t(.*)$/;
	push @names, $1;
	push @seqs, $2;
}
close FH;

my $total = 0;
my $gaps = 0;
my $diffs = 0;
my $length = length ($seqs[0]);

while ($seqs[0]) {
	my @bases = ();
	for (my $i=0; $i<@seqs; $i++) {
		if ($seqs[$i] =~ /^(.)(.*)$/) {
			push @bases, "$1";
			$seqs[$i] = "$2";
		}
	}
# 	print "$total\n";
	if (($bases[0] eq "-") || ($bases[1] eq "-")) {
		$gaps++;
		next;
	} else {
		$total++;
	}
	if ($bases[0] ne $bases[1]) {
		$diffs++;
	}
}

print "file $fastafile:\n";
print "total length: $length\n";
print "total gaps: $gaps\n";
print "total compared: $total\n";
print "total diffs: $diffs\n";
print "percent divergence: " . (($diffs / $total) * 100) . "\n";

sub flattenfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my (undef, $tempfile) = tempfile(UNLINK => 1);
	system ("gawk '{if (NF==0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | gawk '{print \">\" \$1 \"$separator\" \$2}' FS=\",\" > $outfile");
}
