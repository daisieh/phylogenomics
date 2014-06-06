use File::Spec qw (catfile);

my $genelist = shift;
my $genedir = shift;
my $cutoff = shift;

if ($cutoff > 0) {
	$cutoff = $cutoff / 100;
}
my @genes = ();
open FH, "<", $genelist;
foreach my $line (<FH>) {
	chomp $line;
	push @genes, $line;
}
close FH;

my $perc = "0";
$cutoffstr = sprintf("%.2f", $cutoff);
print @genes . " genes\n";
if ($cutoffstr =~ /\d+\.(\d{2})/) {
	$perc = $1;
}
open OUTFH, ">", "cutoff.$perc";
foreach my $gene (@genes) {
	open FH, "<", File::Spec->catfile($genedir, $gene);
	my $totalident = 0;
	my $totallength = 0;
	foreach my $line (<FH>) {
#Potri.001G230000.1.CDS.12	5876	5995	118	123
		my ($name, $start, $end, $ident, $length) = split (/\t/, $line);
		if ($name =~ /\.1\.CDS\./) {
			$totalident += $ident;
			$totallength += $length;
		}
	}
	close FH;
	if (($totallength > 0) && (($totalident/$totallength) > $cutoff)){
		print OUTFH "$gene\t$totalident / $totallength = ".($totalident/$totallength)."\n";
		print "$gene\t".($totalident/$totallength)."\n";
	}
}
close OUTFH;
