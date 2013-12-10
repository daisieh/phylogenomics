my $genicfile = shift;

open FH, "<", $genicfile;

my $currname = "";
my $currstart = 0;
my $currend = 0;

foreach my $line (<FH>) {
	$line =~ /(.+?)\t(.+?)\t(.+)$/;
	my $name = $1;
	my $start = $2 - 1;
	my $end = $3 + 1;
	if ($currname eq "") {
		$currname = $name;
		$currend = $end;
	} else {
		print "$currname-$name\t$currend\t$start\n";
		$currname = $name;
		$currend = $end;
	}
}
