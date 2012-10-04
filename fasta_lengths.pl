my $file = shift;

open FH, "<", $file;

my $line = readline FH;
while ($line) {
	$line =~ />(.*)\s+/;
	my $name = $1;
	$line = readline FH;
	my $length = length ($line) - 1;
	$line = readline FH;
	print "$name\t$length\n";
}

close FH;