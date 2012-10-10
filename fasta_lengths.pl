my $file = shift;

open FH, "<", $file;

my $line = readline FH;
my ($name, $length) = 0;
while ($line) {
    if ($line =~ />(.*)\s+/) {
        if ($length > 0) {
            print "$name\t$length\n";
        }
        $name = $1;
        $length = 0;
    } else {
        $length += length ($line) - 1;
	}
	$line = readline FH;
}
print "$name\t$length\n";
close FH;
