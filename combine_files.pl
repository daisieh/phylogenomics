use File::Basename;

my $usage  = "perl " . basename($0) . " file1 file2...\n";
$usage .= qq{
Takes any number of input lists of the same length and creates a tab-delimited
file with each list as a column.
};

my $num_files = scalar @ARGV;
if ($num_files < 1) { die $usage; }

my @inputs;
for (my $i = 0; $i < $num_files; $i++) {
    open FH, "<", @ARGV[$i] or die $usage;
    my @data = <FH>;
    close FH;
    push @inputs, \@data;
}

my $result = "";
my $num_entries = scalar @{@inputs[0]};
print "$num_entries\n";
for (my $i = 0; $i < $num_files; $i++) {
    if (scalar @{@inputs[$i]} != $num_entries) {
        die "Error: files have different numbers of inputs." . scalar @{@inputs[$i]};
    }
}

for (my $j = 0; $j < $num_entries; $j++) {
    for (my $i = 0; $i < $num_files; $i++) {
        my $entry = @{@inputs[$i]}[$j];
        $entry =~ s/\n//g;
        $result .= $entry . "\t";
    }
    $result .= "\n";
}
    print "$result\n";
