use File::Basename;

my $usage  = "perl " . basename($0) . " file1 file2...\n";
$usage .= qq{
Takes any number of input lists of the same length and creates a tab-delimited
file with each list as a column.
};

use File::Basename;
use Getopt::Long;

my @files = ();
my $has_names = 0;
my $out_file = 0;

GetOptions ('files|input=s{2,}' => \@files,
            'names!' => \$has_names,
            'outputfile:s' => \$out_file) or die "options misspecified";

my $result = combine_files (\@files, $has_names, $out_file);
print $result;

sub combine_files {
    my $fileptr = shift;
    my $has_names = shift;
    my $out_file = shift;

    my @files = @$fileptr;

    if (@files < 2) { die $usage; }

    my @inputs;
    for (my $i=0; $i<@files; $i++) {
        open FH, "<", @files[$i] or die $usage;
        my @data = <FH>;
        close FH;
        push @inputs, \@data;
    }

    my $num_entries = scalar @{@inputs[0]};
    for (my $i = 0; $i < @files; $i++) {
        if (scalar @{@inputs[$i]} != $num_entries) {
            die "Error: files have different numbers of inputs." . scalar @{@inputs[$i]};
        }
    }
    my $result = "";

    for (my $j = 0; $j < $num_entries; $j++) {
        if ($has_names) {
            my $entry = @{@inputs[0]}[$j];
            $entry =~ /(.+?)\t(.*)/;
            $result .= "$1\t";
        }
        for (my $i = 0; $i < @files; $i++) {
            my $entry = @{@inputs[$i]}[$j];
            $entry =~ s/\n//g;
            if ($has_names) {
                $entry =~ /(.+?)\t(.*)/;
                $entry = $2;
            }
            $result .= $entry . "\t";
        }
        $result .= "\n";
    }

    if ($out_file) {
        open FH, ">", $out_file;
        print FH $result;
        return "";
    } else {
        return $result;
    }
}
__END__

=head1 NAME

combine_files

=head1 SYNOPSIS

combine_files [options]


=head1 OPTIONS

    -files|input:   list of files to combine
	-outputfile:    name of output file
	-(no)names:     files have the same column of names (first column of each file)

=head1 DESCRIPTION

Takes any number of input lists of the same length and creates a tab-delimited
file with each list as a column.

=cut
