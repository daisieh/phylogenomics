use File::Basename;
use Getopt::Long;

my @files = ();
my $has_names = 0;
my $has_header = 0;
my $out_file = 0;

GetOptions ('files|input=s{2,}' => \@files,
            'names!' => \$has_names,
            'header!' => \$has_header,
            'outputfile:s' => \$out_file) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

my $result = combine_files (\@files, $has_names, $has_header, $out_file);
print $result;

sub combine_files {
    my $fileptr = shift;
    my $has_names = shift;
    my $has_header = shift;
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

    my $result = "";
    my @labels = ();
    my $num_entries = scalar @{@inputs[0]};
    for (my $i = 0; $i < @files; $i++) {
        if (scalar @{@inputs[$i]} != $num_entries) {
            die "Error: files have different numbers of inputs." . scalar @{@inputs[$i]};
        }
        if ($has_header) {
            my @heads = split /\t/, (@{@inputs[$i]}[0]);
            foreach $head (@heads) {
                $head = @files[$i] . "|" . $head;
            }
            @{@inputs[$i]}[0] = join ("\t", @heads);
        }
    }

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
        if ($has_header && $has_names && ($j==0)) {
            #clean up the header row
            $result =~ s/^.*?\|//;
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
