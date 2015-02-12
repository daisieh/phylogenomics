#!/usr/bin/env perl

use File::Basename;
use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(combine_files);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my @files = ();
my ($has_names, $has_header, $out_file) = 0;

GetOptions ('files|input=s{1,}' => \@files,
            'names!' => \$has_names,
            'header!' => \$has_header,
            'outputfile:s' => \$out_file) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);


if (@files == 1) {
	my $file = pop @files;
	if (-e $file) {
		open FH, "<:crlf", $file;
		foreach my $line (<FH>) {
			chomp $line;
			push @files, $line;
		}
	}
}

my $result = combine_files (\@files, $has_names, $has_header, $out_file);

if ($out_file) {
    open FH, ">", $out_file;
    print FH $result;
} else {
    print $result;
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
	-(no)header:    files have column labels in first row

=head1 DESCRIPTION

Takes any number of input files containing tab-delimited lists of the same length
and creates a tab-delimited file with each list as a column. If outputfile is not
specified, prints the value on STDOUT.

Default is no column names, no common row names.

=cut
