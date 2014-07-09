use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $inputfile = "";
my $outname = "trimmed";
my ($help, $row_thresh, $col_thresh) = 0;

GetOptions ('fasta|input=s' => \$inputfile,
            'outputfile=s' => \$outname,
            'rowthreshold=f' => \$row_thresh,
            'colthreshold=f' => \$col_thresh,
            'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my ($seqmatrix, $seqids)  = parse_fasta ( $inputfile );

my @rows = ();
my @rowids = ();
my @final_rows = ();

# first pass: remove any rows with excessive missing data
my $deleted_rows = 0;
foreach my $seqid (@$seqids) {
	my $ambigs = ($seqmatrix->{$seqid} =~ tr/Nn\-\?//);
	my $missing_frac = $ambigs/length ($seqmatrix->{$seqid});
	if ($missing_frac < $row_thresh) {
		push @rows, $seqmatrix->{$seqid};
		push @rowids, $seqid;
		push @final_rows, "";
	} else {
		$deleted_rows++;
	}
	print "$seqid has $missing_frac missing data\n";
}

if (@rows == 0) {
	die "EMPTY MATRIX: no sequences had less missing data than the specified threshold\n";
}
my $total_length = length($rows[0]);
print "there are " . @final_rows . " final rows\n";
# second pass: remove columns with excessive missing data
# if more than $col_thresh * @final_rows Ns in a column, delete and move on.
my $max_col_ambigs = $col_thresh * @final_rows;

my $deleted_cols = 0;
for (my $i=0;$i<$total_length;$i++) {
	my $column = "";
	foreach my $row (@rows) {
		$row =~ /(.)(.*)/;
		$row = $2;
		$column .= $1;
	}
	my $ambigs = ($column =~ tr/Nn\-\?//);

	if ($ambigs < $max_col_ambigs) {
		foreach my $row (@final_rows) {
			$column =~ /(.)(.*)/;
			$row .= $1;
			$column = $2;
		}
	} else {
		$deleted_cols ++;
		print "deleting col $i\n";
	}
}

print "deleted $deleted_rows rows, $deleted_cols cols\n";
open FH, ">", "$outname.fasta";
for (my $i=0;$i<@rows;$i++) {
	my $row = $final_rows[$i];
	print FH ">$rowids[$i]\n$row\n";
}
close FH;


__END__

=head1 NAME

trim_missing

=head1 SYNOPSIS

trim_missing -input fastafile -output outputfile -row row_missing_frac -col col_missing_frac

=head1 OPTIONS

GetOptions ('fasta|input=s' => \$inputfile,
            'outputfile=s' => \$outname,
            'rowthreshold=f' => \$row_thresh,
            'colthreshold=f' => \$col_thresh,
            'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

  -fasta|input:     fasta file of aligned sequences.
  -outputfile:      output file name.
  -row:	            fraction of missing data allowed per row
  -col:	            fraction of missing data allowed per column

=head1 DESCRIPTION

Takes a fasta file and removes rows and columns that contain more missing data than the specified thresholds.

=cut

