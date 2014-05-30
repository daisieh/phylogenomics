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
my @rowmasks = ();

# first pass: remove any rows with excessive missing data
my $deleted_rows = 0;
foreach my $seqid (@$seqids) {
	my $ambigs = ($seqmatrix->{$seqid} =~ tr/Nn\-\?//);
	if ($ambigs/length ($seqmatrix->{$seqid}) < $row_thresh) {
		push @rows, $seqmatrix->{$seqid};
		push @rowids, $seqid;
		my $rowmask = "$seqmatrix->{$seqid}";
		$rowmask =~ tr/Nn\-\?/1/;
		$rowmask =~ tr/1/0/c;
		push @rowmasks, $rowmask;
	} else {
		$deleted_rows++;
	}
}

if (@rows == 0) {
	die "EMPTY MATRIX: no sequences had less missing data than the specified threshold\n";
}

# second pass: remove columns with excessive missing data
my @col_counts = 0;
for (my $i=0;$i<length($rows[0]);$i++) {
	my $col_count = 0;
	foreach my $row (@rowmasks) {
		$row =~ /(.)(.*)/;
		$col_count += $1;
		$row = $2;
	}
	$col_counts[$i] = $col_count;
}

for (my $j=0; $j< @col_counts; $j++) {
	if (($col_counts[$j]/@rows) > $col_thresh) {
		foreach my $row (@rows) {
			$row =~ /(.{$j})(.)(.*)/;
			$row = "$1#$3";
		}
		if (length ($rows[0]) == 0) {
			die "EMPTY MATRIX: no base positions had less missing data than the specified threshold\n";
		}
	}
}

my $deleted_cols = ($rows[0] =~ tr/#//);
print "deleted $deleted_rows rows, $deleted_cols cols\n";
open FH, ">", "$outname.fasta";
for (my $i=0;$i<@rows;$i++) {
	my $row = $rows[$i];
	$row =~ s/#//g;
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

