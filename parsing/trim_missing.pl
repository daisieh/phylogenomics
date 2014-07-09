use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta pad_seq_ends debug set_debug);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $inputfile = "";
my $outname = "trimmed";
my ($help, $row_frac_missing, $col_frac_missing) = 0;
my $debug = 0;

GetOptions ('fasta|input=s' => \$inputfile,
            'outputfile=s' => \$outname,
            'rowmax=f' => \$row_frac_missing,
            'colmax=f' => \$col_frac_missing,
            'debug' => \$debug,
            'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;
set_debug ($debug);

my ($seqmatrix, $seqids)  = parse_fasta ( $inputfile );

my @rows = ();
my @rowids = ();

# first pass: remove any rows with excessive missing data
my $deleted_rows = 0;
foreach my $seqid (@$seqids) {
	my $ambigs = ($seqmatrix->{$seqid} =~ tr/Nn\-\?//);
	my $missing_frac = $ambigs/length ($seqmatrix->{$seqid});
	if ($missing_frac < $row_frac_missing) {
		push @rows, $seqmatrix->{$seqid};
		push @rowids, $seqid;
	} else {
		$deleted_rows++;
	}
}

if (@rows == 0) {
	die "EMPTY MATRIX: no sequences had less missing data than the specified threshold\n";
}

# second pass: remove columns with excessive missing data
# if more than $col_frac_missing * @rows Ns in a column, delete and move on.

pad_seq_ends (\@rows);
my $total_length = length($rows[0]);
print ("started with " . @$seqids . " rows, $total_length columns\n");

my $max_col_ambigs = int($col_frac_missing * @rows);

my @final_rows = @{check_block(\@rows, $max_col_ambigs)};

my $deleted_cols = $total_length - length (@final_rows[0]);

print ("deleted $deleted_rows rows, $deleted_cols columns\n");
open FH, ">", "$outname.fasta";
for (my $i=0;$i<@rows;$i++) {
	my $row = $final_rows[$i];
	print FH ">$rowids[$i]\n$row\n";
}
close FH;

# recursive function:
# if the number of Ns is less than $max_col_ambigs, this group is fine, return the array.
# if the matrix is a single column with more than $max_col_ambigs, return 0.
# else divide the matrix into halves, recurse on both halves, merge results.

sub check_block {
	my $seq_array = shift;
	my $max_col_ambigs = shift;

	# if we didn't get an array, return 0. (quick exit)
	if (($seq_array == 0) || ($seq_array == ())) {
		debug ("\n");
		return 0;
	}

	my $seqlen = length (@$seq_array[0]);

	# if we did get an array, but the length of the strings is 0, then return 0. (quick exit)
	if ($seqlen == 0) {
		debug ("empty\n");
		return 0;
	}

	# count the number of Ns in this block.
	my $block_missing = 0;
	foreach my $row (@$seq_array) {
		$block_missing += ($row =~ tr/Nn\-\?//);
	}

	debug ("block of $seqlen has $block_missing missing (max $max_col_ambigs): ");

	if ($block_missing <= $max_col_ambigs) {
		# the block is fine. No deletions.
		debug ("KEEP BLOCK\n");
		return $seq_array;
	} elsif ($seqlen == 1) {
		# base case: if we're down to a single column and it still has more than $max_col_ambigs, delete it.
		debug ("DELETE a col\n");
		return 0;
	} else {
		# split this block into two parts: recurse on the first part, then the second part
		# then merge those two finished blocks.
		my $numcols = int ($seqlen / 2);

		# perl regex limit:
		if ($numcols > 32766) {
			debug (" (regex max) ");
			$numcols = 32766;
		}

		debug ("TWO BLOCKS of $numcols, ".($seqlen-$numcols)."\n");
		my $front_block = ();
		foreach my $row (@$seq_array) {
			$row =~ /(.{$numcols})(.*)/;
			push @$front_block, $1;
			$row = $2;
		}
		$front_block = check_block ($front_block, $max_col_ambigs);
		$seq_array = check_block ($seq_array, $max_col_ambigs);

		if ($front_block == 0) {
			return $seq_array;
		}
		if ($seq_array == 0) {
			return $front_block;
		}

		for (my $i=0; $i< @$seq_array; $i++) {
			@$seq_array[$i] = @$front_block[$i] . @$seq_array[$i];
		}
	}
	return $seq_array;
}

__END__

=head1 NAME

trim_missing

=head1 SYNOPSIS

trim_missing -input fastafile -output outputfile -row row_missing_frac -col col_missing_frac

=head1 OPTIONS

GetOptions ('fasta|input=s' => \$inputfile,
            'outputfile=s' => \$outname,
            'rowthreshold=f' => \$row_frac_missing,
            'colthreshold=f' => \$col_frac_missing,
            'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

  -fasta|input:     fasta file of aligned sequences.
  -outputfile:      output file name.
  -row:	            fraction of missing data allowed per row
  -col:	            fraction of missing data allowed per column

=head1 DESCRIPTION

Takes a fasta file and removes rows and columns that contain more missing data than the specified thresholds.

=cut

