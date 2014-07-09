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
}

if (@rows == 0) {
	die "EMPTY MATRIX: no sequences had less missing data than the specified threshold\n";
}
my $total_length = length($rows[0]);
print "there are " . @final_rows . " final rows, $total_length columns\n";
# second pass: remove columns with excessive missing data
# if more than $col_thresh * @final_rows Ns in a column, delete and move on.
my $max_col_ambigs = int($col_thresh * @final_rows);

my $deleted_cols = 0;

@final_rows = @{check_block(\@rows, $total_length/2, $max_col_ambigs)};
$deleted_cols = $total_length - length (@final_rows[0]);
print "deleted $deleted_rows rows, $deleted_cols cols\n";
open FH, ">", "$outname.fasta";
for (my $i=0;$i<@rows;$i++) {
	my $row = $final_rows[$i];
	print FH ">$rowids[$i]\n$row\n";
}
close FH;

# divide the matrix into the first 1000 cols:
# if the number of Ns is less than $max_col_ambigs, this group is fine. No deletions.
# else divide by half, recurse.

sub check_block {
	my $seq_array = shift;
	my $numcols = shift;
	my $max_ambig = shift;

	my $block_missing = 0;
	# if we didn't get an array, return 0. (quick exit)
	if (($seq_array == 0) || ($seq_array == ())) {
		print "\n";
		return 0;
	}

	my $seqlen = length (@$seq_array[0]);

	print "$seqlen, $numcols - ";

	# if we did get an array, but the length of the strings is 0, then return 0. (quick exit)
	if ($seqlen == 0) {
		print "empty\n";
		return 0;
	}

	# perl regex limit:
	if ($numcols > 32766) {
		print " (regex max) ";
		$numcols = 32766;
	}


	# if the block is narrower than the numcols, replace numcols with the width of the block
	if ($seqlen < $numcols) {
		print "REPLACE";
		$numcols = $seqlen;
	}

	# actual base case:
	# if the block is only one col wide, check for Ns and then return either the array or 0.
	if ($seqlen == 1) {
		my $col = join (",", @$seq_array);
		$block_missing = ($col =~ tr/Nn\-\?//);
		print "col has $block_missing missing (max $max_ambig): ";
		if ($block_missing < $max_ambig) {
			print "KEEP the col\n";
			return $seq_array;
		} else {
			print "DELETE a col\n";
			return 0;
		}
	}

	# now for the recursive step:
	# count the number of Ns in this block.
	foreach my $row (@$seq_array) {
		$row =~ /(.{$numcols})(.*)/;
		$block_missing += ($1 =~ tr/Nn\-\?//);
	}

	print "block of $numcols has $block_missing missing (max $max_ambig): ";
	if ($block_missing < $max_ambig) {
		print "KEEP BLOCK\n";
		# the block is fine. No deletions.
	} else {
		# split this block into two parts: recurse on the first part, then the second part
		# then merge those two finished blocks.
		my $new_numcols = int ($numcols/2);
		my $front_block = ();
		foreach my $row (@$seq_array) {
			$row =~ /(.{$numcols})(.*)/;
			push @$front_block, $1;
			$row = $2;
		}
		print "TWO BLOCKS of " . length(@$front_block[0]) . ", ".length(@$seq_array[0])."\n";
		$front_block = check_block ($front_block, $new_numcols, $max_ambig);
		$seq_array = check_block ($seq_array, $seqlen - $new_numcols, $max_ambig);

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

