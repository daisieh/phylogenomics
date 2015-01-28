package Subfunctions;
use strict;
use FindBin;
use lib "$FindBin::Bin/..";
use Nexus qw(parse_nexus);
use Blast qw (parse_xml compare_hsps compare_regions);
use Genbank qw (parse_genbank flatten_features write_features_as_fasta sequence_for_interval parse_regionfile parse_featurefile set_sequence sequin_feature);
use Data::Dumper;


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(timestamp combine_files make_label_lookup sample_list get_ordered_genotypes get_allele_str get_iupac_code reverse_complement parse_fasta write_fasta parse_phylip write_phylip meld_matrices sortfasta meld_sequence_files vcf_to_depth blast_to_alignment blast_short_to_alignment system_call disambiguate_str split_seq line_wrap trim_to_ref align_to_ref align_to_seq subseq_from_fasta translate_seq codon_to_aa pad_seq_ends set_debug debug find_sequences consensus_str blast_to_genbank merge_to_featuretable);
}

my $debug = 0;

sub set_debug {
	$debug = shift;
}

sub debug {
	my $msg = shift;
	if ($debug) {
		print $msg;
	}
}

=head1

B<(String $time, String $date) timestamp ()>

Convenience function to get the current time and date as formatted strings.

=cut

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return ($time, $date);
}

=head1

B<String $result_str combine_files ( Array \@files, Boolean $has_names, Boolean $has_header)>

Takes any number of input files containing tab-delimited lists of the same length
and creates a tab-delimited string with each list as a column.

@files:     a pointer to a list of filenames
$has_names: 1 if files have the same column of names (first column of each file), otherwise 0.
$has_header:1 if files have column labels in first row, otherwise 0.

Default is no column names, no common row names.

=cut

sub combine_files {
    my $fileptr = shift;
    my $has_names = shift;
    my $has_header = shift;
    my $out_file = shift;

    my @files = @$fileptr;

    if (@files < 1) { die "no files provided."; }

    my @inputs;
    for (my $i=0; $i<@files; $i++) {
        open FH, "<", @files[$i] or die "can't open @files[$i]\n";
        my @data = <FH>;
        close FH;
        push @inputs, \@data;
    }

    my $result = "";
    my @labels = ();
    my $num_entries = scalar @{@inputs[0]};
    my $inputhash = {};
    for (my $i = 0; $i < @files; $i++) {
        if (scalar @{@inputs[$i]} != $num_entries) {
            die "Error: files have different numbers of inputs. " . scalar @{@inputs[$i]} . "!=" . $num_entries;
        }
        if ($has_header) {
            my @heads = split /\t/, (@{@inputs[$i]}[0]);
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
    return $result;
}

sub combine_files_old {
    my $fileptr = shift;
    my $has_names = shift;
    my $has_header = shift;
    my $out_file = shift;

    my @files = @$fileptr;

    if (@files < 1) { die "no files provided."; }

    my @inputs;
    for (my $i=0; $i<@files; $i++) {
        open FH, "<", @files[$i] or die "can't open @files[$i]\n";
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
            @{@inputs[$i]}[0] = join ("\t", @heads);
        }
    }

    #if different files have same column names, must disambiguate column names.
#     my @sortedheads = sort(@heads);
#     my $flag = 0;
#     my $common_name = @sortedheads[0];
#     for(my $j=1; $j<@sortedheads; $j++) {
#         if ($common_name eq @sortedheads[$j]) {
#
#         }
#     }
#

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
    return $result;
}

=head1

B<Hashref make_label_lookup ( String $labelfile )>

Given a tab-delimited file of sample ids and human-readable labels, returns
a hash ref for quick lookup.

$labelfile:   tab-delimited file of sample ids and human-readable labels

=cut

sub make_label_lookup {
    my $labelfile = shift;
    my $reverse = shift;
    my %labels;
    if ($labelfile) {
        open FH, "<", "$labelfile" or die "make_label_lookup died: couldn't open $labelfile\n";
        my @items = <FH>;
        close FH;
        foreach my $line (@items) {
            (my $name, my $label) = split (/\t/,$line);
            $label =~ s/\r|\n//;
    		if ($reverse) {
    			$labels{$label} = $name;
    		} else {
				$labels{$name} = $label;
			}
        }
    }
    return \%labels;
}

=head1

B<@String sample_list ( String $samplefile )>

Given a list of samples in a file, return the samples in an array. If there
are file extensions, they are removed.

$samplefile:   samples listed in a file

=cut

sub sample_list {
    my $samplefile = shift;
    my @samples = ();
    if ($samplefile) {
    	if (-e $samplefile) {
			open FH, "<", "$samplefile" or die "sample_list died: couldn't open $samplefile\n";
			my @items = <FH>;
			close FH;
			foreach my $line (@items) {
				chomp $line;
				if ($line =~ /(.*)\.(.*?)$/) {
					if (length($2) < 5) { # only chop off file extensions if they're less than 5 chars long
						$line = $1;
					}
				}
				push @samples, $line;
			}
		} else {
			$samplefile =~ s/\..*?$//;
			push @samples, $samplefile;
		}
    }
    return \@samples;
}

=head1

B<@String get_ordered_genotypes ( String $charstr )>

Returns the ordered list of diploid genotypes possible given the allele string provided.
Genotype ordering: Ref allele is index 0, alt alleles are indexed from there.
For each combination j,k, the ordering is:
F(j,k) = (k*(k+1)/2)+j

$charstr:   string or array of alleles to order into genotypes

=cut


sub get_ordered_genotypes {
	my $charstr = shift;

	my @alleles = split('',$charstr);
	my @genotypes = ();
	for(my $j=0;$j<@alleles;$j++) {
		for (my $k=$j;$k<@alleles;$k++) {
			my $order = ($k*($k+1)/2)+$j;
			@genotypes[$order] = @alleles[$j] . @alleles[$k];
		}
	}
	return \@genotypes;
}

=head1

B<String disambiguate_str ( String/Array $charstr )>

Convenience function to stringify a list or a string of alleles with extra characters.
Returns the string in alphabetical order and in all caps.

$charstr:   string or array of alleles to stringify

=cut


sub disambiguate_str {
	my $arg = shift;

	my $charstr = $arg;
	if (ref($arg) =~ /ARRAY/) {
		$charstr = join ("",@$arg);
	}
	if (length($charstr) == 1) {
		return $charstr;
	}
	if ($charstr =~ /^-+$/) {
		return "-";
	}
	$charstr = uc($charstr);
	$charstr =~ s/N/ACGT/g;
	$charstr =~ s/M/AC/g;
	$charstr =~ s/R/AG/g;
	$charstr =~ s/W/AT/g;
	$charstr =~ s/S/CG/g;
	$charstr =~ s/Y/CT/g;
	$charstr =~ s/K/GT/g;
	$charstr =~ s/V/ACG/g;
	$charstr =~ s/H/ACT/g;
	$charstr =~ s/D/AGT/g;
	$charstr =~ s/B/CGT/g;

	$charstr = join ("",sort(split('',$charstr)));
	return "$charstr";
}

# returns the consensus sequence as a string: wraps the recursive implementation
sub consensus_str {
	my $seq_array = shift;
	pad_seq_ends ($seq_array,"-");

	return shift @{consensus($seq_array)};
}

# recursive implementation: returns an array with one element, the final consensus sequence.
sub consensus {
	my $seq_array = shift;

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

	# check to see if all sequences are identical:
	my $identical = 1;
	my $curr_seq = @$seq_array[0];
	foreach my $seq (@$seq_array) {
		if ($seq ne $curr_seq) {
			$identical = 0;
		}
	}


	if ($identical) {
		debug ("IDENTICAL\n");
		my @result = ();
		push @result, $curr_seq;
		return \@result;
	} elsif ($seqlen == 1) {
		# base case: single character; get consensus.
		my @result = ();
		push @result, get_iupac_code($seq_array);
		debug ("CONSENSUS of " . join(",",@$seq_array) . " is " . $result[0] . "\n");
		return \@result;
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
		foreach my $seq (@$seq_array) {
			$seq =~ /(.{$numcols})(.*)/;
			push @$front_block, $1;
			$seq = $2;
		}
		$front_block = consensus ($front_block);
		$seq_array = consensus ($seq_array);

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


=head1

B<String get_allele_str ( String/Array $charstr )>

Convenience function to stringify a list or a string of alleles with extra characters.
Returns the string in alphabetical order and in all caps.

$charstr:   string or array of alleles to stringify

=cut


sub get_allele_str {
	my $arg = shift;

	my $charstr = $arg;
	if (ref($arg) =~ /ARRAY/) {
		$charstr = join ("",@$arg);
	}
	if (length($charstr) == 1) {
		return $charstr;
	}
	if ($charstr =~ /^-+$/) {
		return "-";
	}
	$charstr = uc($charstr);
	$charstr =~ s/\W//g;
	$charstr =~ s/_//g;
	$charstr =~ s/\s//g;
	$charstr =~ tr/[A-Z]//c;
	$charstr =~ tr/UX/TN/;
	$charstr =~ tr/ABCDGHMNRSTVWY//c;

	$charstr = join ("",sort(split('',$charstr)));
	return "$charstr";
}


=head1

B<String get_iupac_code ( String/Array $charstr )>

Convenience function to return the iupac ambiguity code for whatever alleles are inputted.

$charstr:   string or array of alleles.

=cut

sub get_iupac_code {
	my $arg = shift;

	# regularize the input first, so that it's an alphabetical, no-dups, sorted uc string
	my $charstr = get_allele_str ($arg);
	$charstr =~ tr/A-Z//s;
	if (length($charstr) == 1) {
		return $charstr;
	}
	if ($charstr =~ /-+/) {
		return "-";
	}
	while (length ($charstr) > 1) {
		if ($charstr =~ /N/) {
			 return "N";
		}
		# MRWSYKVHDBN
		$charstr =~ s/ACGT/N/g;
		$charstr =~ s/AC/M/g;
		$charstr =~ s/AG/R/g;
		$charstr =~ s/AT/W/g;
		$charstr =~ s/CG/S/g;
		$charstr =~ s/CT/Y/g;
		$charstr =~ s/GT/K/g;
		$charstr =~ s/MG/V/g;
		$charstr =~ s/MT/H/g;
		$charstr =~ s/RT/D/g;
		$charstr =~ s/ST/B/g;
		$charstr =~ s/MK/N/g;
	}
	return $charstr;
}

=head1

B<String reverse_complement ( String $charstr )>

Convenience function to return the reverse complement of a sequence.

$charstr:   sequence to revcomp.

=cut


sub reverse_complement {
	my $charstr = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($charstr);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}


=head1

B<(\%taxa, \@taxanames) parse_fasta ( Filehandle $inputfile )>

Given a fasta file as input, returns a hash containing all the sequences, keyed by the
values of the taxanames array.

$inputfile:   fasta file to parse.

=cut


sub parse_fasta {
	my $fastafile = shift;
	my $no_substitution = shift;

	my $taxa = {};
	my @taxanames = ();
	my $length = 0;
	open fileIN, "<", "$fastafile" or die "couldn't parse fasta file $fastafile";

	my $input = readline fileIN;
	my $taxonlabel = "";
	my $sequence = "";
	while (defined $input) {
		if ($input =~ /^>(.+)$/) {
			$taxonlabel = $1;
			$taxonlabel =~ s/\s+$//;
			if (!defined $no_substitution) {
				$taxonlabel =~ s/[\s\\\/,;\-]+/_/g;
			}
			if (exists $taxa->{$taxonlabel}) {
				$taxa->{$taxonlabel} = "";
			} else {
				push @taxanames, $taxonlabel;
			}
			if ($length > 0) {
				# we are at the next taxon; push the last one onto the taxon array.
				$taxa->{"length"} = $length;
				$length = 0;
			}
		} else {
			$input =~ /^\s*(.+)\s*$/;
			$taxa->{$taxonlabel} .= $1;
			$length += length($1);
		}
		$input = readline fileIN;
	}

	close (fileIN);
	return $taxa, \@taxanames;
}

sub write_fasta {
	my $fastahash = shift;
	my $fastanames = shift;
	my $result = "";
	if (exists $fastahash->{"characters"}) {
		unless (exists $fastahash->{"names"}) {
			$fastahash->{"names"} = (keys %{$fastahash->{"characters"}});
		}

		foreach my $t (@{$fastahash->{"names"}}) {
			$result .= ">$t\n$fastahash->{characters}->{$t}\n";
		}
	} elsif ((ref $fastanames) eq "ARRAY") {
		foreach my $name (@{$fastanames}) {
			$result .= ">$name\n$fastahash->{$name}\n";
		}
	} else {
		foreach my $name (keys %{$fastahash}) {
			$result .= ">$name\n$fastahash->{$name}\n";
		}
	}
	return $result;
}

sub find_sequences {
	my $fastafile = shift;
	my $names = shift;

	unless (-e $fastafile) {
		die "File $fastafile does not exist.\n";
	}

	my $hashed_seqs = {};
	my ($taxa, $taxanames) = parse_fasta ($fastafile);

	foreach my $name (@$names) {
		if (exists $taxa->{$name}) {
			$hashed_seqs->{$name} = $taxa->{$name};
		}
	}
	return $hashed_seqs;
}


sub pad_seq_ends {
	my $seq_block = shift;
	my $char = shift;

	if (!(defined $char)) {
		$char = "n";
	}

	# find the maximum length of any sequence
	my $max_len = 0;
	foreach my $line (@$seq_block) {
		my $len = length ($line);
		if ($len > $max_len) {
			$max_len = $len;
		}
	}

	foreach my $line (@$seq_block) {
		if (length ($line) < $max_len) {
			$line = $line . $char x ($max_len - length ($line));
		}
	}

	return $max_len;
}

=head1

B<(\%taxa, \@taxanames) parse_phylip ( String $inputfile )>

Given a phylip file as input, returns a hash containing all the sequences, keyed by the
values of the taxanames array.

$inputfile:   phylip file to parse.

=cut

sub parse_phylip {
	my $inputfile = shift;

	my $taxa = {};
	my @taxonlabels = ();

	open FH, "<", $inputfile or die "couldn't find file $inputfile";

	my $line = readline FH;
	while ($line =~ /^\s*$/) {
		$line = readline FH;
	}
	$line =~ /\s*(\d+)\s+(\d+)/;
	my $num_taxa = $1;
	my $num_chars = $2;
	my $count = 0;
	while ($line = readline FH) {
		if ($line =~ /^\s*$/) {
			next;
		}

		# if we haven't seen all of the taxa labels
		if (@taxonlabels < $num_taxa) {
			$line =~ /^(.+?)\s+(.*)$/;
			my $taxon = $1;
			my $chars = $2;
			$chars =~ s/\s//g;
			$taxa->{$taxon} = $chars;
			push @taxonlabels, $taxon;
		} else {
			$line =~ s/\s//g;
			$taxa->{$taxonlabels[$count++]} .= $line;
			if ($count == $num_taxa) {
				$count = 0;
			}
		}
	}
	close FH;
	return $taxa, \@taxonlabels;
}

sub write_phylip {
	my $taxa = shift;
	my $taxanames = shift;

	my $result = "";
	my @workingseqs = ();
	my @trunc_names = ();
	foreach my $t (@$taxanames) {
		push @workingseqs, $taxa->{$t};
		my $working_t = $t;
		$working_t =~ s/\s/_/g;
		if ($working_t =~ /^(.{10}).*$/) {
			$working_t = $1;
		} else {
			$working_t .= $working_t . " " x (10 - length $working_t);
		}
		push @trunc_names, $working_t;
	}

	my $nchar = pad_seq_ends(\@workingseqs);
	my $ntax = @$taxanames;

	#first line is ntax nchar
	$result .= "$ntax $nchar\n";

	my $blocksize = 60;
	#print the first set of chars with names
	for (my $i=0; $i < @$taxanames; $i++) {
		my $seq = $workingseqs[$i];
		if ($seq =~ /^(.{$blocksize})(.*)$/) {
			$seq = $1;
			$workingseqs[$i] = $2;
		} else {
			$workingseqs[$i] = "";
		}
		$seq =~ s/(.{10})/$1 /g;
		$result .= "$trunc_names[$i] $seq\n";
	}
	while (length $workingseqs[0] > 0) {
		$result .= "\n";
		for (my $i=0; $i < @$taxanames; $i++) {
			my $seq = $workingseqs[$i];
			if ($seq =~ /^(.{$blocksize})(.*)$/) {
				$seq = $1;
				$workingseqs[$i] = $2;
			} else {
				$workingseqs[$i] = "";
			}
			$seq =~ s/(.{10})/$1 /g;
			$result .= " "x11 ."$seq\n";
		}
	}
	return $result;
}

sub line_wrap {
	my $sequence = shift;
	my $linelength = shift;

	my $result = "";
	while ($sequence =~ /^(.{$linelength})(.*)$/) {
		$result .= "$1\n";
		$sequence = $2;
	}
	$result .= $sequence;
	return $result;
}

=head1

B<(\%mastertaxa, \%regiontable) meld_matrices ( @matrixnames, %matrices )>

Given a hash of sequence matrices indexed by the values of @matrixnames, melds them into
a single hash of concatenated sequences. The regiontable hash contains the information about
which taxa contained which sequences and where they are in the concatenated supermatrix.

@matrixnames:   Names of the matrices used as keys to the hash.
%matrices       The sequence matrices to be concatenated, indexed by the values of @matrixnames.

=cut

sub meld_matrices {
	my $arg1 = shift;
	my $arg2 = shift;

	my @matrixnames = @$arg1;
	my %matrices = %$arg2;
	my $currlength = 0;

	# start the master matrix: for every taxon in every input file, make a blank entry.
	my %mastertaxa = ();
	my %regiontable = ();
	foreach my $inputfile ( keys (%matrices) ) {
		foreach my $taxon ( keys (%{$matrices{$inputfile}})) {
			$mastertaxa {$taxon} = "";
		}
	}

	# now, in order of the inputted matrix names, add the sequences to the taxa of the master matrix.
	# if a taxon is missing from the matrix, add missing data for that entry.
	foreach my $key (@matrixnames) {
		my $ref = $matrices{$key};
		my @curr_matrix_taxa = keys(%$ref);
		my $total = length($ref->{$curr_matrix_taxa[0]});
		foreach my $v (values %$ref) {
			if (length($v) > $total) {
				$total == length($v);
			}
		}
		$regiontable{"regions"} .= "$key\t";
		my %expandedmatrix = ();
		foreach my $k (keys %mastertaxa) {
			#add entries from this matrix into expandedmatrix
			if (defined $ref->{$k}) {
				$mastertaxa{$k} .= $ref->{ $k };
				$regiontable{$k} = $regiontable{$k} . "x\t";
			} else {
				$mastertaxa{$k} .= "-" x $total;
				$regiontable{$k} = $regiontable{$k} . "\t";
			}
		}
		my $starts_at = $currlength + 1;
		$currlength = $currlength + $total;
		$regiontable{"exclusion-sets"} = $regiontable{"exclusion-sets"} . ($currlength + 1) . "-" . "$currlength\t";

	}
	$mastertaxa{"length"} = $currlength;

	return (\%mastertaxa, \%regiontable);
}

sub sortfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	if ($separator) {
		print "using $separator\n";
	} else {
		$separator = '\n';
	}

	my (undef, $tempfile) = tempfile(UNLINK => 1);
	system ("gawk '{if (NF==0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | sort -t',' -k 1 | gawk '{print \">\" \$1 \"$separator\" \$2}' FS=\",\" > $outfile");
}


=head1

B<(\%mastertaxa, \%regiontable) meld_matrices ( @inputfiles )>

Given a list of sequence files, melds them into a single hash of concatenated sequences.
The regiontable hash contains the information about which taxa contained which sequences
and where they are in the concatenated supermatrix.

@inputfiles:    An array of file names.

=cut

sub meld_sequence_files {
	my $arg = shift;
	my @inputfiles = @$arg;

	my $matrices = {};
	my @matrixnames = ();

	foreach my $inputfile (@inputfiles) {
		push @matrixnames, $inputfile;
		if ($inputfile =~ /\.nex/) {
			my $nexushash = parse_nexus ($inputfile);
			$matrices->{ $inputfile } = $nexushash->{"characters"};
		} elsif ($inputfile =~/\.fa/) {
			($matrices->{ $inputfile }, undef) = parse_fasta ($inputfile);
		} else {
			print "Couldn't parse $inputfile: not nexus or fasta format\n";
		}
	}

	return meld_matrices (\@matrixnames, $matrices);
}


=head1

B<String $out_file meld_matrices ( String $vcf_file, [String $out_file] )>

Generates a summary file with the read depths of each position in the inputted vcf file.

$vcf_file:  A vcf file.

=cut



sub vcf_to_depth {
	my $vcf_file = shift;
	my $out_file = shift;

	$vcf_file =~ /(.+)(\.vcf)/;
	my $basename = $1;

	if ($out_file == 0) {
		$out_file = "$basename.depth";
	}

	system ("awk 'BEGIN {OFS=\"\\t\"} /^.*?\\t(.*?)\\t/ {print \$1,\$2,\$8}' $vcf_file | awk 'BEGIN {OFS=\"\\t\"} {sub(/DP=/,\"\",\$3);sub(/;.*/,\"\",\$3);print \$1,\$2,\$3}' - > $out_file");

	return $out_file;
}

sub blast_to_alignment {
	my $blastfile = shift;

	my @seqids = ();
	my @seqs = ();
	my %result_matrix = ();
	my $curr_query_id = "";
	my $query_seq = "";
	my ($query_end, $subject_end, $ref_length) = 0;
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	# processing the header
	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			$curr_query_id = $1;
		}
		if ($line =~ /Subject/) {
			while ($line = shift @lines) {
				if ($line =~ /Length=(\d+)/) {
					$ref_length = $1;
					last;
				}
			}
			last;
		}
	}

	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			if ($query_seq ne "") {
				$query_seq .= "n" x ($ref_length - length($query_seq));
				$result_matrix{$curr_query_id} = $query_seq;
			}
			$curr_query_id = $1;
			($query_end, $subject_end) = 0;

			$query_seq = "";
		} elsif ($line =~ /Query\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			$curr_query_start = $1;
			$curr_query_seq = $2;
			$curr_query_end = $3;
			if ($curr_query_start < $query_end) {
				while ($line !~ /Query=\s+(.*)$/) {
					$line = shift @lines;
				}
				unshift @lines, $line;
			}
		} elsif ($line =~ /Sbjct\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			my $curr_subject_start = $1;
			my $curr_subject_seq = $2;
			my $curr_subject_end = $3;
			while ($curr_subject_seq =~ /(.*)-(.*)/) {
				$curr_subject_seq = "$1$2";
				$curr_query_seq =~ /(.{length($1)}).(.*)/;
				$curr_query_seq = "$1$2";
			}
			if ($query_seq eq "") {
				$query_seq = "n" x ($curr_subject_start - 1);
				$query_seq .= uc($curr_query_seq);
			} elsif (($query_end + 1) == $curr_query_start) {
				$query_seq .= uc($curr_query_seq);
			} else {
				$query_seq .= "n" x ($curr_subject_start - $subject_end - 1);
				$query_seq .= uc($curr_query_seq);
			}
			$subject_end = $curr_subject_end;
			$query_end = $curr_query_end;
		}
	}
	$query_seq .= "n" x ($ref_length - length($query_seq));
	$result_matrix{$curr_query_id} = $query_seq;

	return \%result_matrix;

}

sub blast_short_to_alignment {
	my $blastfile = shift;

	my @seqids = ();
	my @seqs = ();
	my %result_matrix = ();
	my $curr_query_id = "";
	my $query_seq = "";
	my ($query_end, $subject_end, $ref_length) = 0;
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	# processing the header
	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			$curr_query_id = $1;
		}
		if ($line =~ /Subject/) {
			while ($line = shift @lines) {
				if ($line =~ /Length=(\d+)/) {
					$ref_length = $1;
					last;
				}
			}
			last;
		}
	}

	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			if ($query_seq ne "") {
				my $ns = "n" x ($ref_length - length($query_seq));
				$query_seq .= "n" x ($ref_length - length($query_seq));
				$result_matrix{$curr_query_id} = $query_seq;
			}
			$curr_query_id = $1;
			push @seqids, $curr_query_id;
			$query_seq = "";
		} elsif ($line =~ /Query\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			$curr_query_seq = $2;
			$curr_query_end = $3;
		} elsif ($line =~ /Sbjct\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			my $curr_subject_start = $1;
			my $curr_subject_seq = $2;
			my $curr_subject_end = $3;
			while ($curr_subject_seq =~ /(.*)-(.*)/) {
				$curr_subject_seq = "$1$2";
				$curr_query_seq =~ /(.{length($1)}).(.*)/;
				$curr_query_seq = "$1$2";
			}
			if ($query_seq eq "") {
				$query_seq = "n" x ($curr_subject_start - 1);
				$query_seq .= uc($curr_query_seq);
			}
		}
	}
	$query_seq .= "n" x ($ref_length - length($query_seq));
	$result_matrix{$curr_query_id} = $query_seq;

	return \%result_matrix;

}

sub system_call {
	my $cmd = shift;

	my $exit_val = eval {
		system ($cmd);
	};

	if ($exit_val != 0) {
		print "System call \"$cmd\" exited with $exit_val\n";
		exit;
	}
	return $exit_val;
}

=head1

B<String ($startseq, $regionseq, $endseq) split_seq ( String $seq, int $start, int $end )>

=cut
### taken from aTRAM:

sub split_seq {
    my $seq = shift;
    my $start = shift;
    my $end = shift;
    my $max = 30000;
	my $seqlen = length ($seq);
	my $startseq = "";
	my $regionseq = "";
	my $endseq = "";

	my $currstart = $start-1;
	my $currend = $end;
	while ($currstart > $max) {
			$seq =~ /^(.{$max})(.*)$/;
			$startseq .= $1;
			$seq = $2;
			$currstart -= $max;
			$currend -= $max;
	}
	if ($currstart > 0) {
			$seq =~ /^(.{$currstart})(.*)$/;
			$startseq .= $1;
			$seq = $2;
			$currstart = 1;
			$currend = $end - (length ($startseq));
	}

	my $regionsize = $end - $start + 1;
	while ($regionsize > $max) {
			$seq =~ /^(.{$max})(.*)$/;
			$regionseq .= $1;
			$seq = $2;
			$currstart -= $max;
			$currend -= $max;
			$regionsize -= $max;
	}
	if ($regionsize > 0) {
			$seq =~ /^(.{$regionsize})(.*)$/;
			$regionseq .= $1;
			$endseq = $2;
	}
	return ($startseq, $regionseq, $endseq);
}

sub align_to_ref {
	my $reffile = shift;
	my $contigfile = shift;
	my $outname = shift;
	my $aligner = shift;

	# check that the files we need exist:
	unless ((-e $reffile) && (-e $contigfile)) {
		return undef;
	}

	open REF_FH, "<:crlf", $reffile;
	my $ref = readline REF_FH;
	$ref =~ />(.+)$/;
	my $refname = $1;
	close REF_FH;

	my (undef, $catfile) = tempfile(UNLINK => 1);
	`cat $reffile $contigfile > $catfile`;

	return align_to_seq ($catfile, $refname, $outname, $aligner);
}

sub align_to_seq {
	my $catfile = shift;
	my $refname = shift;
	my $outname = shift;
	my $aligner = shift;

	my $result = 0;
	if ($aligner eq "mafft") {
		$result = run_command (get_bin("mafft"), "$catfile > $outname.align.fasta");
	} else {
		$aligner = "muscle";
		$result = run_command (get_bin("muscle"), "-in $catfile -out $outname.align.fasta");
	}

	if ($result == -1) {
		printlog ("Aligner $aligner could not be found.");
		exit -1;
	}

	if ($result != 0) {
		printlog ("Aligner $aligner failed.");
		return undef;
	}

	my ($seq_hash, undef) = parsefasta ("$outname.align.fasta");
	my $gappedrefseq = delete $seq_hash->{"$refname"};
	$seq_hash->{"reference"} = $gappedrefseq;

	return ($seq_hash);
}

sub trim_to_ref {
	my $seq_hash = shift;
	my $refseq_key = shift;

	my $gappedrefseq = delete $seq_hash->{$refseq_key};

	my $result = 0;
	# as long as there are still gaps in the reference sequence, keep removing the corresponding positions from the contigs.
	while ($gappedrefseq =~ /(\w*)(-+)(.*)/) {
		my $left = $1;
		my $gap = $2;
		my $remainder = $3;
		my $leftlength = length $left;
		my $gaplength = length $gap;

		foreach my $contig (keys %$seq_hash) {
			my $end = $leftlength + $gaplength;
			my $start = $leftlength + 1;
			my ($startseq, $regionseq, $endseq) = split_seq ($seq_hash->{$contig}, $start, $end);
			$seq_hash->{$contig} = "$startseq$endseq";
		}

		$gappedrefseq = "$left$remainder";
	}
	# align the ends of the contig seqs
	foreach my $contig (keys %$seq_hash) {
		if ((length $seq_hash->{$contig}) > (length $gappedrefseq)) {
			# if the contig seq is longer than the gappedrefseq, truncate it
			$seq_hash->{$contig} = substr($seq_hash->{$contig}, 0, length $gappedrefseq);
		} else {
			# if the contig seq is shorter than the gappedrefseq, pad it with gaps.
			$seq_hash->{$contig} .= "-" x ((length $seq_hash->{$contig}) - (length $gappedrefseq));
		}
	}

	$seq_hash->{$refseq_key} = $gappedrefseq;
	return $seq_hash;
}

### end aTRAM functions

=head1

B<String $out_file subseq_from_fasta ( String $fastafile, int $start, int $end )>

For use in dealing with extremely large fasta files, so we don't have to actually read in
the whole sequence into memory before finding a subseq.

=cut


sub subseq_from_fasta {
	my $fastafile = shift;
	my $start = shift;
	my $end = shift;

	my $sequence = "";
	my $pos = 0;
	my $newstart = 0;
	my $length = $end - $start;
	my $newend = 0;
	open FH, "<", $fastafile or die "couldn't open $fastafile";

	my $line = readline FH; # first line is the name
	$line = readline FH;
	while (defined $line) {
		$line =~ s/\s//g;
		my $linelen = length ($line);
		if (($pos + $linelen) >= $start) {
			if ($newstart == 0) {
				$newstart = $start - $pos;
				$newend = $newstart + $length;
			}
			$sequence .= "$line";
			if ($pos >= $end) {
				last;
			}
		}
		$pos += $linelen;

		$line = readline FH;
	}
	close FH;
	my (undef, $finalseq, undef) = split_seq ($sequence, $newstart, $newend);

	return $finalseq;
}

sub translate_seq {
	my $nucl_seq = shift;

	my $aa_seq = "";
	while ($nucl_seq =~ /^(\w\w\w)(.*)$/) {
		$nucl_seq = $2;
		$aa_seq .= codon_to_aa ($1);
	}
	if ($nucl_seq ne "") {
		# this wasn't div by 3, so there's a problem.
		return "";
	}
	return $aa_seq;
}

# reference table:
# F	TTT	TTC
# L	TTA	TTG
# S	TCT	TCC	TCA	TCG
# Y	TAT	TAC
# *	TAA	TAG	TGA
# C	TGT	TGC
# W	TGG
# L	CTT	CTC	CTA	CTG
# P	CCT	CCC	CCA	CCG
# H	CAT	CAC
# Q	CAA	CAG
# R	CGT	CGC	CGA	CGG
# I	ATT	ATC	ATA
# M	ATG
# T	ACT	ACC	ACA	ACG
# N	AAT	AAC
# K	AAA	AAG
# S	AGT	AGC
# R	AGA	AGG
# V	GTT	GTC	GTA	GTG
# A	GCT	GCC	GCA	GCG
# D	GAT	GAC
# E	GAA	GAG
# G	GGT	GGC	GGA	GGG


sub codon_to_aa {
	my $codon = shift;
	$codon = uc($codon);
	if ($codon !~ /.../) { return "-"; }
	if ($codon =~ /AC./) { return "T"; }
	if ($codon =~ /CT./) { return "L"; }
	if ($codon =~ /CC./) { return "P"; }
	if ($codon =~ /CG./) { return "R"; }
	if ($codon =~ /GT./) { return "V"; }
	if ($codon =~ /GC./) { return "A"; }
	if ($codon =~ /GG./) { return "G"; }
	if ($codon =~ /TC./) { return "S"; }
	if ($codon =~ /TT[CTY]/) { return "F"; }
	if ($codon =~ /TT[AGR]/) { return "L"; }
	if ($codon =~ /TA[CTY]/) { return "Y"; }
	if ($codon =~ /TA[AGR]|TGA/) { return "*"; }
	if ($codon =~ /TG[CTY]/) { return "C"; }
	if ($codon =~ /TGG/) { return "W"; }
	if ($codon =~ /CA[CTY]/) { return "H"; }
	if ($codon =~ /CA[AGR]/) { return "Q"; }
	if ($codon =~ /AT[CATMYWH]/) { return "I"; }
	if ($codon =~ /ATG/) { return "M"; }
	if ($codon =~ /AA[CTY]/) { return "N"; }
	if ($codon =~ /AA[AGR]/) { return "K"; }
	if ($codon =~ /AG[CTY]/) { return "S"; }
	if ($codon =~ /AG[AGR]/) { return "R"; }
	if ($codon =~ /GA[CTY]/) { return "D"; }
	if ($codon =~ /GA[AGR]/) { return "E"; }


# 	$charstr =~ s/N/ACGT/g;
# 	$charstr =~ s/M/AC/g;
# 	$charstr =~ s/R/AG/g;
# 	$charstr =~ s/W/AT/g;
# 	$charstr =~ s/S/CG/g;
# 	$charstr =~ s/Y/CT/g;
# 	$charstr =~ s/K/GT/g;
# 	$charstr =~ s/V/ACG/g;
# 	$charstr =~ s/H/ACT/g;
# 	$charstr =~ s/D/AGT/g;
# 	$charstr =~ s/B/CGT/g;
# TAA TAG TGA TAR
	return "X";
}

# returns an array of gene hashes:
# gene->{name}
# gene->{strand}
# gene->{type}
# gene->{contains} = an array of regions
#	[ (start, end), (start, end) ]

sub blast_to_genbank {
	my $gbfile = shift;
	my $fastafile = shift;
	my $outfile = shift;

	my $gene_array = parse_genbank($gbfile);
	my ($ref_hash, $ref_array) = flatten_features($gene_array);
	print Dumper ($gene_array);

	# look for regions too small to blast accurately:
	my $tiny_regions = {};
	my $tiny_region_extension_length = 20;
	foreach my $region (@$ref_array) {
		my $start = $ref_hash->{$region}->{'start'};
		my $end = $ref_hash->{$region}->{'end'};
		if ($end - $start < 10) {
			$tiny_regions->{$region}->{'characters'} = $ref_hash->{"$region"}->{'characters'};
			$start -= $tiny_region_extension_length;
			$end += $tiny_region_extension_length;
			$ref_hash->{"$region"}->{'characters'} = sequence_for_interval ("$start..$end");
		}
	}

	open FAS_FH, ">", "$gbfile.fasta";
	my @new_ref_array = ();
	foreach my $ref (@$ref_array) {
		push @new_ref_array, "$ref\t$ref_hash->{$ref}->{'start'}\t$ref_hash->{$ref}->{'end'}";
		print FAS_FH ">$ref\t$ref_hash->{$ref}->{'start'}\t$ref_hash->{$ref}->{'end'}\n$ref_hash->{$ref}->{'characters'}\n";
	}
	close FAS_FH;

	system("blastn -query $fastafile -subject $gbfile.fasta -outfmt 5 -out $outfile.xml -word_size 10");

	print "parsing results\n";

	# choose the best hits:
	my $hit_array = Blast::parse_xml ("$outfile.xml");
	my $hits = {};
	foreach my $hit (@$hit_array) {
		my $subject = $hit->{"subject"}->{"name"};
		my @hsps = sort Blast::compare_hsps @{$hit->{"hsps"}};
		my $best_hit = shift @hsps;
		$hits->{$subject}->{"hsp"} = $best_hit;
		if ($best_hit->{"hit-from"} < $best_hit->{"hit-to"}) {
			$hits->{$subject}->{"orientation"} = 1;
		} else {
			$hits->{$subject}->{"orientation"} = -1;
		}
	}
	my @result_array = ();
	foreach my $subj (@new_ref_array) {
		$subj =~ s/\t.*$//;
		if (exists $tiny_regions->{$subj}) {
			$hits->{$subj}->{hsp}->{'query-from'} += $tiny_region_extension_length;
			$hits->{$subj}->{hsp}->{'query-to'} -= $tiny_region_extension_length;
		}
		# gene->{name}
# gene->{strand}
# gene->{feature}
# gene->{contains} = an array of regions
#	[ (start, end), (start, end) ]
		$subj =~ /$(\d+)_(\d+)_(.+?)_(.+)$/;
		my $gene_id = $1;
		my $feat_id = $2;
		my $gene_name = $3;
		my $feat_name = $4;
# 		if (!(exists $result_hash->{$gene_name})) {
# 			my $gene_hash = {};
# 			push @result_array, $gene_hash;
# 		}
	}

	open OUTFH, ">", "$outfile.regions" or die "couldn't create $outfile";
	foreach my $subj (@new_ref_array) {
		print OUTFH "$subj($ref_hash->{$subj}->{'strand'})\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\n";
	}

	close OUTFH;

}

sub merge_to_featuretable {
	# a regionfile is the output of parse_blast.pl comparing the fastafile to the reference fasta file from genbank.pl
	my $regionfile = shift;

	# fastafile is the new fasta file to be aligned: same as fasta file used in parse_blast.pl
	my $fastafile = shift;

	# featurefile is the feature table file generated by genbank.pl
	my $featurefile = shift;
	my $outfile = shift;
	my $name = shift;

	if (!defined $name) {
		$name = "";
	}

	my ($gene_array, $gene_index_array) = Genbank::parse_regionfile($regionfile);
	my ($destination_gene_array, $destination_gene_index_array) = Genbank::parse_featurefile($featurefile);
	my ($fastahash, undef) = parse_fasta($fastafile);
	my $seqlen = 0;
	foreach my $k (keys $fastahash) {
		# there should be only one key, so just one name.
		$name = $k;
		$seqlen = length ($fastahash->{$k});
		Genbank::set_sequence($fastahash->{$k});
	}

	my $gene_hash = {};
	foreach my $id (@$gene_index_array) {
		my $gene = shift $gene_array;
		$gene_hash->{$id} = $gene;
	}

	my $dest_gene_hash = {};
	foreach my $id (@$destination_gene_index_array) {
		my $gene = shift $destination_gene_array;
		$dest_gene_hash->{$id} = $gene;
	}

	# fill in the genes from the regionfile with the info from the destination gene array
	my @final_gene_array = ();
	foreach my $id (@$gene_index_array) {
		my $gene = $gene_hash->{$id};
		my $dest_gene = $dest_gene_hash->{$id};

		foreach my $q (keys %{$dest_gene->{"qualifiers"}}) {
			if (!(exists $gene->{"qualifiers"}->{$q})) {
				$gene->{"qualifiers"}->{$q} = $dest_gene->{"qualifiers"}->{$q}
			}
		}
		my @new_contains = ();
		$gene->{"id"} = $id;
	# 	print Dumper($gene);
		foreach my $destcontains (@{$dest_gene->{"contains"}}) {
			my $genecontains = shift $gene->{"contains"};
			$destcontains->{"region"} = $genecontains->{"region"};
			push @new_contains, $destcontains;
		}
		$gene->{"contains"} = \@new_contains;
		push @final_gene_array, $gene;
	}

	open FH, ">", "$outfile.tbl";

	# print header
	print FH ">Features\t$name\n";

	open FASTA_FH, ">", "$outfile.fasta";
	# start printing genes
	foreach my $gene (@final_gene_array) {
		# first, print overall gene information
		my $genename = $gene->{"qualifiers"}->{"gene"};
		my $gene_id = $gene->{"id"};
		my $feat_id = 0;
		foreach my $r (@{$gene->{'region'}}) {
			$r =~ /(\d+)\.\.(\d+)/;
			print FH "$1\t$2\tgene\n";
		}
		foreach my $q (keys %{$gene->{'qualifiers'}}) {
			print FH "\t\t\t$q\t$gene->{qualifiers}->{$q}\n";
		}

		# then, print each feature contained.
		foreach my $feat (@{$gene->{'contains'}}) {
			foreach my $reg (@{$feat->{"region"}}) {
				my $strand = "+";
				my ($start, $end) = split (/\.\./, $reg);
				if ($end < $start) {
					$strand = "-";
					my $oldend = $start;
					$start = $end;
					$end = $oldend;
				}
				my $regseq = Genbank::sequence_for_interval ($reg);
				my $featname = $feat->{"type"};
				print FASTA_FH ">$gene_id"."_$feat_id"."_$genename"."_$featname($strand)\t$start\t$end#$regseq\n";
				$feat_id++;
			}
			print FH Genbank::sequin_feature ($feat->{'region'}, $feat);
		}
	}

	close FH;

}

# must return 1 for the file overall.
1;
