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
    for (my $i = 0; $i < @files; $i++) {
        if (scalar @{@inputs[$i]} != $num_entries) {
            die "Error: files have different numbers of inputs." . scalar @{@inputs[$i]};
        }
        if ($has_header) {
            my @heads = split /\t/, (@{@inputs[$i]}[0]);
#             foreach $head (@heads) {
#                 $head = @files[$i] . "|" . $head;
#             }
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
    my %labels;
    if ($labelfile) {
        open FH, "<", "$labelfile" or die "make_label_lookup died: couldn't open $labelfile\n";
        my @items = <FH>;
        close FH;
        foreach my $line (@items) {
            (my $name, my $label) = split (/\t/,$line);
            $label =~ s/\r|\n//;
            $labels{$name} = $label;
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

sub get_iupac_code {
	my $arg = shift;

	my $charstr = $arg;
	if (ref($arg) =~ /ARRAY/) {
		$charstr = join ("",@$arg);
	}
	$charstr =~ s/\W//g;
	$charstr =~ s/_//g;
	$charstr =~ s/\s//g;
	$charstr = uc($charstr);
	$charstr = join ("",sort(split('',$charstr)));

	if ($charstr eq "AC") {
		return "M";
	} elsif ($charstr eq "AG") {
		return "R";
	} elsif ($charstr eq "AT") {
		return "W";
	} elsif ($charstr eq "CG") {
		return "S";
	} elsif ($charstr eq "CT") {
		return "Y";
	} elsif ($charstr eq "GT") {
		return "K";
	} elsif ($charstr eq "ACG") {
		return "V";
	} elsif ($charstr eq "ACT") {
		return "H";
	} elsif ($charstr eq "AGT") {
		return "D";
	} elsif ($charstr eq "CGT") {
		return "B";
	} elsif ($charstr eq "ACGT") {
		return "N";
	} else {
		return $charstr;
	}
}


# must return 1 for the file overall.
1;
