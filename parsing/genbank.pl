use Data::Dumper;

my $gbfile = shift;

open FH, "<", $gbfile;

my $gene_hash = {};
my $line = "";
my $id_count = 0;
my $gene_features = "";
my $gene_id = 0;
my $features = 0;
while (defined $line) {
	$line = readline FH;

	# if the line is blank, skip to the next line
	if ($line =~ /^\s+$/) {
		next;
	}

	# if we hit a line that says FEATURES, we're starting the features.
	if ($line =~ /FEATURES/) {
		print "starting features\n";
		$features = 1;
		next;
	}

	if ($features == 1) {
		if ($line =~ /source/) {
			# source feature. Do nothing.
		} elsif ($line =~ /^\s*gene\s+(.+)$/) {
			# start of the next gene
			if ($gene_features ne "") {
				# don't parse features if the feature happens before we start a gene:
				if ($gene_id != 0) {
					# we already were processing a gene: finish it up.
					$gene_hash->{$gene_id}->{"features"} = parse_gene_features ($gene_id, $gene_features);
				}

				# reset for new gene:
				$gene_features = "";
				$id_count++;
			}

			# ok, we've processed any previous genes. Start the next gene.
			$gene_id = $id_count;
			my $gene_regions = parse_interval($1);
			if (@$gene_regions == 1) {
				# should be a single interval.
				$gene_hash->{$gene_id}->{"region"} = @$gene_regions[0];
			} else {
				# this is a badly-formed gene.
			}
		} elsif ($line =~ /ORIGIN/) {
			# this is the end of the features
			last;
		} else {
			$gene_features .= $line;
		}
	}
}

#finish up for the last gene feature.
$gene_hash->{$gene_id}->{"features"} = parse_gene_features ($gene_id, $gene_features);

close FH;
print Dumper($gene_hash) . "\n";

sub parse_interval {
	my $intervalstr = shift;
	my @regions = ();
	if ($intervalstr =~ /^complement\s*\((.+)\)/) {
		# this is a complementary strand feature.
		my $subregions = parse_interval($1);
		foreach my $subreg (@$subregions) {
			if ($subreg =~ /(\d+)\.\.(\d+)/) {
				push @regions, "$2..$1";
			}
		}
	} elsif ($intervalstr =~ /^join\s*\((.+)\)$/) {
		# this is a series of intervals
		my @subintervals = split(/,/, $1);
		foreach my $subint (@subintervals) {
			my $subregions = parse_interval($subint);
			push @regions, @$subregions;
		}
	} elsif ($intervalstr =~ /(\d+)\.\.(\d+)/) {
		push @regions, "$intervalstr";
	}
	return \@regions;
}

sub parse_gene_features {
	my $id = shift;
	my $gene_features = shift;
	print "Hey, $id has some features!\n$gene_features";
	return $gene_features;
}
