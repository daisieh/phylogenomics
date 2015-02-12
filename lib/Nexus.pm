package Nexus;
use strict;
use Data::Dumper;
use Subfunctions qw(debug set_debug pad_seq_ends);


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(parse_nexus write_nexus_character_block write_nexus_taxa_block write_nexus_trees_block write_nexus_sets_block clean_nexus_taxa_names);
}

=head1

Nexus data structure:
$nexushash->{"characters"}: a hash of character strings keyed by taxon names
$nexushash->{"taxa"}: an array of taxon names
$nexushash->{"trees"}: an array of newick trees
$nexushash->{"charlabels"}: if specified, the names of the characters
$nexushash->{"ntax"}: number of taxa
$nexushash->{"nchar"}: number of characters
$nexushash->{"alttaxa"}: an array of alternate taxon name arrays
$nexushash->{"sets"}: a hash of sets, which are hashed by their set type.
	->{"charset"}: an array of hashes, each representing a charset.
		->{"name"}
		->{"start"}
		->{"end"}

=cut

sub parse_nexus {
	my $inputfile = shift;

	my $nexushash = {};
	my $gapchar = "-";
	my $interleave = 0;
	my $nchar = 0;

	open fileIN, "<:crlf", "$inputfile" or die "no file named $inputfile";
	my @inputs = <fileIN>;

	my $input = "";
	if ($inputs[1] eq "") {
		$input = $inputs[0];
		$input =~ s/\r/\n/gs;
	} else {
		foreach my $line (@inputs) {
			$input .= "$line";
		}
	}
	close (fileIN);

	#remove comment blocks
	$input =~ s/\[.*?\]//sg;

	#find blocks:
	my $blocks = {};
	while ($input =~ /Begin (.+?);(.+?)End;/isg) {
		my $blockname = uc($1);
		my $block = $2;
		$blocks->{$blockname} = $block;
	}

	my @blocks_to_process = qw (TAXA CHARACTERS DATA);
	foreach my $blockname (@blocks_to_process) {
		if (exists $blocks->{$blockname}) {
			my $statements = {};
			while ($blocks->{$blockname} =~ /\s*(.+?)\s+(.+?);/isg) {
				my $statementname = uc($1);
				my $statement = $2;
				$statements->{$statementname} = $statement;
			}
			$blocks->{$blockname} = $statements;
		}
	}

	# look for TAXA block:
	if (exists $blocks->{"TAXA"}) {
		@{$nexushash->{"taxa"}} = split (/\s+/, $blocks->{"TAXA"}->{"TAXLABELS"});
		$nexushash->{"ntax"} = scalar @{$nexushash->{"taxa"}};
		Subfunctions::debug ("processing TAXA block, found $nexushash->{ntax} taxa\n");
	}

	# look for CHARACTERS or DATA block:
	my $datablock;
	if (exists $blocks->{"CHARACTERS"}) {
		$datablock = $blocks->{"CHARACTERS"};
	} elsif (exists $blocks->{"DATA"}) {
		$datablock = $blocks->{"DATA"};
	}

	if ($datablock) {
		$nexushash->{"ntax"} = scalar @{$nexushash->{"taxa"}};
		if (exists $datablock->{"FORMAT"}) {
			my @params = ();
			my $paramline = "$datablock->{'FORMAT'}";
			while ($paramline =~ /(.+?)\s*=\s*(\S+)\s*(.*)/) {
				push @params, (uc($1)."=".uc($2));
				$paramline = $3;
			}
			foreach my $param (@params) {
				if ($param =~ /GAP=(.*)/) {
					$gapchar = $1;
				} elsif ($param =~ /INTERLEAVE=(.*)/) {
					if ($1 =~ /YES/) {
						$interleave = 1;
					}
				}
			}
		}
		if (exists $datablock->{"DIMENSIONS"}) {
			my @params = ();
			my $paramline = "$datablock->{'DIMENSIONS'}";
			while ($paramline =~ /(.+?)\s*=\s*(\S+)\s*(.*)/) {
				push @params, (uc($1)."=".uc($2));
				$paramline = $3;
			}
			foreach my $param (@params) {
				if ($param =~ /NTAX=(.*)/) {
					if (($nexushash->{"ntax"} > 0) && ($1 != $nexushash->{"ntax"})) {
						die "ntax in dimensions does not match ntax in taxa block.\n";
					}
					$nexushash->{"ntax"} = $1;
				} elsif ($param =~ /NCHAR=(.*)/) {
					$nchar = $1;
				}
			}
		}
		if (exists $datablock->{"MATRIX"}) {
			my $taxa = $nexushash->{"characters"};
			my $matrix = "$datablock->{'MATRIX'}";
			my @lines = split (/\n/, $matrix);
			foreach my $line (@lines) {
				if ($line =~ /\s*(\S+)\s+(\S+)/) {
					my $currtaxon = $1;
					my $currdata = $2;
					$currtaxon =~ s/\'//g;
					$currtaxon =~ s/\"//g;
					$currdata =~ s/$gapchar/-/g;
					if (exists $nexushash->{"characters"}->{$currtaxon}) {
						$nexushash->{"characters"}->{$currtaxon} .= $currdata;
					} else {
						$nexushash->{"characters"}->{$currtaxon} = $currdata;
					}
				}
			}
			foreach my $taxon (keys %$taxa) {
				if (length($taxa->{$taxon}) != $nchar) {
					die "Characters specified for $taxon do not match $nchar.";
				}
			}
		}
	}

	$nexushash->{"nchar"} = $nchar;
	Subfunctions::debug (Dumper($nexushash));
	return $nexushash;
}

sub clean_nexus_taxa_names {
	my $taxarray = shift;

	my $cleanedtaxarray = ();
	foreach my $taxon (@$taxarray) {
		my $cleanedtaxon = $taxon =~ s/-/_/gr;
		push @$cleanedtaxarray, $cleanedtaxon;
	}

	return $cleanedtaxarray;
}

sub write_nexus_character_block {
	my $nexushash = shift;
	my $blocksize=2000;
	my $nexblock = "";
	my $nchar = 1;
	my $i = 1;
	my $flag = 1;
	my $len;

	unless (exists $nexushash->{"characters"}) {
		print "write_nexus_character_block: no characters specified.\n";
		exit;
	}

	unless (exists $nexushash->{"taxa"}) {
		$nexushash->{"taxa"} = (keys %{$nexushash->{"characters"}});
	}

	my $cleanedtaxarray = clean_nexus_taxa_names($nexushash->{"taxa"});

	#copy working versions
	my @working_seqs = ();
	foreach my $t (@{$nexushash->{"taxa"}}) {
		push @working_seqs, "$nexushash->{characters}->{$t}";
	}

	Subfunctions::pad_seq_ends (\@working_seqs, "-");

	my $nchar = length($working_seqs[0]);
	while ((length $working_seqs[0]) >= $blocksize) {
		for (my $i=0; $i < @$cleanedtaxarray; $i++) {
			$working_seqs[$i] =~ /^(.{$blocksize})(.*)$/;
			$working_seqs[$i] = $2;
			$nexblock .= "" . @$cleanedtaxarray[$i] . "\t";
			$nexblock .= "$1\n";
		}
	}
	for (my $i=0; $i < @$cleanedtaxarray; $i++) {
		$nexblock .= "" . @$cleanedtaxarray[$i] . "\t";
		$nexblock .= "$working_seqs[$i]\n";
	}

	my $charstatelabels = "";
	if (exists $nexushash->{"charlabels"}) {
		my @charstate_arr = ();
		my $taxa_names = $nexushash->{"charlabels"};
		for (my $i=1; $i<= @{$nexushash->{"charlabels"}}; $i++) {
			my $taxon = @{$nexushash->{"charlabels"}}[$i-1];
			push @charstate_arr, "$i $taxon";
		}
		$charstatelabels = "CharStateLabels\n" . join (", ", @charstate_arr) . ";\n";
	}

	my $result = "Begin CHARACTERS;\nDimensions nchar=$nchar;\n";
	$result .= "$charstatelabels";
	$result .= "Format datatype=dna gap=- interleave=yes;\n";
 	$result .= "Matrix\n$nexblock;\nEnd;\n\n";
	return $result;
}

sub write_nexus_taxa_block {
	my $nexushash = shift;

	unless (exists $nexushash->{"taxa"}) {
		print "write_nexus_taxa_block: no taxa specified.\n";
		exit;
	}

	my $cleanedtaxarray = clean_nexus_taxa_names($nexushash->{"taxa"});

	my $taxlabels = "";
	for (my $i=1; $i<= @$cleanedtaxarray; $i++) {
		$taxlabels .= "[$i " . @$cleanedtaxarray[$i-1] . "]\n";
	}
	$taxlabels .= "TaxLabels " . join(" ", @$cleanedtaxarray) . ";";

	my $result = "begin TAXA;\n";
	$result .= "Dimensions ntax=" . @$cleanedtaxarray . ";\n";
	$result .= "$taxlabels\n";
	$result .= "End;\n\n";
	return $result;
}

sub write_nexus_sets_block {
	my $nexushash = shift;

	unless (exists $nexushash->{"sets"}) {
		print "write_nexus_taxa_block: no taxa specified.\n";
		exit;
	}

	my $setblock = "";
	foreach my $settype (keys %{$nexushash->{"sets"}}) {
		my $sets = $nexushash->{"sets"}->{$settype};
		foreach my $set (@$sets) {
			$setblock .= "$settype " . $set->{"name"} . " = " . $set->{"start"} . "-" . $set->{"end"} . ";\n";
		}
	}

	my $result = "begin SETS;\n";
	$result .= "$setblock";
	$result .= "End;\n\n";
	return $result;
}



sub write_nexus_trees_block {
	my $nexushash = shift;
	my $tree_array = $nexushash->{"trees"}; # an array of trees

	unless (exists $nexushash->{"trees"}) {
		print "write_nexus_trees_block: no trees specified.\n";
		exit;
	}

	my @name_blocks = ();
	if (exists $nexushash->{"taxa"}) {
		# push the one primary name block
		push @name_blocks, $nexushash->{"taxa"};
	}
	if (exists $nexushash->{"alttaxa"}) {
		# push all of the possible alternate name blocks
		push @name_blocks, @{$nexushash->{"alttaxa"}};
	}

	my $treeblock = "";
	foreach my $tree (@$tree_array) {
		my @keys = keys %$tree;
		my $k = shift @keys;
		$treeblock .= "Tree $k = $tree->{$k}";
		if ($treeblock !~ /;$/) {
			$treeblock .= ";";
		}
	}

	my $translate = "";
	if (@name_blocks > 0) {
		my @trans_arr = ();
		my $taxa_names = $nexushash->{"taxa"};
		for (my $i=1; $i<= @$taxa_names; $i++) {
			my $taxon = @$taxa_names[$i-1];
			if (!(exists $trans_arr[$i-1])) {
				push @trans_arr, "$i $taxon";
			}
		}

		foreach my $taxa_names (@name_blocks) {
			for (my $i=1; $i<= @$taxa_names; $i++) {
				my $taxon = @$taxa_names[$i-1];
				$treeblock =~ s/$taxon/$i/g;
			}
		}
		$translate = "Translate\n" . join (",\n", @trans_arr) . ";\n";
	}

	my $result = "begin TREES;\n";
	$result .= "$translate$treeblock";
	$result .= "\nEnd;\n\n";
	return $result;
}

# must return 1 for the file overall.
1;
