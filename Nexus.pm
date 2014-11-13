package Nexus;
use strict;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
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
	our @EXPORT_OK   = qw(parse_nexus write_nexus_character_block write_nexus_taxa_block write_nexus_trees_block);
}


=head1

B<(\%taxa, \@taxanames) parse_nexus ( String $inputfile )>

Given a NEXUS file as input, returns a hash containing all the sequences, keyed by the
values of the taxanames array.

$inputfile:   NEXUS file to parse.

=cut

sub parse_nexus {
	my $inputfile = shift;

	my $taxa = {};
	my @taxonlabels = ();
	my $gapchar = "-";
	my $interleave = 0;
	my $nchar = 0;
	my $ntax = 0;

	open fileIN, "<", "$inputfile" or die "no file named $inputfile";
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
		@taxonlabels = split (/\s+/, $blocks->{"TAXA"}->{"TAXLABELS"});
	}

	# look for CHARACTERS or DATA block:
	my $datablock;
	if (exists $blocks->{"CHARACTERS"}) {
		$datablock = $blocks->{"CHARACTERS"};
	} elsif (exists $blocks->{"DATA"}) {
		$datablock = $blocks->{"DATA"};
	}

	if ($datablock) {
		$ntax = scalar @taxonlabels;
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
					if (($ntax > 0) && ($1 != $ntax)) {
						die "ntax in dimensions does not match ntax in taxa block.\n";
					}
					$ntax = $1;
				} elsif ($param =~ /NCHAR=(.*)/) {
					$nchar = $1;
				}
			}
		}
		if (exists $datablock->{"MATRIX"}) {
			my $matrix = "$datablock->{'MATRIX'}";
			my @lines = split (/\n/, $matrix);
			foreach my $line (@lines) {
				if ($line =~ /\s*(\S+)\s+(\S+)/) {
					my $currtaxon = $1;
					my $currdata = $2;
					$currtaxon =~ s/\'//g;
					$currtaxon =~ s/\"//g;
					$currdata =~ s/$gapchar/-/g;
					if (exists $taxa->{$currtaxon}) {
						$taxa->{$currtaxon} = $taxa->{$currtaxon} . $currdata;
					} else {
						$taxa->{$currtaxon} = $currdata;
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

	$taxa->{"length"} = $nchar;

	return $taxa, \@taxonlabels;
}

=head1

B<String $nexus_str convert_aln_to_nexus ( SimpleAlign $aln )>

Takes a taxon hash (and optional taxon array) and returns a NEXUS-formatted string representing the alignment.

=cut

sub write_nexus_character_block {
	my $taxa_hash = shift;
	my $taxarray = shift;
	my $blocksize=2000;
	my $nexblock = "";
	my $result = "";
	my $nchar = 1;
	my $i = 1;
	my $flag = 1;
	my $len;

	if ($taxarray == undef) {
		$taxarray = (keys %$taxa_hash);
	}

	#copy working versions
	my @working_seqs = ();
	foreach my $t (@$taxarray) {
		push @working_seqs, "$taxa_hash->{$t}";
	}

	pad_seq_ends (\@working_seqs, "-");

	my $nchar = length($working_seqs[0]);
	while ((length $working_seqs[0]) >= $blocksize) {
		for (my $i=0; $i < @$taxarray; $i++) {
			$working_seqs[$i] =~ /^(.{$blocksize})(.*)$/;
			$working_seqs[$i] = $2;
			$nexblock .= "" . @$taxarray[$i] . "\t";
			$nexblock .= "$1\n";
		}
	}
	for (my $i=0; $i < @$taxarray; $i++) {
		$nexblock .= "" . @$taxarray[$i] . "\t";
		$nexblock .= "$working_seqs[$i]\n";
	}

	$result .= "Begin CHARACTERS;\nDimensions nchar=$nchar;\n";
	$result .= "Format datatype=dna gap=- interleave=yes;\n";
 	$result .= "Matrix\n$nexblock;\nEnd;\n\n";
	return $result;
}

sub write_nexus_taxa_block {
	my $taxa_names = shift;

	my $result = "begin TAXA;\n";
	$result .= "Dimensions ntax=" . @$taxa_names . ";\n";
	my $taxlabels = "";
	for (my $i=1; $i<= @$taxa_names; $i++) {
		$result .= "[$i " . @$taxa_names[$i-1] . "]\n";
	}
	$result .= "TaxLabels " . join(" ", @$taxa_names) . ";\nEnd;\n\n";

	return $result;
}

sub write_nexus_trees_block {
	my $trees = shift;

	my @name_blocks = ();
	while (my $block = shift) {
		push @name_blocks, $block;
	}
	print Dumper(@name_blocks);
#  	my $taxa_names = shift @name_blocks;
	my $result = "begin TREES;\n";

	my $treeblock = "";
	foreach my $k (keys %$trees) {
		$treeblock .= "Tree $k = $trees->{$k}";
		if ($treeblock !~ /;$/) {
			$treeblock .= ";";
		}
		$treeblock .= "\n";
	}

	my $translate = "";
	my @trans_arr = ();

	foreach my $taxa_names (@name_blocks) {
		for (my $i=1; $i<= @$taxa_names; $i++) {
			my $taxon = @$taxa_names[$i-1];
			$treeblock =~ s/$taxon/$i/g;
			if (!(exists $trans_arr[$i-1])) {
				push @trans_arr, "$i $taxon";
			} else {
				$trans_arr[$i-1] .= " $taxon";
			}
		}
	}
	$translate = "Translate\n" . join (",\n", @trans_arr) . ";\n";
	$result .= "$translate$treeblock";
	$result .= "\nEnd;\n\n";
	return $result;
}



# must return 1 for the file overall.
1;
