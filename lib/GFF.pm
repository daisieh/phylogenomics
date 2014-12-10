#!/usr/bin/perl

package GFF;
use strict;
use FindBin;
use Data::Dumper;
use lib "$FindBin::Bin/../";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code line_wrap subseq_from_fasta);


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(feature_to_seq parse_gff_file parse_gff_block parse_attributes export_gff_block read_gff_block write_gff_file set_gff_sequence);
}

=head1

GFF data structure:
$gff_array: an array of genes in a GFF file
$gff_gene_hash: for each gene, a hash
	->{"ID"}
	->{"source"}: 	where the annotation is from
	->{"score"}:
	->{'Name'}
	->{'start'}: 	gene start
	->{'end'}:		gene end
	->{'seqid'}:	scaffold/LG name
	->{'phase'}
	->{'strand'}
	->{'type'}:		gene, CDS, etc
	->{'attributes'}
	->{"mRNA"}: 	an array of hashes, keyed by number, each representing a transcript
		 'ID' => 'Potri.001G000300.1.v3.0',
		 'score' => '.',
		 'CDS' => {
					'1' => {
							 'phase' => '0',
							 'strand' => '+',
							 'score' => '.',
							 'attributes' => {
											   'pacid' => '27045442',
											   'Parent' => 'Potri.001G000300.1.v3.0'
											 },
							 'end' => 228,
							 'start' => 19
						   }
				  },
		 'end' => 470,
		 'phase' => '.',
		 'three_prime_UTR' => {
								'1' => {
										 'phase' => '.',
										 'strand' => '+',
										 'score' => '.',
										 'attributes' => {
														   'pacid' => '27045442',
														   'Parent' => 'Potri.001G000300.1.v3.0'
														 },
										 'end' => 470,
										 'start' => 229
									   }
							  },
		 'strand' => '+',
		 'five_prime_UTR' => {
							   '1' => {
										'phase' => '.',
										'strand' => '+',
										'score' => '.',
										'attributes' => {
														  'pacid' => '27045442',
														  'Parent' => 'Potri.001G000300.1.v3.0'
														},
										'end' => 18,
										'start' => 1
									  }
							 },
		 'attributes' => {
						   'pacid' => '27045442',
						   'longest' => '1',
						   'Parent' => 'Potri.001G000300.v3.0'
						 },
		 'Name' => 'Potri.001G000300.1',
		 'start' => 1
	   }
},

=cut


sub parse_gff_file {
	my $gff_file = shift;

	my $gff_array = ();
	open my $gff_fh, "<", $gff_file;
	my $line = "";
	while ($line = readline $gff_fh) {
		if ($line =~ /^##/) {
			next;
		}
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $line);
		if ($type eq "gene") {
			$attributes =~ /ID=(.+?);/;
			my $gene = $1;
			seek($gff_fh, -length($line), 1);
			my $gff_block = parse_gff_block (read_gff_block ($gff_fh, $gene));

			push @$gff_array, $gff_block;
		}
	}
	close $gff_fh;
	return $gff_array;
}

sub read_gff_block {
	my $gff_fh = shift;
	my $gene = shift;

	my $gff_block = "";
	my $in_gene = 0;
	if (!(defined $gene)) {
		my $line = readline $gff_fh;
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $line);
		$attributes =~ /ID=(.+?);/;
		$gene = $1;
		seek($gff_fh, -length($line), 1);
	}
	while (my $line = readline $gff_fh) {
		if ($line =~ /##gff-version 3/) {
			next;
		}
		if ($line =~ /^##/) {
			seek($gff_fh, -length($line), 1);
			last;
		}
		# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
		if ($line =~ /gene.*$gene/) {
			$gff_block = $line;
			$in_gene = 1;
		} elsif (($line =~ /gene/) && ($in_gene == 1)) {
			seek($gff_fh, -length($line), 1);
			last;
		} elsif ($in_gene == 1) {
			$gff_block .= $line;
		}
	}
	return $gff_block;
}

sub feature_to_seq {
	my $sequence = shift;
	my $feature = shift;
	my $params = shift;

	my ($padded, $separate) = 0;
	if ((ref $params) =~ /HASH/) {
		if (defined $params->{"padded"}) {
			$padded = $params->{"padded"};
		}
		if (defined $params->{"separate"}) {
			$separate = $params->{"separate"};
		}
		if (defined $params->{"offset"}) {

		}

	}

	my $finalseq = "";
	my @seqarray = ();
	if (! (defined $feature)) {
		return "";
	}
	for (my $i = 1; $i <= keys %$feature; $i++) {
		my $feat = $feature->{$i};
		my ($startseq, $seq, $endseq) = split_seq ($sequence, $feat->{"start"}, $feat->{"end"});
		if ($padded == 1) {
			push @seqarray, "x" x length($startseq) . $seq . "x" x length($endseq);
		} else {
			push @seqarray, $seq;
		}
	}

	if ($separate == 0) {
		while ($seqarray[0] ne "") {
			my $currchars = "";
			foreach my $seq (@seqarray) {
				if ($seq =~ m/(.)(.*)$/) {
					$currchars .= $1;
					$seq = $2;
				}
			}
			$currchars =~ s/x//g;
			if ($currchars eq "") {
				$finalseq .= "N";
			} else {
				$finalseq .= get_iupac_code ($currchars);
			}
		}
		return \($finalseq);
	} else {
		return \@seqarray;
	}


}

sub set_gff_sequence {
	my $gff_hash = shift;
	my $sequence = shift;

	$gff_hash->{"sequence"} = $sequence;
	$gff_hash->{"start"} = 1;
	$gff_hash->{"end"} = length $gff_hash->{"sequence"};
}

sub parse_gff_block {
	# must pass in a block corresponding to a single gene.
	# should create a dictionary where the "gene" key contains the actual full location info
	# and all subsequent type entries contain offset information.
	my $gff_block = shift;

	my @lines = split(/\n/, $gff_block);
	my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = 0;
	my $offset = 0;
	my $gff_hash = {};

	my $gene_line = shift @lines;
	# this first line is the gene line. We are going to need to use this to get the whole sequence and the offset info.
	# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
	($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $gene_line);
	if ($type =~ /gene/) {
		my $attr_hash = parse_attributes($attributes);
		$offset = $start - 1;
		$gff_hash->{"seqid"} = $seqid;
		$gff_hash->{"source"} = $source;
		$gff_hash->{"type"} = $type;
		$gff_hash->{"start"} = $start;
		$gff_hash->{"end"} = $end;
		$gff_hash->{"score"} = $score;
		$gff_hash->{"strand"} = $strand;
		$gff_hash->{"phase"} = $phase;
		$gff_hash->{"ID"} = delete $attr_hash->{"ID"};
		$gff_hash->{"Name"} = delete $attr_hash->{"Name"};
		$gff_hash->{"attributes"} = $attr_hash;
	}
	my $gene_name = $gff_hash->{"Name"};

	# read in the block and hash it:
	foreach my $line (@lines) {
		if ($line =~ /^\s*$/) {
			next;
		}
		($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $line);
		# the first two are the same as its parent gene, so can be discarded.
		# type is the key under which this should be filed, with a subkey for the unique id of this feature.
		#	the id is named after the parent: #PAC:27046906.exon.3
		# the start and end are offsets from the parent gene.
		# we don't care about the next three, usually, but we might as well hash them too.
		# keep the attributes as is.

		$start = $start - $offset;
		$end = $end - $offset;

		my $attr_hash = parse_attributes($attributes);
		my $id = delete $attr_hash->{"ID"};
		my $name = delete $attr_hash->{"Name"};

		my $hash_ptr = {};
		$hash_ptr->{"start"} = $start;
		$hash_ptr->{"end"} = $end;
		$hash_ptr->{"score"} = $score;
		$hash_ptr->{"strand"} = $strand;
		$hash_ptr->{"phase"} = $phase;
		$hash_ptr->{"attributes"} = $attr_hash;


		if ($type eq "mRNA") {
			# Chr01	phytozome9_0	mRNA	8391	12209	.	-	.	ID=PAC:27046907;Name=Potri.001G000400.3;pacid=27046907;longest=0;Parent=Potri.001G000400
			if ($name =~ /$gene_name\.(\d+)/) {
				$gff_hash->{"mRNA"}->{$id} = $hash_ptr;
				$hash_ptr->{"ID"} = $id;
				$hash_ptr->{"Name"} = $name;
			}
		} else {
			# Chr01	phytozome9_0	CDS	11082	11166	.	-	0	ID=PAC:27046907.CDS.1;Parent=PAC:27046907;pacid=27046907
			if ($id =~ /(.*?)\.$type\.(\d+)/) {
				$gff_hash->{"mRNA"}->{$1}->{$type}->{$2} = $hash_ptr;
			}
		}
	}
	# now that we've finished hashing, we can rename the mRNAs with the numerical index.
	my $mRNA_hash = {};
	foreach my $k (keys %{$gff_hash->{"mRNA"}}) {
		my $name = $gff_hash->{"mRNA"}->{$k}->{"Name"};
		my $this_hash = delete $gff_hash->{"mRNA"}->{$k};
		$name =~ /$gene_name\.(\d+)/;
		$mRNA_hash->{$1} = $this_hash;
	}
	$gff_hash->{"mRNA"} = $mRNA_hash;

	return $gff_hash;
}

sub parse_attributes {
	my $attributes = shift;
	# ID=Potri.019G067600;Name=Potri.019G067600
	# ID=PAC:27027285;Name=Potri.T085300.1;pacid=27027285;longest=1;Parent=Potri.T085300
	my $attr_hash = {};
	my @attrs = split (/;/, $attributes);
	foreach my $attr (@attrs) {
		if ($attr =~ /(.*?)=(.*)/) {
			$attr_hash->{$1} = $2;
		}
	}
	return $attr_hash;
}

sub export_attributes {
	my $attr_hash = shift;

	my $attributes = "";
	foreach my $k (keys %$attr_hash) {
		$attributes .= "$k=$attr_hash->{$k};";
	}
	$attributes =~ s/;$//;
	return $attributes;
}

sub export_gff_block {
	my $gff_hash = shift;

	my $gff_string = "";
	my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = 0;

	# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
	$seqid = $gff_hash->{"seqid"};
	$source = $gff_hash->{"source"};
	$type = $gff_hash->{"type"};
	$start = $gff_hash->{"start"};
	$end = $gff_hash->{"end"};
	$score = $gff_hash->{"score"};
	$strand = $gff_hash->{"strand"};
	$phase = $gff_hash->{"phase"};
	$attributes = "ID=$gff_hash->{ID};Name=$gff_hash->{Name};";
	$attributes .= export_attributes ($gff_hash->{"attributes"});
	my @line = ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes);
	$gff_string = join("\t",@line) . "\n";

	for (my $i=1; exists $gff_hash->{"mRNA"}->{$i}; $i++) {
		$type = "mRNA";
		$start = $gff_hash->{"mRNA"}->{$i}->{"start"};
		$end = $gff_hash->{"mRNA"}->{$i}->{"end"};
		my $id = $gff_hash->{"mRNA"}->{$i}->{"ID"};
		my $name = $gff_hash->{"mRNA"}->{$i}->{"Name"};
		$attributes = "ID=$id;Name=$name;";
		$attributes .= export_attributes ($gff_hash->{"mRNA"}->{$i}->{"attributes"});
		my @line = ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes);
		$gff_string .= join("\t",@line) . "\n";
		foreach my $type (keys %{$gff_hash->{"mRNA"}->{$i}}) {
			my $mRNA_hash = $gff_hash->{"mRNA"}->{$i};
			if ((ref $mRNA_hash->{$type}) =~ /HASH/) {
				for (my $j=1; exists $mRNA_hash->{$type}->{$j}; $j++) {
				# Chr01	phytozome9_0	CDS	11082	11166	.	-	0	ID=PAC:27046907.CDS.1;Parent=PAC:27046907;pacid=27046907
					$start = $mRNA_hash->{$type}->{$j}->{"start"};
					$end = $mRNA_hash->{$type}->{$j}->{"end"};
					$attributes = "ID=".$id.".".$type.".".$j.";";
					$attributes .= export_attributes ($mRNA_hash->{$type}->{$j}->{"attributes"});
					my @line = ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes);
					$gff_string .= join("\t",@line) . "\n";
				}
			}
		}
	}

	return $gff_string;
}

sub write_gff_file {
	my $gff_hash = shift;
	my $outfile = shift;

	open OUTFH, ">", $outfile;
	print OUTFH "##gff-version 3\n";
	print OUTFH export_gff_block ($gff_hash);
	if (exists $gff_hash->{"sequence"}) {
		print OUTFH "##FASTA\n";
		my $seq = line_wrap ($gff_hash->{"sequence"}, 80);
		print OUTFH ">$gff_hash->{ID}\n$seq\n";
	}
	close OUTFH;
}

# a non-destructive way to walk through a gff_hash
# 	for (my $i=1; exists $gff_hash->{"mRNA"}->{$i}; $i++) {
# 		print "mRNA $i\n";
# 		foreach my $type (keys $gff_hash->{"mRNA"}->{$i}) {
# 			my $mRNA_hash = $gff_hash->{"mRNA"}->{$i};
# 			if ((ref $mRNA_hash->{$type}) =~ /HASH/) {
# 				for (my $j=1; exists $mRNA_hash->{$type}->{$j}; $j++) {
# 					print "type $type $j\n";
# 				}
# 			}
# 		}
# 	}


return 1;
