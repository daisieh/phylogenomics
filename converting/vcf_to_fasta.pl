#!/usr/bin/env perl

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(sample_list get_ordered_genotypes get_iupac_code);
use Data::Dumper;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $samplefile = 0;
my $help = 0;
my $cov_thresh = 5;
my $pl_thresh = 0;
my $outfile = "";
my $indels = 0;
my $multiple = 0;
my $start_pos = 0;
my $end_pos = 0;
my $range = "";

GetOptions ('samples|input|vcf=s' => \$samplefile,
            'outputfile=s' => \$outfile,
            'threshold|pl=i' => \$pl_thresh,
            'minimum|coverage|reads=i' => \$cov_thresh,
            'indels' => \$indels,
            'multiple' => \$multiple,
            'start=i' => \$start_pos,
            'end=i' => \$end_pos,
            'range=s' => \$range,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my $name;
my @samples = ();
my @samplefiles = ();
my %sample_positions = ();

#			GT	AD	DP	GQ	PL
#			0/0	20,0	20	60.2	0,60,816
#type:		s	i	i	i	i
#number:		1	.	1	1	G

if ($samplefile =~ /(.*?)\.vcf$/) {
	@samplefiles = ($samplefile);
} else {
	@samplefiles = @{sample_list ($samplefile)};
}


if ($samplefile =~ /recode\.vcf/) {
	$multiple = 1;
}

if ($range =~ /(\d+)\.\.(\d+)/) {
	$start_pos = $1;
	$end_pos = $2;
}

my $in_header = 1;
my $format_hash = {};
my $info_hash = {};
my @format_array = ();
my @unknown_array = ();

foreach my $samplefile (@samplefiles) {
	$in_header = 1;
	if ($samplefile =~ /\.vcf$/) {
		$samplefile .= ".vcf";
	}
	open VCF_FH, "<:crlf", $samplefile or die "couldn't open input file $samplefile.";
	my $i = 0;
	while (my $line = readline VCF_FH) {
		if ($in_header == 1) {
			# process header:
			if ($line =~ m/^##/) {
				if ($line =~ m/^##(FORMAT|INFO)=<ID=(.+?),Number=(.+?),Type=(.+?),Description=\"(.+)\">/) {
					my $prophash = {};
					$prophash->{'ID'} = $2;
					$prophash->{'Number'} = $3;
					$prophash->{'Type'} = $4;
					$prophash->{'Description'} = $5;
					$format_hash->{$prophash->{'ID'}} = $prophash;
					
				}
				next;
			}
	
			if ($line =~ m/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t*(FORMAT\t.+)$/) {
				@samples = split(/\t/,$1);
				if ((shift @samples) !~ /FORMAT/) {
					die "$samplefile does not have genotype data.\n";
				}
				$in_header = 0;
				print "found ".@samples." samples\n";
				foreach my $sample (@samples) {
					$sample_positions{$sample} = ();
				}
			}
		} else {
		
			# process sample lines
			#Chr19	15894518	.	G	A	113.75	.	.	GT:AD:DP:GQ:PL	0/0:20,0:20:60.2:0,60,816	0/0:12,0:12:36.11:0,36,454	0/0:15,0:15:45.09:0,45,555	0/0:10,0:10:30.1:0,30,417	0/0:17,0:17:51.11:0,51,642	0/1:13,3:16:52.21:52,0,434	./.:6,0:6:18.06:0,18,235	0/0:15,0:15:45.14:0,45,578	./.:5,0:5:15.05:0,15,207	0/0:19,0:19:57.03:0,57,634	0/0:10,0:10:30.04:0,30,352	0/0:14,0:15:39.09:0,39,483	0/0:16,0:16:48.15:0,48,621	0/0:14,0:14:42.14:0,42,552	0/0:24,0:24:69.21:0,69,878	0/0:11,0:11:33.11:0,33,437	0/0:14,0:14:42.07:0,42,502	0/0:33,0:33:99:0,99,1054	0/0:23,0:23:69.07:0,69,784	./.:7,0:7:21.07:0,21,285	0/0:21,0:21:63.17:0,63,804	0/0:13,0:13:39.07:0,39,471	0/1:20,4:24:61.58:62,0,652	0/1:18,5:23:74.47:74,0,623	0/0:19,0:19:57.19:0,57,772	0/0:27,0:27:78.1:0,78,904	0/0:22,0:22:66.21:0,66,861	0/0:12,0:12:36.1:0,36,443	0/0:12,0:12:36.05:0,36,415	0/0:11,0:11:30.09:0,30,370	0/0:37,0:37:99:0,111,1304	0/0:17,0:17:51.17:0,51,686	0/0:16,0:16:48.15:0,48,612	0/0:12,0:12:36.08:0,36,433	0/0:14,0:14:36.1:0,36,441	./.:8,0:8:24.08:0,24,328	0/0:19,0:19:57.15:0,57,734	0/0:28,0:28:84.21:0,84,1054	0/0:12,0:12:36.11:0,36,456
			my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$samples) = split (/\t/,$line, 10);

			# check if a positional range was specified:
			if ($pos < $start_pos) {
				next;
			}
			if (($end_pos > 0) && ($pos > $end_pos)) {
				last;
			}

			chomp $samples;
			my @samples_at_pos = split (/\t/,$samples);
			my $info_depth = 0;
			my $depth_pos = -1;
			my $pl_pos = -1;

			# for each position, I need to know ref allele, alt allele, read depth (DP, int), and phred score (PL, string)
			my @AMBIGUOUS_POS = ($ref,$alt,0,"");
			if ($i == 0) {
				$i = $pos;
			}
			
			# sometimes the depth is stored in the INFO field (for single sample files)
			if ($info =~ /.*DP=(.*?);.*/) {
				$info_depth = $1;
			}
			
			# other times the depth and PL are fields per sample: find out which from the FORMAT field
			@format_array = split (/:/,$format);
			
# 			genotype is always the first specified, so shift it off
# 			shift @format_array;
			
			for (my $j=0; $j<@format_array; $j++) {						
				if ($format_array[$j] eq "DP") {
					$depth_pos = $j;
				} elsif ($format_array[$j] eq "PL") {
					$pl_pos = $j;
				}
			}
			
			# if we're missing calls for a position, fill in ambiguities for all samples at this position
			while ($i < $pos) {
				for (my $j=0;$j<@samples_at_pos;$j++) {
					push @{$sample_positions{$samples[$j]}}, \@AMBIGUOUS_POS;
				}
				$i++;
			}

			# now process the calls for the specified position for all samples.
			for (my $j=0;$j<@samples_at_pos;$j++) {
				my $depth = 0;
				my $pl = "";
				my @sample_fields = split (/:/, $samples_at_pos[$j]);
				
# 				genotype is always the first specified, so shift it off
# 				shift @sample_fields;
				
				if (@sample_fields > 0) {
					if ($depth_pos != -1) {
						$depth = $sample_fields[$depth_pos];
					} else {
						$depth = $info_depth;
					}
					if ($pl_pos != -1) {
						$pl = $sample_fields[$pl_pos];
					}
				} 
				
				if ($depth == undef) { $depth = 0; }
				my @this_pos = ($ref, $alt, $depth, $pl);
				push @{$sample_positions{$samples[$j]}}, \@this_pos;
			}
			$i++;
		}
	}
}
print Dumper ($sample_positions) . "\n";
my $result_str = "";

if (@samples == 1) {
	unless ($outfile ne "") {
		$outfile = $name;
	}
}

my $fh;

if ($outfile eq "") {
	open $fh, ">&STDOUT";
} else {
	open $fh, ">", $outfile.".fasta";
}

foreach my $sample (@samples) {
	my $seq = "";
	my @positions = @{$sample_positions{$sample}};
# 	my @indels = @{$sample_positions{$name."_indels"}};
# 	if ($indels == 1) {
# 		foreach my $line (@indels) {
# 	# FRE13_HapA.plastome.3.aln.fasta	2	.	TGGG	TGGTTCATGGG	143	.	INDEL;DP=21;VDB=0.0000;AF1=1;AC1=2;DP4=0,0,11,0;MQ=36;FQ=-67.5	PL	184,33,0
# 			my ($chrom,$pos,undef,$ref,$alt,$qual,$filter,$info,undef,$pl) = split (/\t/,$line);
# 			chomp $pl;
# 			my @pls = split(/,/, $pl);
# 			my @alts = split(/,/, $alt);
# 			my $alleles = "0";
# 			for (my $i=1;$i<=@alts;$i++) {
# 				$alleles = $alleles.$i;
# 			}
# 			unshift @alts, $ref;
# 			my @genotypes = @{get_ordered_genotypes($alleles)};
# 			my $genotype = "";
# 			for (my $i=0;$i<@pls;$i++) {
# 				if ($pls[$i] == 0) {
# 					$genotype = $genotypes[$i];
# 					$genotype =~ /(.)(.)/;
# 					if ($1 != $2) {
# 						$genotype = "";
# 					} else {
# 						$genotype = $alts[$genotype];
# 					}
# 				}
# 			}
# 
# 			my ($snppos, $snpqual, $depth, $snpref, $snpalt, $snppl) = @{@positions[$pos-1]};
# 			print "$pos\t$depth\t$qual\t$snpqual\tindels: $genotype\n";
# 			if ($snppl =~ /0/) {
# 				$genotype = $snpref;
# 			} else {
# 				$genotype = "N";
# 			}
# 			my @pls = split(/,/, $snppl);
# 			my $alleles = "$ref$alt";
# 			$alleles =~ s/,//g;
# 			my @genotypes = @{get_ordered_genotypes($alleles)};
# 			my @newpos = ($pos, $qual, $depth, $genotype, "N", 0);
# 			@positions[$pos-1] = \@newpos;
# 		}
# 	}
	for (my $pos=0;$pos<@positions; $pos++) {
		my ($ref, $alt, $depth, $pl) = @{@positions[$pos]};
		if ($depth > $cov_thresh) {
			if ($pl =~ /,/) {
				my @pls = split(/,/, $pl);
				my $alleles = "$ref$alt";
				$alleles =~ s/,//g;
				my @alts = ($ref,$alt);
				my @genotypes = @{get_ordered_genotypes($alleles)};
				my $genotype = "";
				for (my $i=0;$i<@pls;$i++) {
					if ($pls[$i] == 0) {
						$genotypes[$i] =~ /(.)(.)/;
						if ($1 eq $2) {
							$genotype = $1;
						} else {
							$genotype = $genotypes[$i];
						}
					}
				}

				$alt = "N";
				my $genotype = "";
				my @sort_pls = sort @pls;
				my $smallest_pl = shift @sort_pls;
				for (my $g=0;$g<@pls;$g++) {
					if ($pls[$g] == $smallest_pl) {
						$genotype = $genotypes[$g];
						$genotype =~ s/-//g;
						$alt = get_iupac_code($genotype);
					}
				}
				if ($alt !~ /[AGCT]/) {
# 					print "$depth\t$genotype = $alt\n";
				}
				$seq .= $alt;
			} elsif ($pl <= $pl_thresh) {
				$seq .= $ref;
			} else {
				$seq .= "N";
			}
		} else {
			if ($indels != 1) {
				$seq .= "n";
			}
		}
	}
	print $fh ">$sample\n$seq\n";

}

close $fh;

__END__

=head1 NAME

vcf2fasta

=head1 SYNOPSIS

vcf2fasta -samplefile -output [-threshold]

=head1 OPTIONS

  -samples|input|vcf:   name of sample or list of samples to convert
  -outputfile:      optional: prefix of output fasta file
  -threshold|pl:	optional: maximum threshold for Phred-scaled likelihood (lower is best, default is 0)
  -min|coverage:    optional: read depth required to call base (default is 5)
  -multiple
            'multiple' => \$multiple,
            'start=i' => \$start_pos,
            'end=i' => \$end_pos,
            'range=s' => \$range,

=head1 DESCRIPTION

Converts a VCF file (or files) to a fasta file, using the Phred-scaled likelihood of each genotype call.
=cut

