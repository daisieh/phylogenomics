#!/usr/bin/perl

use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Bioperl qw (get_locations_from_genbank_file);
use CircleGraph;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($samplefile, $gb_file, $labelfile, $outfile, $label_samples) = 0;
my $samplesize = 1000;
my $keepfiles = 0;
my $help = 0;
my $min_coverage = 0;
my $circle_size = 0;


GetOptions ('samples|input=s' => \$samplefile,
            'outputfile=s' => \$outfile,
            'window:i' => \$samplesize,
            'labels:s' => \$labelfile,
            'genbank|gb:s' => \$gb_file,
            'keepfiles!' => \$keepfiles,
            'coverage:i' => \$min_coverage,
            'size:i' => \$circle_size,
            'legend' => \$label_samples,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my @names;
if ($samplefile =~ /(.+?)\.depth/) {
    $samplefile = $1;
    push @names, $samplefile;
} else {
    open FH, "<", "$samplefile" or die "no such file";
    my @rawnames = <FH>;
    close FH;

    # don't process blank names
    foreach my $name (@rawnames) {
        $name =~ s/\n//;
        if ($name !~ /^\s*$/) {
            my $samplename = basename($name);
            if ($samplename =~ /(.*?)\.depth/) {
            	$samplename = $1;
            }
            push @names, $samplename;
        }
    }
}


# process the gb file
if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
        my ($fh, $filename) = tempfile(UNLINK => 1);
        my $gb_locs = get_locations_from_genbank_file($gb_file);
		print $fh $gb_locs;
		close $fh;
	    $gb_file = "$filename";

	    $gb_locs =~ /(.+)\n(.+?)$/;
	    my $source = $2;
	    $source =~ /.+?\t.+?\t(.+)$/;
	    $circle_size = $1;
	}
} else {
	if ($circle_size == 0) {
		die "No genbank file was provided: must specify circle size in base pairs.\n";
	}
}


my $labels;
if ($labelfile) {
    $labels = make_label_lookup ($labelfile);
}

my %samples;
my @sample_list;
foreach my $samplename (@names) {
    my $depthfile = "$samplename".".depth";
    $samplename = basename($samplename);
    print "processing $samplename\n";
    open VCF_FH, "<", "$depthfile" or die "couldn't open $depthfile";
    my @depths = ();
    my @indels = ();
    my $i = 1;
    my $in_indel = 0;
    my $line = readline VCF_FH;
    my $curr_indel_start = 0;

    #eat any header lines:
    while ($line =~ /#/) {
    	$line = readline VCF_FH;
    }

	$line =~ /\t(.*?)\t(.*?)\s/;
	if ($2 == 0) {
		$curr_indel_start = 1;
	}
    for ($i=1; $i <= $circle_size; $i++) {
        $line =~ /\t(.*?)\t(.*?)\s/;
        my $depth = $2;
        if ($i != $1) { # this position isn't listed in the depths file; assume it is 0
        	$depth = 0;
        } else {
			$line = readline VCF_FH;
        }
		if ($curr_indel_start) { # are we in an indel already?
			if ($depth > $min_coverage) { # we're not in an indel anymore; end the indel.
				push @indels, "$curr_indel_start\t$i";
				$curr_indel_start = 0;
			}
		} else { # we're not in an indel. Start one if the depth is 0.
			if ($depth <= $min_coverage) {
				$curr_indel_start = $i;
			}
		}
		push @depths, $depth;
    }
    close VCF_FH;
    if ($curr_indel_start) {
        push @indels, "$curr_indel_start\t".($i-1);
    }
    $samples{$samplename}{"depths"} = \@depths;
    $samples{$samplename}{"indels"} = \@indels;
}


while (my ($key, $value) = each %samples) {
    my @depth_array = @{$value->{"depths"}};
    open FH, ">", "$outfile"."_$key"."_depths.txt";
    print FH "pos\t$key\n";
    my $sum = 0;
    my $i;
    for ($i=1; $i<=@depth_array; $i++) {
        if ($samplesize) {
            if (($i%$samplesize)==0) {
                print FH "$i\t".@depth_array[$i-1]."\n";
            }
        } else {
            print FH "$i\t".@depth_array[$i-1]."\n";
        }
    }
    close FH;

    my @indel_array = @{$value->{"indels"}};
    open FH, ">", "$outfile"."_$key"."_indels.txt";
    my $sum = 0;
    my $i;
    for ($i=1; $i<=@indel_array; $i++) {
        print FH "$i\t".@indel_array[$i-1]."\n";
    }
    print FH "size\t1\t$circle_size\n";
    close FH;
}

print "drawing graphs...\n";
my $circlegraph_obj = new CircleGraph();
$circlegraph_obj->inner_radius($circlegraph_obj->inner_radius - 100);
my $max=1;
my %graphs;
while (my ($key, $value) = each %samples) {
    # process each depths file
    open FH, "<", "$outfile"."_$key"."_depths.txt";
    my @items = <FH>;
    close FH;
    my $header = shift @items;
    my (@xvals, @yvals);
    foreach my $line (@items) {
        $line =~ /(.*?)\t(.*?)\s/;
        push @xvals, $1;
        push @yvals, $2;
        if ($2 > $max) {
            $max = $2;
        }
    }
    if (!exists($graphs{"x"})) {
        $graphs{"x"} = \@xvals;
    }
    $graphs{$key} = \@yvals;
}

#process yvals so that they are scaled to 1
while (my ($key, $value) = each %samples) {
    my @yvals = @{$graphs{$key}};
    for (my $y=0;$y<@yvals;$y++) {
        @{$graphs{$key}}[$y] = @yvals[$y]/$max;
    }
}

# draw the gene map
if ($gb_file) {
	draw_gene_map ($gb_file, $circlegraph_obj, {direction=>"OUT",width=>15});
} else {
	$circlegraph_obj->set_font("Helvetica", 6, "black");
	my $radius = $circlegraph_obj->outer_radius + 10;
	# whatever the size, probably should draw about 100 labels around the perimeter.
	my $ideal_interval = $circle_size/100;
	my $factor = 5;
	for (my $i=1;$i < (log($ideal_interval)/log(10)); $i++) {
		$factor = $factor * 10;
	}
	my $interval = 0;
	while ($interval < $ideal_interval) {
		$interval += $factor;
	}
	for (my $i = 0; $i < $circle_size; $i=$i+($interval/5)) {
		my $angle = ($i/$circle_size) * 360;
		if (($i % $interval) == 0) {
			$circlegraph_obj->circle_label($angle, $radius, "$i");
		} else {
			$circlegraph_obj->circle_label($angle, $radius, "-");
		}
	}
}

# if the minimum coverage is specified as > 0, draw a threshold bar for that.
if ($min_coverage > 0) {
	my $threshold = (($min_coverage/$max)*($circlegraph_obj->outer_radius - $circlegraph_obj->inner_radius)) + $circlegraph_obj->inner_radius;
	$circlegraph_obj->draw_circle($threshold, {filled=>1, color=>"light_grey"});
	$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius, {filled=>1, color=>"white"});
}

# draw the maps
my $j = 0;
my @xvals = @{$graphs{"x"}};
for (my $j=0; $j< keys(%samples); $j++) {
	my $key = $names[$j];
    my @yvals = @{$graphs{$key}};
    # plot the coverage map
    $circlegraph_obj->plot_line(\@xvals, \@yvals, {color=>$j});

    # plot the indel map
    $circlegraph_obj = draw_regions ( "$outfile"."_$key"."_indels.txt", $circlegraph_obj, $j, $circlegraph_obj->inner_radius + ($j*3) -150);
}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius);
$circlegraph_obj->draw_circle($circlegraph_obj->outer_radius);
$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius-155);

# draw the labels
if ($label_samples) {
	for ($j=0; $j< keys(%samples); $j++) {
		my $label = $names[$j];
		if (exists $labels->{$names[$j]}) {
			$label = $labels->{$names[$j]};
		}
		$circlegraph_obj->append_to_legend($label,$j);
	}
}

$circlegraph_obj->append_to_legend("Maximum coverage was $max, scaled to 1");
$circlegraph_obj->append_to_legend("Minimum coverage was $min_coverage");
$circlegraph_obj->append_to_legend("Sampling frequency was $samplesize");
$circlegraph_obj->draw_legend_text({size => 9, height => 12});

open FH, ">", "$outfile.ps" or die "couldn't open $outfile.ps";
print FH $circlegraph_obj->output_ps();
close FH;

if (!$keepfiles) {
    while (my ($key, $value) = each %samples) {
        unlink "$outfile"."_$key"."_depths.txt" or warn "could not delete $outfile"."_$key"."_depths.txt";
        unlink "$outfile"."_$key"."_indels.txt" or warn "could not delete $outfile"."_$key"."_indels.txt";
    }
}

__END__

=head1 NAME

analyze_depths

=head1 SYNOPSIS

analyze_depths -samples -window -output [-genbank] [-labels]

=head1 OPTIONS

  -samples|input:   input file of samples with depth files to analyze
  -window:          sampling frequency (if not specified, 1000 bp)
  -outputfile:      prefix of output files
  -genbank|gb:      optional: genbank file to generate a map along the graph
  -labels:          optional: tab-delimited list of samples and their corresponding labels
  -keepfiles:       optional: keep data files that are created (default is no)
  -coverage:		optional: minimum level of coverage to call as missing
  -legend:			optional: if present, include list of samples in the legend

=head1 DESCRIPTION

Graphs coverage depths using depth data generated from samtools mpileup.

=cut

