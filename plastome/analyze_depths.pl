use CircleGraph;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";
require "circlegraphs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($samplefile, $gb_file, $labelfile, $outfile) = 0;
my $samplesize = 1000;
my $keepfiles = 0;
my $help = 0;

GetOptions ('samples|input=s' => \$samplefile,
            'outputfile=s' => \$outfile,
            'window:i' => \$samplesize,
            'labels:s' => \$labelfile,
            'genbank|gb:s' => \$gb_file,
            'keepfiles!' => \$keepfiles,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

open FH, "<", "$samplefile" or die "no such file";
my @names = <FH>;
close FH;

my $circle_size = 157033;

my $labels;
if ($labelfile) {
    $labels = make_label_lookup ($labelfile);
}

my %samples;
foreach my $name (@names) {
    $name =~ s/\n//;
    my $depthfile = "$name".".depth";
    print "processing $depthfile\n";
    open VCF_FH, "<", "$depthfile";
    my @depths = ();
    my @indels = ();
    my $i = 1;
    my $in_indel = 0;
    my $line = readline VCF_FH;
    my $curr_indel_start = 0;
	$line =~ /\t(.*?)\t(.*?)\s/;
	if ($2 == 0) {
		$curr_indel_start = 1;
	}
    for ($i=1; $i <= $circle_size; $i++) {
        $line =~ /\t(.*?)\t(.*?)\s/;
        if ($i == $1) {
            push @depths, $2;
            if ($2 == 0) {
				if ($curr_indel_start) {
					push @indels, "$curr_indel_start\t$i";
					$curr_indel_start = 0;
				} else {
					$curr_indel_start = $i;
				}
			}
            $line = readline VCF_FH;
        } else {
            push @depths, "0";
        }
#         $line = readline VCF_FH;
    }
    close VCF_FH;
    if ($curr_indel_start) {
        push @indels, "$curr_indel_start\t".($i-1);
    }
    $samples{$name."_depths"} = \@depths;
    $samples{$name."_indels"} = \@indels;
}


while (($key, $value) = each %samples) {
    my @array = @$value;
    if ($key =~ /depths/) {
        open FH, ">", "$outfile"."_$key.txt";
        print FH "pos\t$key\n";
        my $sum = 0;
        my $i;
        for ($i=1; $i<=@array; $i++) {
            if ($samplesize) {
                if (($i%$samplesize)==0) {
                    print FH "$i\t".@array[$i-1]."\n";
                }
            } else {
                print FH "$i\t".@array[$i-1]."\n";
            }
        }
        close FH;
    } else {
        open FH, ">", "$outfile"."_$key.txt";
        my $sum = 0;
        my $i;
        for ($i=1; $i<=@array; $i++) {
            print FH "$i\t".@array[$i-1]."\n";
        }
        print FH "size\t1\t$circle_size\n";
        close FH;
    }
}

print "drawing graphs...\n";
my $circlegraph_obj = new CircleGraph();
my $max=0;
my %graphs;
foreach my $name (@names) {
    #open the depths file
    open FH, "<", "$outfile"."_$name"."_depths.txt";
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
    $graphs{$name} = \@yvals;
}

#process yvals so that they are scaled
foreach my $name (@names) {
    my @yvals = @{$graphs{$name}};
    for (my $y=0;$y<@yvals;$y++) {
        @{$graphs{$name}}[$y] = @yvals[$y]/$max;
    }
}

my $j = 0;

# draw the coverage maps
my @xvals = @{$graphs{"x"}};
foreach my $name (@names) {
    my @yvals = @{$graphs{$name}};
    $circlegraph_obj->set_color($j);
    $circlegraph_obj->plot_line(\@xvals, \@yvals);
    $j++;
}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius);
$circlegraph_obj->draw_circle($circlegraph_obj->outer_radius);
$circlegraph_obj->set_font("Helvetica", 6, "black");
for (my $i = 0; $i < @xvals; $i++) {
    my $angle = (@xvals[$i]/$circle_size) * 360;
    my $radius = $circlegraph_obj->outer_radius + 10;
    my $label = @xvals[$i];
    $circlegraph_obj->circle_label($angle, $radius, "$label");
}

# draw the gene map
if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
        my ($fh, $filename) = tempfile(UNLINK => 1);
		print $fh get_locations_from_genbank_file($gb_file);
		close $fh;
	    $gb_file = "$filename";
	}
	draw_gene_map ($gb_file, $circlegraph_obj);
}

print "drawing reference mismatch maps...\n";
# draw the indel maps
$j = 0;
$circlegraph_obj->inner_radius($circlegraph_obj->inner_radius - 150);
foreach my $name (@names) {
    $circlegraph_obj = draw_regions ( "$outfile"."_$name"."_indels.txt", $circlegraph_obj, $j, $circlegraph_obj->inner_radius + ($j*5) );
    $j++;
}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius-5);

# draw the labels
$j = 0;
foreach my $name (@names) {
    my $label = $name;
    if (exists $labels->{$name}) {
        $label = $labels->{$name};
    }
    $circlegraph_obj->append_to_legend($label,$j++);
}

$circlegraph_obj->append_to_legend("Maximum coverage was $max, scaled to 1");
$circlegraph_obj->draw_legend_text({size => 10, height => 15});

open FH, ">", "$outfile.ps" or die "couldn't open $outfile.ps";
print FH $circlegraph_obj->output_ps();
close FH;

if (!$keepfiles) {
    while (($key, $value) = each %samples) {
        unlink "$outfile"."_$key.txt" or warn "could not delete $outfile"."_$key.txt";
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

=head1 DESCRIPTION

Graphs coverage depths using depth data generated from samtools mpileup.

=cut

