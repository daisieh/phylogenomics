#!/usr/bin/perl
use strict;
use Pod::Usage;
use File::Basename;
use Getopt::Long;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($iconsfile, $inputfile, $out_file, $help) = 0;
GetOptions ('icons|colors=s' => \$iconsfile,
            'samples|input=s' => \$inputfile,
            'output=s' => \$out_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

open fileIN, "<:crlf", "$iconsfile" or die "no file named $iconsfile";
my @inputs = <fileIN>;

my %colors = ();

foreach my $color (@inputs) {
	(my $name, my $label) = split (/\s*\t/, $color);
	$label =~ s/\n|\r//;
	$colors{"$name"} = "$label";
}

open fileIN, "<:crlf", "$inputfile" or die "no file named $inputfile";
@inputs = <fileIN>;

my $result = "";

#print headers
$result .= "<\?xml version=\"1.0\" encoding=\"UTF-8\"\?>\n<kml xmlns=\"http://earth.google.com/kml/2.2\">\n<Document>
\t<name>localities.kml</name>\n\n";

foreach my $key ( keys %colors ) {
	my $label = $colors{$key};
	if ($label eq "") {
		$label = "brown";
	}
 	my $url = "http://www.mutantdaisies.com/images/icons/colors/$label.png";
# 	my $url = "http://http://maps.google.com/mapfiles/kml/pal3/icon12.png";
	$result .= "\t<StyleMap id=\"$label\">\n\t\t<Pair>\n\t\t\t<key>normal</key>\n\t\t\t<styleUrl>#$label</styleUrl>\n\t\t</Pair>\r";
	$result .= "\t\t<Pair>\n\t\t\t<key>highlight</key>\n\t\t\t<styleUrl>#$label</styleUrl>\n\t\t</Pair>\n\t</StyleMap>\n";
	$result .= qq{\t<Style id=\"$label\">\n\t\t<IconStyle>\n\t\t\t<Icon>\n\t\t\t\t<href>$url</href>\n\t\t\t\t<scale>1</scale>\n\t\t\t</Icon>\n\t\t\t<hotSpot x=\"0.5\" y=\"0.5\" xunits=\"fraction\" yunits=\"fraction\"/>\n\t\t</IconStyle>\n\t</Style>\n};
}

$result .= "\t<Folder>\n<name>localities</name>\n";
foreach my $line (@inputs) {
	print $line;
	(my $species, my $name, my $lat, my $long) = split (/\t/, $line);
	$long =~ s/\n|\r//;
	my $label = $colors{$species};
	if ($label eq "") {
		$label = "alb_alb";
	}
	if ($lat ne "") {
		$result .= "\t\t<Placemark>\n\t\t\t<name>$name</name>\n\t\t\t<styleUrl>#$label</styleUrl>\n";
		$result .= "\t\t\t<Point>\n\t\t\t\t<coordinates>$long,$lat</coordinates>\n\t\t\t</Point>\n\t\t</Placemark>\n";
	}
}

$result .= "\t</Folder>\n</Document>\n</kml>";

open (fileOUT, ">$out_file.kml") or die "couldn't make $out_file.kml\n";
truncate fileOUT, 0;
print fileOUT "$result\n";
close fileOUT;

__END__

=head1 NAME

make_kml

=head1 SYNOPSIS

make_kml [-colors colorfile] [-samples inputfile] [-output out_file]

=head1 OPTIONS

  -colors:          tab-delimited list of labels and colors
  -samples:         tab-delimited list of localities: label, name, lat, long
  -output:          prefix of output files

=head1 DESCRIPTION

Makes a kml file out of a list of sample localities, tab-delimited:
label, name, lat, long


=cut

