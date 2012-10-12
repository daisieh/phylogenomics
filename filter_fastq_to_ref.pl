use Getopt::Long;
use Pod::Usage;
use IPC::Open3;
use File::Spec;
use IO::File;
use Symbol qw(gensym);

my ($fastqfile1, $fastqfile2, $reference, $resultfile) = 0;

GetOptions ('1=s' => \$fastqfile1,
            '2=s' => \$fastqfile2,
            'reference=s' => \$reference,
            'outputfile:s' => \$resultfile) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

open my $refFH, "<", "$reference.fai";
if ((-e $refFH) != 1) {
    print "samtools reference index didn't exist, making reference index.\n";
    system("samtools faidx $reference");
}

system("bwa aln -t 8 $reference $fastqfile1 > $fastqfile1.sai");
system("bwa aln -t 8 $reference $fastqfile2 > $fastqfile2.sai");

my $command = qq{bwa sampe $reference $fastqfile1.sai $fastqfile2.sai $fastqfile1 $fastqfile2 | \
awk '$2 != 77 && $2 != 141' | \
samtools view -S -u -h -t $reference.fai - | \
samtools sort - $resultfile.sorted};

system($command);
system("samtools index $resultfile.sorted.bam");


__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

filter_cp [options]

=head1 OPTIONS

    -fasta:     fasta file to filter
    -reference: fasta file of the reference to be filtered against
    -outputfile: output file name
    -evalue:    sets the evalue for blastn (default is 10)
    --blast/noblast:    runs blastn or uses previously generated filtered list (default is to run blastn)

=head1 DESCRIPTION

=cut
