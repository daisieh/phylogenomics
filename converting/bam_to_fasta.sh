samtools view $1 | gawk '{if (and($2,0x0200)) next; if (and($2,0x0040)) print ">"$1"/1\n"$10; if (and($2,0x0080)) print ">"$1"/2\n"$10; else print ">"$1"\n"$10;}' - > $2.fasta
