samtools view $1 | awk '{if ($2 == 77) print ">"$1"/1\n"$10; if ($2 == 141) print ">"$1"/2\n"$10; if (($2 != 77) && ($2 != 141)) print ">"$1"\n"$10;}' - > $2.fasta
