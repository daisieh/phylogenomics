samtools view $1 | awk '{print ">"$1"\n"$10}' - > $2.fasta
