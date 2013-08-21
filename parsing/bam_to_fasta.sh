samtools view $1.bam | awk '{print ">$1\n$10"}' - > $1.fasta
