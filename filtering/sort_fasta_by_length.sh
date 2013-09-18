# sorts a fasta file by size of seqs, largest first.

infile=$1
gawk '{if (NF == 0) next; s = ""; for (i=2;i<=NF;i++) s = s$i; print length(s)","$1","s}' RS=">" $infile | sort -n -r | gawk '{print ">"$2"\n"$3}' FS=","
