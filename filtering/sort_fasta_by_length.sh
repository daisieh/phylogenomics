infile=$1
outfile=$1.sorted.fasta
gawk '{if ($0 - /^$/) next; s = ""; for (i=2;i<=NF;i++) s = s$i; print $1","s}' RS=">" $infile | gawk '{s = ""; for (i=1;i<=NF;i++) s = s"_"$i;print $4 "," s}' FS="_" | sort -n -r | gawk '{sub ("_","",$2); print ">" $2 "\n" $3 "\n"}' FS="," > $outfile
