infile=$1
outfile=$1.sorted.fasta
gawk '{if ($0 - /^$/) next; s = ""; for (i=2;i<=NF;i++) s = s$i; print length(s)","$1","s}' RS=">" $infile | sort -n -r | gawk '{sub ("_","",$2); print ">" $2 "\n" $3 "\n"}' FS=","> $outfile
#| sort -n -r | gawk '{sub ("_","",$2); print ">" $2 "\n" $3 "\n"}' FS=","
