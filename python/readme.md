To subset raw bam reads and use bwa to align them:
`bwa_to_bam.py -r reference.fasta -i input.txt -p 8`
where input.txt has a tab-delim list of `host	sample	path`
