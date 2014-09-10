#!/bin/bash

samtools view $1 | gawk '{print "@"$1"\n"$10"\n+\n"$11}' > $2.fastq
