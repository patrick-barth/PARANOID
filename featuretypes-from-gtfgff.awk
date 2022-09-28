#!/usr/bin/gawk -f

BEGIN { FS="\t"; feattypes[""]=0; }
NF > 2 && $1 !~ /^#/ { feattypes[$3]=1; }
$1 ~ /^##FASTA/ { exit }
END { delete feattypes[""]; for (i in feattypes) { print i; } }
