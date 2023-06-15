#!/usr/bin/gawk -f
#' Extract exons from UTRome GTF as BED
BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
    if ($3 ~ /exon/) {
	match($9, /transcript_id "([^"]+)"/, txid);
	print $1, $4, $5, txid[1], 0, $7;
    }
}
