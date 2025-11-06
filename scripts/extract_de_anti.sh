#!/bin/bash

# Usage: bash extract.sh transcript_list.txt input.gtf input.bed

if [ "$#" -ne 3 ]; then
    echo "Usage: bash extract.sh transcript_list.txt input.gtf input.bed"
    exit 1
fi

TRANSCRIPTS="$1"
GTF="$2"
BED="$3"

GTF_OUT="$(basename "$GTF" .gtf)_significant.gtf"
BED_OUT="$(basename "$BED" .bed)_significant.bed"

# Extract matching lines from GTF
grep -wFf "$TRANSCRIPTS" "$GTF" > "$GTF_OUT"

# Extract matching lines from BED
grep -wFf "$TRANSCRIPTS" "$BED" > "$BED_OUT"