#!/usr/bin/env bash

IN_DIR=$1 # directory with (.vcf).gz files
OUT_DIR=$2

for i in "$IN_DIR"/*.gz; do
	echo "$i"
	
	BASENAME=$(basename "$i")

	gunzip -c "$i" > "$OUT_DIR"/"${BASENAME%.gz}" & 
done
