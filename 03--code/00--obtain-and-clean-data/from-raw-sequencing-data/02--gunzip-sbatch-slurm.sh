#!/usr/bin/env bash

DIR_OF_GZS=$1
OUTPUT_DIR=$2

for i in "$DIR_OF_GZ"/*.gz; do
	echo "$i"
    sbatch 03--gunzip-single-sbatch.sh "$i" "$OUTPUT_DIR"
done