#!/usr/bin/env zsh

# USAGE ####################
# ./move-rename-unzip-fq.sh OUT_DIR
#
# run from within folder that contains sample folders that contain .fq files, e.g.:
#
# .
# ├── SAMPLE_DIR_1
# │   ├── fwd_1.fq.gz
# │   └── rev_2.fq.gz
# ├── SAMPLE_DIR_2
# │   ├── fwd_1.fq.gz
# │   └── rev_2.fq.gz
# └── move-rename-unzip-fq.sh

OUT_DIR=$1

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$OUT_DIR"

for i in * ; do
	# get sample name for each pair of fq files
	SAMPLE_NAME=$(basename "$i")

	# make new sample dir in output dir for the unzipped & renamed files
	OUT_SAMPLE_DIR="$OUT_DIR"/"$SAMPLE_NAME"
	mkdir "$OUT_SAMPLE_DIR"

	# rename and unzip FWD fq into new file
	gunzip -c "$i"/*_1.fq.gz > "$OUT_SAMPLE_DIR"/"${SAMPLE_NAME}_1.fq"
	# rename and unzip REV fq into new file
	gunzip -c "$i"/*_2.fq.gz > "$OUT_SAMPLE_DIR"/"${SAMPLE_NAME}_2.fq"
done
