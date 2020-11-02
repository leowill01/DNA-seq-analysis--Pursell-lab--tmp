#!/usr/bin/env bash

# ABOUT ####################

# !THIS SCRIPT MODIFIES FILES IN-PLACE, SO MAKE SURE THERE IS ANOTHER COPY SOMEWHERE!

# This script takes a dir that contains a dir for each sequencing sample, each of which has both fwd and rev .fq.gz read files split into multiple parts and combines and decompresses them into two single FWD and REV .fq read files for each sample.

# USAGE ####################
# ./move-rename-unzip-fq.sh IN_DIR
#
# run from anyhere, just specify the IN_DIR bc it modifies files relative to only that.
#
# IN_DIR
# ├── SAMPLE_DIR_1
# │   ├── fwd_1.fq.gz
# │   └── rev_2.fq.gz
# └── SAMPLE_DIR_2
#     ├── fwd_1.fq.gz
#     └── rev_2.fq.gz


IN_DIR=$1

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$IN_DIR"

for SAMPLE_DIR in "$IN_DIR"/*/ ; do
	# get sample name for each pair of fq files
	SAMPLE_NAME=$(basename "$SAMPLE_DIR")
	echo "Decompressing, combining, and renaming reads for:"
	echo "$SAMPLE_NAME"

	# unzip and concatenate both FWD read files ending in _1.fq.gz into a single FWD read .fq file and rename with sample name
	gunzip -c "$IN_DIR"/"$SAMPLE_NAME"/*_1.fq.gz > "$IN_DIR"/"$SAMPLE_NAME"/"${SAMPLE_NAME}-1fwd.fq"

	# same for the reverse reads
	gunzip -c "$IN_DIR"/"$SAMPLE_NAME"/*_2.fq.gz > "$IN_DIR"/"$SAMPLE_NAME"/"${SAMPLE_NAME}-2rev.fq"

	# remove .gz files
	rm "$IN_DIR"/"$SAMPLE_NAME"/*.gz
done
