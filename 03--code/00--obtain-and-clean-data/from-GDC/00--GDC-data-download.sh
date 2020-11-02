#!/bin/zsh
# Author: Leo Williams
# Date: 2019-11-06
# Platform: macOS

# About ########################################

	# This script takes as input a TCGA data download manifest and a NIH/GDC secure access token, and downloads the manifest files

# Usage ########################################

	# zsh <script.sh> <manifest_file.txt> <gdc_token.txt> <output_dir>

# Setup ########################################

	# arguments
	MANIFEST_FILE=$1
	DOWNLOAD_DIR=$2
	AUTHORIZATION_TOKEN=$3

# Code ########################################

	# First, obtain a download manifest from the Genomic Data Commons (GDC) using the web portal

	# Download files in a manifest from TCGA using the gdc-client program
	gdc-client download \
	-m "$MANIFEST_FILE" \
	-d "$DOWNLOAD_DIR" \
	-t "$AUTHORIZATION_TOKEN" \
	--retry-amount 99999 \
	--wait-time 5 \
	> "$DOWNLOAD_DIR"/out-gdc-dl.log 2>&1 &

	# --log-file "$DOWNLOAD_DIR"/out-gdc-dl.log \
