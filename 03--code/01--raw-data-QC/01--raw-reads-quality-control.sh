#!/usr/bin/env zsh

# FILENAME:		[FILENAME]
# AUTHOR:		[AUTHOR]
# DATE:			[DATE]
# USAGE:		[USAGE]
# DESCRIPTION:	See https://stackoverflow.com/questions/430078/shell-script-templates for template ideas
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	[DEV PLATFORM]

# ==============================================================================

# ARGUMENTS ########################################

IN_DIR=$1 # has FASTQC
# OUTPUT_DIR=$2

# SETUP ########################################

# # Make timestamped dir for results output
# RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
# mkdir $RESULTS_DIR

# # Make a copy of the script run and put it into the results dir
# SCRIPT_PATH=$(echo "$0")
# cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

