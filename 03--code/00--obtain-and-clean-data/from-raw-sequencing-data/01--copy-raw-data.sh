#!/bin/bash

# copy all sample folders (containing raw data) in current folder to a destination project folder for working data

# navigate to folder with raw data sample folders

# define destination project data folder
DEST=$1

# cp all data folders
cp -nRv * $1

