#!/bin/bash

# Check if file path arguments are provided
if [ $# -ne 1 ]; then
	echo "Usage: $0 <file_path1>"
	exit 1
fi

# Get the file paths from command line arguments
FILE_PATH1="$1"

# Check if first file exists
if [ ! -f "$FILE_PATH1" ]; then
	echo "Error: File '$FILE_PATH1' does not exist"
	exit 1
fi

BASE_PATH="/usr1/aidanz/databases"
OUTPUT_PATH="/usr1/aidanz/projects/read_sim/data/profiles"

# Get basename of the input file
BASE_NAME=$(basename "$FILE_PATH1" .fastq)
BASE_NAME=$(basename "$BASE_NAME" .fq)
BASE_NAME=$(basename "$BASE_NAME" .fastq.gz)
BASE_NAME=$(basename "$BASE_NAME" .fq.gz)

kraken2 --db $BASE_PATH/kraken/bfdb --threads 32 --report $OUTPUT_PATH/kraken/${BASE_NAME}_bfdb.kreport --output $OUTPUT_PATH/kraken/${BASE_NAME}_bfdb.kraken $FILE_PATH1
