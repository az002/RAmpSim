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

kraken2 --db $BASE_PATH/kraken/zymo_genomes --threads 32 --report $OUTPUT_PATH/kraken/${BASE_NAME}_zymodb.kreport --output $OUTPUT_PATH/kraken/${BASE_NAME}_zymodb.kraken $FILE_PATH1
ganon classify -d $BASE_PATH/ganon/zymo -s $FILE_PATH1 -o $BASE_PATH/ganon/${BASE_NAME}_zymo
ganon report -d $BASE_PATH/ganon/zymo -i $BASE_PATH/ganon/${BASE_NAME}_zymo.rep -o $BASE_PATH/ganon/${BASE_NAME}_zymo_profile
ganon table -i $BASE_PATH/ganon/${BASE_NAME}_zymo_profile.tre -o $OUTPUT_PATH/ganon/${BASE_NAME}_zymo_profile.tsv --rank species --output-value percentage
sylph profile $BASE_PATH/sylph/zymo_genomes.syldb $FILE_PATH1 -t 32 -o $OUTPUT_PATH/sylph/${BASE_NAME}.tsv -c 50
centrifuger -x $BASE_PATH/centrifuger/zymo -u $FILE_PATH1 -t 32 > $OUTPUT_PATH/cfr/${BASE_NAME}_classification.tsv
centrifuger-quant -x $BASE_PATH/centrifuger/zymo -c $OUTPUT_PATH/cfr/${BASE_NAME}_classification.tsv > $OUTPUT_PATH/cfr/${BASE_NAME}_report.tsv