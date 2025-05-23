#!/bin/bash
# This script merges files from two directories by calculating the mean of corresponding columns.
# It assumes that the first 6 lines of the files are headers and should be printed as is.

# Directories with files to merge
DIR1="OUTPUT_210425"
DIR2="OUTPUT_220425"
# Output directory
DIR_OUT="OUTPUT_means_230425"
# The output directory will be created if it does not exist
mkdir -p "$DIR_OUT"

# Loop in one of the input folders (arbitrarily chosen as the first one)
for file in "$DIR1"/*; do
    filename=$(basename "$file")
    # Verify if the file exists in the second folder
    if [[ -f "$DIR2/$filename" ]]; then
        # Output file
        output_file="$DIR_OUT/$filename"

        # First lines of the output file and calculate means
        paste "$DIR1/$filename" "$DIR2/$filename" | awk '
        BEGIN { OFS="\t" } 
        NR <= 5 { print; next }  # Print first lines (the seed will be the one of the first file)
        NR == 6 { print $1, $2, $3, $4, $5, $6; next }
        { 
            print ($1 + $7) / 2, ($2 + $8) / 2, ($3 + $9) / 2, ($4 + $10) / 2, ($5 + $11) / 2, ($6 + $12) / 2
        }' > "$output_file"
        # Informative messagges
        echo "Processed: $filename"
    else
        echo "Waring: $filename not found in $DIR2"
    fi
done

# Informative message that the merging was completed.
echo "Final output files saved in $DIR_OUT"