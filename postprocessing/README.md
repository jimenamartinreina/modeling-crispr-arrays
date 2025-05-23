# Postprocessing

This folder contains a Bash script designed to calculate the mean between two files with identical formats. These files are outputs from simulation runs conducted as part of this project.
The script is made to extract the files from two different directories.

## Purpose

Due to hardware limitations, it was not feasible to perform 10,000 simulation replicates in a single run. Instead, two separate rounds of 5,000 replicates were executed. The provided script merges the results by computing the mean between the corresponding outputs of these two runs.
This approach enables analysis of results equivalent to a single 10,000-replicate simulation, overcoming computational constraints.

## Usage

1. Ensure both simulation output directories are present in this directory.
2. Run the Bash script: 
```bash
./merge_files.sh
```
3. See the output files in the output directory you have defined in the script.