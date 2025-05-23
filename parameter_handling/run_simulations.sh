#!/bin/bash

output_dir="OUTPUT_010225"
mkdir -p "$output_dir"

# Read each line of parameters.txt
while read param1 param2 param3; do
  # Create an unique name for the output file, and keep it in the created directory
  output_file="${output_dir}/OutputFile_${param1}_${param2}_${param3}.txt"

  # Execute the executable binary script and redirect the output to a file.
  ./simulate_spacers1.2 "$param1" "$param2" "$param3" > "$output_file"

  # Informative message the simulation was completed.
  echo "Simulation completed with parameters (alpha, gamma and pEndemic): $param1 $param2 $param3. Saved in $output_file"
  
done < parameters.txt

# Remember to give the shell script permission to execute: chmod +x run_simulations.sh

# In case you want to run the script in the background, you can use:
# nohup ./run_simulations.sh &

# To run in parallel:
# nohup ./run_simulations1.sh & nohup ./run_simulations2.sh & nohup ./run_simulations3.sh & nohup ./run_simulations4.sh & nohup ./run_simulations5.sh &
# In my case, I used 5 different scripts that run in parallel, each with a different set of 100 input parameters and output files.