#!/bin/bash

# List of different values for spall.M
spall_M_values=("1.0" "0.5" "0.01")

# Loop over each value of spall.M
for M in "${spall_M_values[@]}"; do
  # Create a temporary input file
  temp_input="temp_input_file.txt"
  cp tests/void2/input $temp_input
  
  # Modify the temporary input file with the current value of spall.M
  sed -i "s/^spall.M = .*/spall.M = $M/" $temp_input
  
  # Run the simulation with the current value of spall.M
  ./bin/voidpf2-2d-g++ $temp_input
  
  # Rename the automatically created 'output' folder to 'output_M_$M'
  output_dir="tests/void2/output_M_$M"
  mv tests/void2/output $output_dir
  
  echo "Simulation complete for M=$M. Results saved in $output_dir"
  
  # Delete the temporary input file
  rm $temp_input
done

echo "All simulations for spall.M completed."

