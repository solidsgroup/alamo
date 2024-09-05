#!/bin/bash

# List of different values for spall.cp
spall_cp_values=("1.0" "0.5" "0.2" "0.0")

# Loop over each value of spall.cp
for cp in "${spall_cp_values[@]}"; do
  # Create a temporary input file
  temp_input="temp_input_file.txt"
  cp tests/void2/input $temp_input
  
  # Modify the temporary input file with the current value of spall.cp
  sed -i "s/^spall.cp = .*/spall.cp = $cp/" $temp_input
  
  # Run the simulation with the current value of spall.cp
  ./bin/voidpf2-2d-g++ $temp_input
  
  # Rename the automatically created 'output' folder to 'output_cp_$cp'
  output_dir="tests/void2/output_cp_$cp"
  mv tests/void2/output $output_dir
  
  echo "Simulation complete for cp=$cp. Results saved in $output_dir"
  
  # Delete the temporary input file
  rm $temp_input
done

echo "All simulations for spall.cp completed."

