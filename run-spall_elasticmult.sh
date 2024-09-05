#!/bin/bash

# List of different values for spall.elastic_mult
spall_elastic_mult_values=("-2.0" "-1.5" "-1.0" "0.0" "1.0" "2.0")

# Loop over each value of spall.elastic_mult
for elastic_mult in "${spall_elastic_mult_values[@]}"; do
  # Create a temporary input file
  temp_input="temp_input_file.txt"
  cp tests/void2/input $temp_input
  
  # Modify the temporary input file with the current value of spall.elastic_mult
  sed -i "s/^spall.elastic_mult = .*/spall.elastic_mult = $elastic_mult/" $temp_input
  
  # Run the simulation with the current value of spall.elastic_mult
  ./bin/voidpf2-2d-g++ $temp_input
  
  # Rename the automatically created 'output' folder to 'output_elastic_mult_$elastic_mult'
  output_dir="tests/void2/output_elastic_mult_$elastic_mult"
  mv tests/void2/output $output_dir
  
  echo "Simulation complete for elastic_mult=$elastic_mult. Results saved in $output_dir"
  
  # Delete the temporary input file
  rm $temp_input
done

echo "All simulations for spall.elastic_mult completed."

