#!/bin/bash

# List of different values for spall.kappa
spall_kappa_values=("-0.0005" "-0.001" "-0.002" "0.001")

# Loop over each value of spall.kappa
for kappa in "${spall_kappa_values[@]}"; do
  # Create a temporary input file
  temp_input="temp_input_file.txt"
  cp tests/void2/input $temp_input
  
  # Modify the temporary input file with the current value of spall.kappa
  sed -i "s/^spall.kappa = .*/spall.kappa = $kappa/" $temp_input
  
  # Run the simulation with the current value of spall.kappa
  ./bin/voidpf2-2d-g++ $temp_input
  
  # Rename the automatically created 'output' folder to 'output_kappa_$kappa'
  output_dir="tests/void2/output_kappa_$kappa"
  mv tests/void2/output $output_dir
  
  echo "Simulation complete for kappa=$kappa. Results saved in $output_dir"
  
  # Delete the temporary input file
  rm $temp_input
done

echo "All simulations for spall.kappa completed."

