#!/bin/bash

# List of different values for spall.po
spall_po_values=("1.0" "0.8" "0.6" "0.0" "-1.0")

# Loop over each value of spall.po
for po in "${spall_po_values[@]}"; do
  # Create a temporary input file
  temp_input="temp_input_file.txt"
  cp tests/void2/input $temp_input
  
  # Modify the temporary input file with the current value of spall.po
  sed -i "s/^spall.po = .*/spall.po = $po/" $temp_input
  
  # Run the simulation with the current value of spall.po
  ./bin/voidpf2-2d-g++ $temp_input
  
  # Rename the automatically created 'output' folder to 'output_po_$po'
  output_dir="tests/void2/output_po_$po"
  mv tests/void2/output $output_dir
  
  echo "Simulation complete for po=$po. Results saved in $output_dir"
  
  # Delete the temporary input file
  rm $temp_input
done

echo "All simulations for spall.po completed."

