from visit import *
import os
import sys
import time
import shutil

# Ensure the script is executed with an output directory argument
if len(sys.argv) < 2:
    print("Usage: {} <output_directory>".format(sys.argv[0]))
    sys.exit(1)

# Define base directory within current working directory
base_dir = os.path.join(os.getcwd(), "extracted_curves")

# Ensure the directory exists
if not os.path.exists(base_dir):
    os.makedirs(base_dir)

# Open the .visit database from the current directory
OpenDatabase("localhost://{}/celloutput.visit".format(os.getcwd()))

# Define common lineout points
start_point = (0, 0.15)
end_point   = (1, 0.15)

# Dictionary to store variable names and corresponding filenames
variables = {
    "Pressure": "pressuredata",
    "PrimitiveVec001": "densitydata",
    "PrimitiveVec002": "velocitydata",
    "PrimitiveVec004": "internaldata"
}

for var, filename in variables.items():
    # --- Clear plots from both windows ---
    # Clear plots in main window (usually window 1)
    SetActiveWindow(1)
    DeleteAllPlots()
    # Clear plots in the lineout (curve) window (usually window 2)
    SetActiveWindow(2)
    DeleteAllPlots()
    
    # --- Create a pseudocolor plot for the variable in window 1 ---
    SetActiveWindow(1)
    AddPlot("Pseudocolor", var)
    DrawPlots()
    time.sleep(1)  # Allow VisIt time to process the operation
    
    # Set the time slider to the latest state
    n = TimeSliderGetNStates()
    SetTimeSliderState(n - 1)
    
    # --- Create the lineout which generates a curve in window 2 ---
    Lineout(start_point, end_point)
    DrawPlots()
    time.sleep(1)  # Ensure VisIt updates the visualization
    
    # Activate the lineout window
    SetActiveWindow(2)
    RedrawWindow()
    
    # Check that the lineout produced a plot
    if GetNumPlots() == 0:
        print(f"Error: No plots found for {var} after Lineout operation. Export aborted.")
        continue
    
    # Ensure the curve plot is active
    SetActivePlots(0)
    
    # --- Set up export attributes for Curve2D ---
    e = ExportDBAttributes()
    e.db_type = "Curve2D"
    e.db_type_fullname = "Curve2D"
    e.filename = filename
    e.variables = (var,)
    e.opts.types = (0,)
    
    try:
        ExportDatabase(e)
        print(f"Successfully exported {var} to {filename}")
    except Exception as ex:
        print(f"Export failed for {var}: {ex}")
        continue
        
    # --- Move the exported file to the base directory ---
    source_path = os.path.join(os.getcwd(), filename + ".curve")
    destination_path = os.path.join(base_dir, os.path.basename(source_path))
    
    try:
        shutil.move(source_path, destination_path)
        print(f"Moved {source_path} to {destination_path}")
    except Exception as move_ex:
        print(f"Failed to move {source_path}: {move_ex}")

