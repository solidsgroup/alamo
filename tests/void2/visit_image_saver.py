# Import VisIt functions
from visit import *
import os
import sys
import glob

# Function to find all directories with "node" in their names
def find_node_directories(base_path):
    node_dirs = []
    for root, dirs, files in os.walk(base_path):
        for dir_name in dirs:
            if "node" in dir_name:
                node_dirs.append(os.path.join(root, dir_name))
    return node_dirs

# Correct base directory
base_dir = "/home/thoopul/alamo/tests/void2"

# Loop through each directory starting with "output"
for output_folder in glob.glob(os.path.join(base_dir, "output")):
    
    # Find all directories containing "node" within the current output folder
    node_directories = find_node_directories(output_folder)

    # Create a directory to save images in the current output folder
    image_output_dir = os.path.join(output_folder, "visit_images")
    if not os.path.exists(image_output_dir):
        os.makedirs(image_output_dir)
        print(f"Created image output directory: {image_output_dir}")
    else:
        print(f"Image output directory already exists: {image_output_dir}")

    # Loop through each "node" directory
    for node_dir in node_directories:
        header_file = os.path.join(node_dir, "Header")
        
        if os.path.exists(header_file):
            print(f"Processing: {header_file}")
            
            OpenDatabase("localhost:" + header_file, 0)

            # Add the desired plots
            AddPlot("Pseudocolor", "d", 1, 1)
            AddPlot("Contour", "d", 1, 1)
            DrawPlots()

            # Save window settings
            SaveWindowAtts = SaveWindowAttributes()
            SaveWindowAtts.outputToCurrentDirectory = 0
            SaveWindowAtts.outputDirectory = image_output_dir
            SaveWindowAtts.fileName = os.path.basename(node_dir)  # Save with the directory name
            SaveWindowAtts.family = 0
            SaveWindowAtts.format = SaveWindowAtts.PNG
            SaveWindowAtts.width = 1024
            SaveWindowAtts.height = 1024
            SaveWindowAtts.screenCapture = 0
            SaveWindowAtts.saveTiled = 0
            SaveWindowAtts.quality = 80
            SaveWindowAtts.progressive = 0
            SaveWindowAtts.binary = 0
            SaveWindowAtts.stereo = 0
            SaveWindowAtts.compression = SaveWindowAtts.NONE
            SaveWindowAtts.forceMerge = 0
            SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions
            SaveWindowAtts.pixelData = 1
            SaveWindowAtts.advancedMultiWindowSave = 0
            
            # Save the window
            SetSaveWindowAttributes(SaveWindowAtts)
            SaveWindow()

            print(f"Image saved for {header_file}")

            DeleteAllPlots()

# Close the VisIt session
CloseComputeEngine()

sys.exit("All images saved successfully.")

