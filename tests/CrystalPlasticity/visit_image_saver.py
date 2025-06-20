# Import visit functions here
from visit import *
import os
import sys
import glob

# Find all files with the word node
def find_node_directories(base_path):
    node_dirs = []
    for root, dirs, files in os.walk(base_path):
        for dir_name in dirs:
            if "node" in dir_name:
                node_dirs.append(os.path.join(root, dir_name))
    return node_dirs

# the base directory path
base_dir = "/home/thoopul/alamo/tests/CrystalPlasticity"

# loop through the output directory that you need to use
for output_folder in glob.glob(os.path.join(base_dir, "output_sr")):
    
    # Find all directories containing "node" within the output directory
    node_directories = find_node_directories(output_folder)

    # Create a directory to save images in the current output folder
    image_output_dir = os.path.join(output_folder, "visit_images_gamma")
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
            AddPlot("Pseudocolor", "gamma1", 1, 1)
           
            pc = PseudocolorAttributes()
            pc.minFlag = 1           # use a user-defined minimum
            pc.maxFlag = 1           # use a user-defined maximum
            pc.min      = 0.00     # <-- set your desired lower bound 
            pc.max      =  0.011      # <-- set your desired upper 
            pc.scaling  = pc.Linear  # or pc.Log for logarithmic scaling
            SetPlotOptions(pc)       # apply the settings

            #AddPlot("Contour", "d", 1, 1)
            DrawPlots()

            # Save window settings
            SaveWindowAtts = SaveWindowAttributes()
            SaveWindowAtts.outputToCurrentDirectory = 0
            SaveWindowAtts.outputDirectory = image_output_dir
            SaveWindowAtts.fileName = os.path.basename(node_dir) 
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

# Close visit
CloseComputeEngine()

sys.exit("All images saved successfully.")

