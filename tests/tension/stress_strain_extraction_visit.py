import os
from PIL import Image
import sys
import shutil  # Import shutil to move files
import re  # For extracting the image ID

# Function to extract the numeric part from the directory name
def extract_image_id(directory_name):
    match = re.search(r'\d+_\d+_\d+', directory_name)
    if match:
        return match.group(0)
    else:
        return "unknown_id"

# Function to get min and max values for a given variable
def get_min_max_values(variable):
    Query("MinMax")
    min_max_vals = GetQueryOutputValue()
    min_val = min_max_vals[0]  # Min value is the first element
    max_val = min_max_vals[1]  # Max value is the second element
    return min_val, max_val

# Function to save the min and max values to a single file
def save_min_max_values(variable, min_value, max_value, filename="all_micrographs_min_max.txt"):
    with open(filename, "a") as file:  # Using "a" to append to the file
        file.write(f"{image_id}_Minimum {variable}: {min_value}\n")
        file.write(f"{image_id}_Maximum {variable}: {max_value}\n")

# Function to process the image: crop to stress distribution, remove background, and resize
def process_image(image_path, final_size, bg_color=(255, 255, 255)):
    with Image.open(image_path) as img:
        img = img.convert("RGBA")
        datas = img.getdata()

        # Detecting edges/boundaries of the stress distribution
        left, top, right, bottom = img.width, img.height, 0, 0
        for y in range(img.height):
            for x in range(img.width):
                pixel = img.getpixel((x, y))
                if pixel[0] < bg_color[0] or pixel[1] < bg_color[1] or pixel[2] < bg_color[2]:  # Non-background pixel
                    if x < left:
                        left = x
                    if x > right:
                        right = x
                    if y < top:
                        top = y
                    if y > bottom:
                        bottom = y

        # Crop the image to the detected boundaries
        if left < right and top < bottom:
            img = img.crop((left, top, right, bottom))

        # Make remaining white background pixels transparent
        new_data = [(255, 255, 255, 0) if (pixel[0] > 200 and pixel[1] > 200 and pixel[2] > 200) else pixel for pixel in img.getdata()]
        img.putdata(new_data)

        # Resize the image to the final desired size
        img = img.resize(final_size, Image.ANTIALIAS)

        # Save the processed image
        img.save(image_path, "PNG")

# List of 2D stress and strain components
components = ['stress_xx', 'stress_xy', 'stress_yy', 'strain_xx', 'strain_xy', 'strain_yy']

# Define the base directory where your output folders are located
base_directory = "/home/thoopul/alamo/tests/tension/new"

# Loop through all folders starting with 'output'
for folder_name in os.listdir(base_directory):
    if folder_name.startswith("output"):
        output_directory = os.path.join(base_directory, folder_name)

        # Open the database 
        visit_file = os.path.join(output_directory, "celloutput.visit")
        if not os.path.exists(visit_file):
            print(f"Visit file not found in {output_directory}, skipping this folder.")
            continue

        OpenDatabase(visit_file)

        # Extract the image ID from the directory name
        image_id = extract_image_id(output_directory)

        # Create results folder inside this output directory
        results_folder = os.path.join(base_directory, "results_final")
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)

        for component in components:
            # Add a pseudocolor plot for each stress and strain component
            AddPlot("Pseudocolor", component)

            # Draw the plot
            DrawPlots()

            # Step 1: Move to the final time step
            TimeSliderSetState(TimeSliderGetNStates() - 1)

            # Step 2: Remove annotations (legend, axes, title, user info, etc.)
            a = AnnotationAttributes()
            a.axes3D.visible = 0  # Turn off 3D axes
            a.axes2D.visible = 0  # Turn off 2D axes
            a.axesArray.visible = 0  # Turn off array axes
            a.legendInfoFlag = 0  # Turn off the legend
            a.userInfoFlag = 0  # Turn off user info
            a.databaseInfoFlag = 0  # Turn off database info
            SetAnnotationAttributes(a)

            # Step 3: Get min and max values for the component
            min_val, max_val = get_min_max_values(component)
            print(f"Original Min for {component}: {min_val}")
            print(f"Original Max for {component}: {max_val}")

            # Optionally save min and max values to a file
            save_min_max_values(component, min_val, max_val)

            # Step 4: Normalize the component values using min and max
            DefineScalarExpression("normalized_" + component, f"({component} - {min_val}) / ({max_val} - {min_val})")
            ChangeActivePlotsVar("normalized_" + component)

            # Step 5: Query min and max after normalization
            Query("MinMax")
            normalized_min_max_vals = GetQueryOutputValue()
            normalized_min_val = normalized_min_max_vals[0]
            normalized_max_val = normalized_min_max_vals[1]
            print(f"Normalized Min for {component}: {normalized_min_val}")
            print(f"Normalized Max for {component}: {normalized_max_val}")

            # Step 6: Query the data bounds to fit the domain exactly
            Query("SpatialExtents")  # Query the spatial extents of the data
            extents = GetQueryOutputValue()
            xmin, xmax, ymin, ymax = extents[:4]

            # Step 7: Set the window coordinates to fit the exact bounds of the data
            v = View2DAttributes()
            v.windowCoords = (xmin, xmax, ymin, ymax)  # Use the data bounds to set window coordinates
            v.viewportCoords = (0, 1, 0, 1)  # Fill the entire plot area with the domain
            v.fullFrameActivationMode = v.Off  # Disable full-frame mode to prevent extra space
            SetView2D(v)

            # Step 8: Get the dimensions of the existing PNG image in the output directory
            existing_image = os.path.join(output_directory, f"euler.ic.eulerangles.filename___home_thoopul_alamo_tests_tension_images_{image_id}_AnglesRepairQuatBlur.png")
            if os.path.exists(existing_image):
                with Image.open(existing_image) as img:
                    img_width, img_height = img.size
                    print(f"Using dimensions from existing PNG image: {img_width}x{img_height}")
            else:
                img_width, img_height = 800, 800  # Default size if no image is found
                print(f"Using bad dimensions")

            # Step 9: Set SaveWindowAttributes to match the size of the existing PNG
            s = SaveWindowAttributes()
            s.outputDirectory = results_folder  # Ensure images are saved in the results folder
            s.fileName = f"{image_id}_{component}"  # Image filename based on component
            s.family = 0  # Don't append numbers to filenames
            s.format = s.PNG  # Save as PNG
            s.width = img_width  # Set width from existing PNG
            s.height = img_height  # Set height from existing PNG
            s.quality = 100  # Max quality for PNG
            s.resConstraint = s.NoConstraint  # No resolution constraint to match domain
            SetSaveWindowAttributes(s)

            # Step 10: Save the plot as an image at the final time step
            SaveWindow()

            # Step 11: Process the image (crop and resize)
            current_image_path = f"{image_id}_{component}.png"
            process_image(current_image_path, final_size=(img_width, img_height))

            # Step 12: Move the image to the results folder (if needed)
            target_image_path = os.path.join(results_folder, current_image_path)
            if os.path.exists(current_image_path):
                shutil.move(current_image_path, target_image_path)
                print(f"Moved {current_image_path} to {target_image_path}")

            # Optionally close the current plot to prepare for the next component
            DeleteAllPlots()

        # Close the database 
        CloseDatabase(visit_file)

# Exit VisIt after completing the script
sys.exit()

