import os
from PIL import Image

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


# Open the database 
OpenDatabase("/home/thoopul/alamo/tests/tension/output_6_009_00008_AnglesRepairQuatBlur/celloutput.visit")  

# Add a pseudocolor plot for 'stress_xx'
AddPlot("Pseudocolor", "stress_xx")

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

# Step 3: Query the data bounds to fit the domain exactly
Query("SpatialExtents")  # Query the spatial extents of the data
extents = GetQueryOutputValue()
xmin, xmax, ymin, ymax = extents[:4]

# Step 4: Set the window coordinates to fit the exact bounds of the data
v = View2DAttributes()
v.windowCoords = (xmin, xmax, ymin, ymax)  # Use the data bounds to set window coordinates
v.viewportCoords = (0, 1, 0, 1)  # Fill the entire plot area with the domain
v.fullFrameActivationMode = v.Off  # Disable full-frame mode to prevent extra space
SetView2D(v)

# Step 4: Get the dimensions of the existing PNG image in the output directory
output_directory = "/home/thoopul/alamo/tests/tension/output_6_009_00008_AnglesRepairQuatBlur"  
existing_image = os.path.join(output_directory, "euler.ic.eulerangles.filename___home_thoopul_alamo_tests_tension_images_6_009_00008_AnglesRepairQuatBlur.png") 
if os.path.exists(existing_image):
    with Image.open(existing_image) as img:
        img_width, img_height = img.size
        print(f"Using dimensions from existing PNG image: {img_width}x{img_height}")


# Step 5: Set SaveWindowAttributes to match the size of the existing PNG
s = SaveWindowAttributes()
s.outputDirectory = output_directory
s.fileName = "stress_xx_final_time"  # Image filename
s.family = 0  # Don't append numbers to filenames
s.format = s.PNG  # Save as PNG
s.width = img_width  # Set width from existing PNG
s.height = img_height  # Set height from existing PNG
s.quality = 100  # Max quality for PNG
s.resConstraint = s.NoConstraint  # No resolution constraint to match domain
SetSaveWindowAttributes(s)

# Step 6: Save the plot as an image at the final time step
SaveWindow()

# Step 7: Process the saved image to crop, remove background, and resize
output_image_path = os.path.join(output_directory, "stress_xx_final_time.png")  # Path of the saved image
process_image(output_image_path, final_size=(img_width, img_height))  # Crop, remove background, resize

# Optionally close the plot to clean up
DeleteAllPlots()

# Close the database 
CloseDatabase("/home/thoopul/alamo/tests/tension/output_6_009_00008_AnglesRepairQuatBlur/celloutput.visit")
EXit()

