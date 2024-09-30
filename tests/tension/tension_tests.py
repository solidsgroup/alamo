import os
import shutil
import re
from PIL import Image

# Paths to the folder containing images and the input file template
image_folder = '/home/thoopul/alamo/tests/tension/images/'
input_file = '/home/thoopul/alamo/tests/tension/input'  # Original input file
tmp_file = '/home/thoopul/alamo/tests/tension/tmp_input'  # Temporary file for modification

# Step 1: Find all images in the image folder
image_files = [f for f in os.listdir(image_folder) if f.endswith('.png')]

# Step 2: Process each image one by one
for image_file in image_files:
    # Load the image and get its dimensions
    image_path = os.path.join(image_folder, image_file)
    with Image.open(image_path) as img:
        width, height = img.size
        print(f"Processing {image_file} with dimensions: {width}x{height}")
    
    # Step 3: Copy the original input file to a temp file for modification
    shutil.copy(input_file, tmp_file)
    
    # Step 4: Read and modify the temp file content
    with open(tmp_file, 'r') as file:
        modified_input = file.read()

    # Step 5: Modify the geometry lines using regex for precise replacement
    prob_lo_new = f"geometry.prob_lo    = -{width/2:.1f} -{height/2:.1f} 0.0"
    prob_hi_new = f"geometry.prob_hi    = {width/2:.1f} {height/2:.1f} 0.0"

    # Use regex to find and replace the geometry lines
    modified_input = re.sub(r'geometry\.prob_lo\s*=\s*[-\d\.]+\s*[-\d\.]+\s*0\.0', prob_lo_new, modified_input)
    modified_input = re.sub(r'geometry\.prob_hi\s*=\s*[-\d\.]+\s*[-\d\.]+\s*0\.0', prob_hi_new, modified_input)

    # Step 6: Modify boundary condition values by multiplying 0.1 with the width dimension
    modified_input = modified_input.replace('mechanics.bc.constant.val.yhi    = 0.0 0.1',
                                            f'mechanics.bc.constant.val.yhi    = 0.0 {0.1 * width}')
    modified_input = modified_input.replace('mechanics.bc.constant.val.xloyhi = 0.0 0.1',
                                            f'mechanics.bc.constant.val.xloyhi = 0.0 {0.1 * width}')
    modified_input = modified_input.replace('mechanics.bc.constant.val.xhiyhi = 0.0 0.1',
                                            f'mechanics.bc.constant.val.xhiyhi = 0.0 {0.1 * width}')

    # Step 7: Modify image path in euler.ic.eulerangles.filename line
    modified_input = re.sub(r'euler\.ic\.eulerangles\.filename\s*=\s*.+', 
                            f'euler.ic.eulerangles.filename = {image_path}', 
                            modified_input)

    # Step 8: Modify plot_file to `output_imagename`
    output_name = f"output_{os.path.splitext(image_file)[0]}"
    modified_input = modified_input.replace('plot_file=tests/tension/output',
                                            f'plot_file=./new/{output_name}')
                                            
    if width>height:
       modified_input = modified_input.replace('euler.ic.eulerangles.fit      = fitwidth',f'euler.ic.eulerangles.fit      = stretch')
    
    
    # Step 9: Write the modified input to the temp file
    with open(tmp_file, 'w') as file:
        file.write(modified_input)
    
    # Step 10: Run the code using the modified temp file
    os.system(f"/home/thoopul/alamo/bin/tension-2d-g++ {tmp_file}")
    
    print(f"Modified input file for {image_file}, updated geometry, boundary conditions, plot_file, and ran the code.")

