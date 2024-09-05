import os
import imageio

# Base directory where the "output" folders are located
base_dir = "/home/thoopul/alamo/tests/void2"

# Loop through each folder starting with "output"
output_dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if d.startswith("output")]

for output_dir in output_dirs:
    print(f"Processing directory: {output_dir}")

    # Directory where your images are saved
    image_output_dir = os.path.join(output_dir, "visit_images")

    # Check if the directory exists
    if not os.path.exists(image_output_dir):
        print(f"Image directory does not exist: {image_output_dir}")
        continue

    # Get a sorted list of all PNG files in the directory
    image_files = sorted([os.path.join(image_output_dir, file) for file in os.listdir(image_output_dir) if file.endswith(".png")])

    if not image_files:
        print(f"No PNG files found in {image_output_dir}")
        continue

    # Name of the output GIF file
    gif_filename = os.path.join(output_dir, "output_animation.gif")

    # Create the GIF
    with imageio.get_writer(gif_filename, mode='I', duration=0.2) as writer:
        for filename in image_files:
            image = imageio.imread(filename)
            writer.append_data(image)

    print(f"GIF created successfully: {gif_filename}")

