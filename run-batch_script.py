import subprocess

# List of all batch scripts
batch_scripts = [
    "run-spall_M.sh",
    "run-spall_kappa.sh",
    "run-spall_po.sh",
    "run-spall_cp.sh",
    "run-spall_elasticmult.sh",
]

# Loop through and run each batch script
for script in batch_scripts:
    print(f"Running {script}...")
    try:
        # Run the batch script
        subprocess.run(["bash", script], check=True)
        print(f"{script} completed successfully.\n")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running {script}: {e}\n")

print("All batch scripts have been executed.")

