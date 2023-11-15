import os
import shutil
from decimal import Decimal

source_directory = '../airfoil_OFJ_static/static_URANS'
target_directory = './'
new_directory_name = '0'

# Function to extract the number from each directory name
def extract_number(directory_name):
    try:
        # Try to parse the directory name as a Decimal number
        return Decimal(directory_name)
    except:
        # Return a default value for non-numeric directory names
        return Decimal('-Infinity')

# List all directories in the source directory
directories = os.listdir(source_directory)

if directories:
    # Sort directories using the number extraction key
    largest_directory = max(directories, key=extract_number)

    # Build absolute source and target paths
    source_path = os.path.join(source_directory, largest_directory)
    target_path = os.path.join(target_directory, new_directory_name)

    # Copy the largest source directory to the target directory with the new name
    shutil.copytree(source_path, target_path)

    print(f"Directory '{largest_directory}' successfully copied to '{target_directory}' with the name '{new_directory_name}'.")
else:
    print("No directories found in the source directory.")
