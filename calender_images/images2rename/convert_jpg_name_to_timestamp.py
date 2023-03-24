import os
from PIL import Image

input("Hit any Key to rename all image files in this directory to the date they were taken")

# Sets path to the directory program is opened from
path = os.path.dirname(__file__)

# Pull date from .jpg header
# JPG TIMESTAMP ID:36867
def get_date_taken(path):
    return Image.open(path)._getexif()[36867]

date_counts = {}
# Return all files as a list
for file in os.listdir(path):
    # Check the files which end with specified extension
    if file.endswith(".jpg"):
        # Renames the files in the list that end with .jpg
        # Removes special characters and spaces and renames as a .jpg file
        date_taken = get_date_taken(file).split(' ')[0]  # Extract date in YYYY:MM:DD format
        date_taken = date_taken.replace(':', '-')  # Replace ':' with '-' to get YYYY-MM-DD format
        
        if date_taken not in date_counts:
            date_counts[date_taken] = 0
            
        date_counts[date_taken] += 1
        
        if date_counts[date_taken] == 1:
            new_filename = f"{date_taken}.jpg"
        else:
            new_filename = f"{date_taken}_{date_counts[date_taken]}.jpg"
        
        os.rename(file, new_filename)
