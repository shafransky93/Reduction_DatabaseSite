#!/usr/bin/env python
# coding: utf-8
print('#--------------------------------------------------------------#')
print('#                  Bens Astronomical Data                      #')
print('#         Addition Subtraction Stacking Software               #')
print('#                                                              #')
print('#                                                              #')
print('#                                                              #')
print('# !BACKUP RAW DATA INTO SEPERATE FOLDER PRIOIR TO MANIPULATION!#')
print('#Data will not be overwriten, but good practice incase of error#')
print('################################################################')



import os
import numpy as np
from astropy.stats import mad_std
from astropy.nddata import CCDData
from astropy.io import fits
from pathlib import Path
import ccdproc
import warnings

print('Python directories sucessfully imported...')

warnings.filterwarnings("ignore")
print('Warnings suppressed')

##Direct input field for directory selection 
x_dir = os.path.dirname(__file__)
x_reduced_dir = os.path.dirname(__file__)+'/reduced'
print('File directories sucessuflly built...')

#Build file structure propigation
#Bild colections for:
#                          lights, flats, bias', overscan trimmed bias',
#                          reduced flats and reduced light images.
x_collection = ccdproc.ImageFileCollection(x_dir)
x_flat_collection = ccdproc.ImageFileCollection(x_flat_dir)
bias_collection = ccdproc.ImageFileCollection(bias_dir)
bias_trimmed_collection = ccdproc.ImageFileCollection(bias_trimmed_dir)
x_reduced_collection = ccdproc.ImageFileCollection(x_reduced_dir)
x_flat_reduced_collection = ccdproc.ImageFileCollection(x_flat_reduced_dir)



# Stack images and save
x_reduced_dir = os.path.dirname(__file__)
x_reduced_collection = ccdproc.ImageFileCollection(x_reduced_dir)
print('Reduced Single channel light directory and image file collection sucessfully built...')

x_light_to_combine = x_reduced_collection.files_filtered(imagetyp="OBJECT", include_path=True)
x_combined_light = ccdproc.combine(x_light_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

x_combined_light.meta['combined'] = True
x_combined_light.write(os.path.join(x_reduced_dir, "stacked_x.fits"), overwrite=True)

print('Single channel light stacked sucessfully,')
print('Single channel master light saved as: stacked_x.fits')



# Confermation of science being done
for x in range(10):
   print('Science!')
