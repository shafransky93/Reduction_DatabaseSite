#!/usr/bin/env python
# coding: utf-8
print('#--------------------------------------------------------------#')
print('#                  Bens Astronomical Data                      #')
print('#         Addition Subtraction Stacking Software               #')
print('#                                                              #')
print('#                                                              #')
print('######################PATH STRUCTURE############################')
print('################################################################')
print('##                             data                           ##')
print('##                     _________|____________                 ##')
print('##                    /                      \                ##')
print('##                   x                      bias              ##')
print('##              _____|______                 |                ##')
print('##             /            \                |                ##')
print('##          light          flat           trimmed             ##')
print('##            |              |               |                ##')
print('##         reduced        reduced         reduced             ##')
print('##                                                            ##')
print('################################################################')
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
x_dir = os.path.dirname(__file__)+'/light'
x_flat_dir = os.path.dirname(__file__)+'/flat'
bias_dir =os.path.dirname(__file__)+'/bias'
x_reduced_dir = os.path.dirname(__file__)+'/light/reduced'
bias_trimmed_dir = os.path.dirname(__file__)+'/bias/trimmed'
x_flat_reduced_dir = os.path.dirname(__file__)+'/flat/reduced'
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



print('Image file collections sucessfully built...')


#x_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#x_flat_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']


for ccd, file_name in bias_collection.ccds(imagetyp='BIAS',    # Just get the bias frames
                                return_fname=True              # Provide the file name too.
                               ):    
   # Trim the overscan
   ccd = ccdproc.trim_image(ccd[:, :3078])    
   print('Ovserscan sucessfully trimmed from bias,')
   # Save the result
   ccd.write(os.path.join(bias_trimmed_dir, file_name+"trimmed.fits"), overwrite=True)
   print('Trimmed bias saved as ',os.path.join(bias_trimmed_dir, file_name+"trimmed.fits"))


#combine bias images
bias_trimmed_dir = os.path.dirname(__file__)+'/bias/trimmed'
bias_trimmed_collection = ccdproc.ImageFileCollection(bias_trimmed_dir)
print('Trimmed bias directory and image file collection sucessfully built...')

calibrated_biases = bias_trimmed_collection.files_filtered(imagetyp='BIAS', include_path=True)

combined_bias = ccdproc.combine(calibrated_biases,
                                method='average',
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5,
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median,
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

combined_bias.meta['combined'] = True
print('Bias sucessfully stacked...')
#Save combined bias in bias directory
combined_bias.write(os.path.join(bias_dir, 'combined_bias.fits'), overwrite=True)
print('Master bias sucessfully saved as: combined_bias.fits')

# Reduce red flats
for ccd, file_name in x_flat_collection.ccds(imagetyp='DOME FLAT',       # Just get the bias frames
                                             return_fname=True           # Provide the file name too.
                                             ):
    # Trim the overscan
    ccd = ccdproc.trim_image(ccd[:, :3078])
    # Subtract the overscan ##Causes weird error where it will delete data if left in 
    ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, :3095], median=True)
    # Subtract the bias 
    ccd = ccdproc.subtract_bias(ccd, combined_bias)
    print('Single channel overscan trimmed and master bias subtracted sucessfully,')
    # Save the result    
    ccd.write(os.path.join(x_flat_dir, 'reduced', 'x_flat_reduced_'+file_name), overwrite=True)    
    print('Reduced Single channel flat saved as: ',os.path.join(x_flat_dir, 'reduced', 'x_flat_reduced_'+file_name))

print('Single channel flats sucessfully reduced...')



x_flat_reduced_dir = os.path.dirname(__file__)+'/flat/reduced'
x_flat_reduced_collection = ccdproc.ImageFileCollection(x_flat_reduced_dir)
print('Reduced flat directories and image file collections sucessfully built...')


def inv_median(a):
    return 1/np.median(a)

#Stack Single channel flats
x_to_combine = x_flat_reduced_collection.files_filtered(imagetyp='DOME FLAT', include_path=True)
x_combined_flat = ccdproc.combine(x_to_combine,
                                method='average', 
                                scale= inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

x_combined_flat.meta['combined'] = True
print('Single channel flats stacked sucessfully,')
x_combined_flat.write(os.path.join(x_flat_reduced_dir, 'stacked_x_flat_reduced.fits'), overwrite=True)
print('Single channel master flat saved as: stacked_x_flat_reduced.fits')


# Reduce Single channel images and save
for light, file_name in x_collection.ccds(imagetyp='OBJECT', return_fname=True):
    x_light_ccds = []
    x_light_ccds.append(light)

    x_reduced = ccdproc.trim_image(light[:, :3078])
    x_reduced = ccdproc.subtract_bias(x_reduced, combined_bias)
    x_reduced = ccdproc.flat_correct(x_reduced, x_combined_flat)

    all_x = []
    all_x.append(x_reduced)

    print('Single channel light overscan trimmed, master bias and master Single channel flat subtraced' )
    x_reduced.write(os.path.join(x_reduced_dir, "x_light_reduced_"+file_name), overwrite=True)
    print('Reduced Single channel light saved as: ',os.path.join(x_reduced_dir, "x_light_reduced_"+file_name))

print('Single channel sucessfully reduced...')


# Stack Single channel images and save
x_reduced_dir = os.path.dirname(__file__)+'/light/reduced'
x_reduced_collection = ccdproc.ImageFileCollection(x_reduced_dir)

home_dir = os.path.dirname(__file__)
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
x_combined_light.write(os.path.join(home_dir,'stacked_x_light_reduced.fits'), overwrite=True)

print('Single channel light stacked sucessfully,')
print('Single channel master light saved as: stacked_x_light_reduced.fits')



# Confermation of science being done
for x in range(10):
   print('Science!')
