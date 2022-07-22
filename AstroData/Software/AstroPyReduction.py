#!/usr/bin/env python
# coding: utf-8
print('#--------------------------------------------------------------#')
print('#                  Bens Astronomical Data                      #')
print('#         Addition Subtraction Stacking Software               #')
print('#                                                              #')
print('#                                                              #')
print('######################PATH STRUCTURE############################')
print('################################################################')
print('####                           Data                          ###')
print('##         _____________________|______________________       ##')
print('##        /               |                |           \      ##')
print('##       R                B                V           bias   ##')
print('##      _|___            _|___            _|__           |    ##')
print('##     /     \          /     \          /    \          |    ##')
print('##   light   flat    light   flat      light  flat   trimmed  ##')
print('##    |       |        |       |         |      |             ##')
print('## reduced reduced  reduced reduced  reduced reduced          ##')
print('##                                                           ###')
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

##User input field for directory selection if desired
##(comment out if using user inuput)
#os.mkdir(os.path.dirname(__file__)+'/data/R/light')

#r_dir = input('Please input directory for the red filter lights:')
#r_flat_dir = input('Please input directory for the red filter flats:')
#v_dir = input('Please input directory for the green filter lights:')
#v_flat_dir = input('Please input directory for the green filter flats:')
#b_dir = input('Please input directory for the blue filter lights:')
#b_flat_dir = input('Please input directory for the blue filter flats:')
#bias_dir = input('Please input directory for bias frames:')

##Direct input field for directory selection 
##(comment out if using user inuput)
r_dir = os.path.dirname(__file__)+'/R/light'
r_flat_dir = os.path.dirname(__file__)+'/R/flat'
v_dir = os.path.dirname(__file__)+'/V/light'
v_flat_dir = os.path.dirname(__file__)+'/V/flat'
b_dir = os.path.dirname(__file__)+'/B/light'
b_flat_dir = os.path.dirname(__file__)+'/B/flat'
bias_dir =os.path.dirname(__file__)+'/bias'
r_reduced_dir = os.path.dirname(__file__)+'/R/light/reduced'
v_reduced_dir = os.path.dirname(__file__)+'/V/light/reduced'
b_reduced_dir = os.path.dirname(__file__)+'/B/light/reduced'
bias_trimmed_dir = os.path.dirname(__file__)+'/bias/trimmed'
r_flat_reduced_dir = os.path.dirname(__file__)+'/R/flat/reduced'
v_flat_reduced_dir = os.path.dirname(__file__)+'/V/flat/reduced'
b_flat_reduced_dir = os.path.dirname(__file__)+'/B/flat/reduced'
print('File directories sucessuflly built...')

#Build file colections for:
#                          lights, flats, bias', overscan trimmed bias',
#                          reduced flats and reduced light images.
r_collection = ccdproc.ImageFileCollection(r_dir)
r_flat_collection = ccdproc.ImageFileCollection(r_flat_dir)
v_collection = ccdproc.ImageFileCollection(v_dir)
v_flat_collection = ccdproc.ImageFileCollection(v_flat_dir)
b_collection = ccdproc.ImageFileCollection(b_dir)
b_flat_collection = ccdproc.ImageFileCollection(b_flat_dir)
bias_collection = ccdproc.ImageFileCollection(bias_dir)
bias_trimmed_collection = ccdproc.ImageFileCollection(bias_trimmed_dir)
r_reduced_collection = ccdproc.ImageFileCollection(r_reduced_dir)
v_reduced_collection = ccdproc.ImageFileCollection(v_reduced_dir)
b_reduced_collection = ccdproc.ImageFileCollection(b_reduced_dir)
r_flat_reduced_collection = ccdproc.ImageFileCollection(r_flat_reduced_dir)
v_flat_reduced_collection = ccdproc.ImageFileCollection(v_flat_reduced_dir)
b_flat_reduced_collection = ccdproc.ImageFileCollection(b_flat_reduced_dir)



print('Image file collections sucessfully built...')

##collection summary input selection
#print('Please state which data directory collection you wish to examine:')
#print('Folder selection shoud be in the format of "r","g","b","bias" ')
#print('Header will contain File, Image Type, X Pixel Range, Y Pixel Range, Biassec and Trimsec')
#fs = input('Folder Selection: ')
#fs_collection = str(fs+'_collection')
#fs_dir = (fs+'_dir')
#fs_collection = ccdproc.ImageFileCollection(fs_dir)
#fs_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']


#r_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#r_flat_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#v_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#v_flat_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#b_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#b_flat_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']
#bias_collection.summary['file', 'imagetyp', 'filter2','naxis1','naxis2','biassec','trimsec']


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
print('Biass sucessfully stacked...')
#Save combined bias in bias directory
combined_bias.write(os.path.join(bias_dir, "combined_bias.fit"), overwrite=True)
print('Maseter bias sucessfully saved as: combined_bias.fits')

# Reduce red flats
for ccd, file_name in r_flat_collection.ccds(imagetyp='DOME FLAT',       # Just get the bias frames
                                             return_fname=True           # Provide the file name too.
                                             ):
    # Subtract the overscan
    ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, 3095:], median=True)
    # Trim the overscan
    ccd = ccdproc.trim_image(ccd[:, :3078])
    # Subtract the bias 
    ccd = ccdproc.subtract_bias(ccd, combined_bias)
    print('Red flat overscan trimmed and master bias subtracted sucessfully,')
    # Save the result    
    ccd.write(os.path.join(r_flat_dir, "reduced", "r_flat_reduced_"+file_name), overwrite=True)    
    print('Reduced red flat saved as: ',os.path.join(r_flat_dir, "reduced", "r_flat_reduced_"+file_name))

print('Red flats sucessfully reduced...')

# Reduce green / v flats
for ccd, file_name in v_flat_collection.ccds(imagetyp='DOME FLAT',       # Just get the bias frames
                                             return_fname=True           # Provide the file name too.
                                             ):
    ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, 3095:], median=True)
    ccd = ccdproc.trim_image(ccd[:, :3078])
    ccd = ccdproc.subtract_bias(ccd, combined_bias)
    print('Green flat overscan trimmed and master bias subtracted sucessfully,')
    ccd.write(os.path.join(v_flat_dir, "reduced" , "v_flat_reduced_"+file_name), overwrite=True)
    print('Reduced green flat saved as: ',os.path.join(v_flat_dir, "reduced", "v_flat_reduced_"+file_name))

print('Green flats sucessfully reduced...')

# Reduce blue flats
for ccd, file_name in b_flat_collection.ccds(imagetyp='DOME FLAT',       # Just get the bias frames
                                             return_fname=True           # Provide the file name too.
                                             ):        
    ccd = ccdproc.subtract_overscan(ccd, overscan=ccd[:, 3095:], median=True)
    ccd = ccdproc.trim_image(ccd[:, :3078])
    ccd = ccdproc.subtract_bias(ccd, combined_bias)
    print('Blue flat overscan trimmed and master bias subtracted sucessfully,')
    ccd.write(os.path.join(b_flat_dir, "reduced" , "b_flat_reduced_"+file_name), overwrite=True)
    print('Reduced blue flat saved as: ',os.path.join(b_flat_dir, "reduced", "b_flat_reduced_"+file_name))

print('Blue flats sucessfully reduced...')




r_flat_reduced_dir = os.path.dirname(__file__)+'/R/flat/reduced'
v_flat_reduced_dir = os.path.dirname(__file__)+'/V/flat/reduced'
b_flat_reduced_dir = os.path.dirname(__file__)+'/B/flat/reduced'
r_flat_reduced_collection = ccdproc.ImageFileCollection(r_flat_reduced_dir)
v_flat_reduced_collection = ccdproc.ImageFileCollection(v_flat_reduced_dir)
b_flat_reduced_collection = ccdproc.ImageFileCollection(b_flat_reduced_dir)
print('Reduced flat directories and image file collections sucessfully built...')


def inv_median(a):
    return 1 / np.median(a)

#Stack red flats
r_to_combine = r_flat_reduced_collection.files_filtered(imagetyp='DOME FLAT', include_path=True)
r_combined_flat = ccdproc.combine(r_to_combine,
                                method='average', 
                                scale= inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

r_combined_flat.meta['combined'] = True
print('Red flats stacked sucessfully,')
r_combined_flat.write(os.path.join(r_flat_reduced_dir, "stacked_r_flat_reduced.fits"), overwrite=True)
print('Red master flat saved as: stacked_r_flat_reduced.fits')

# Stack green / v flats
v_to_combine = v_flat_reduced_collection.files_filtered(imagetyp='DOME FLAT', include_path=True)
v_combined_flat = ccdproc.combine(v_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

v_combined_flat.meta['combined'] = True
print('Green flats stacked sucessfully,')
v_combined_flat.write(os.path.join(v_flat_reduced_dir, "stacked_v_flat_reduced.fits"), overwrite=True)
print('Green master flat saved as: stacked_v_flat_reduced.fits')


# Stack blue flats
b_to_combine = b_flat_reduced_collection.files_filtered(imagetyp='DOME FLAT', include_path=True)
b_combined_flat = ccdproc.combine(b_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

b_combined_flat.meta['combined'] = True
print('Blue flats stacked sucessfully,')
b_combined_flat.write(os.path.join(b_flat_reduced_dir, "stacked_b_flat_reduced.fits"), overwrite=True)
print('Blue master flat saved as: stacked_b_flat_reduced.fits')


# Reduce red images and save
all_reds = []
red_light_ccds = []
for light, file_name in r_collection.ccds(imagetyp='OBJECT', return_fname=True):
    red_light_ccds.append(light)

    r_reduced = ccdproc.trim_image(light[:, :3078])
    r_reduced = ccdproc.subtract_bias(r_reduced, combined_bias)
    r_reduced = ccdproc.flat_correct(r_reduced, r_combined_flat)
    all_reds.append(r_reduced)
    print('Red light overscan trimmed, master bias and master red flat subtraced' )
    r_reduced.write(os.path.join(r_reduced_dir, "r_light_reduced_"+file_name), overwrite=True)
    print('Reduced red light saved as: ',os.path.join(r_reduced_dir, "r_light_reduced_"+file_name))

print('Red light sucessfully reduced...')

# Reduce green / v images and save
all_greens = []
green_light_ccds = []
for light, file_name in v_collection.ccds(imagetyp='OBJECT', return_fname=True):
    green_light_ccds.append(light)

    v_reduced = ccdproc.trim_image(light[:, :3078])
    v_reduced = ccdproc.subtract_bias(v_reduced, combined_bias)
    v_reduced = ccdproc.flat_correct(v_reduced, v_combined_flat)
    all_greens.append(v_reduced)
    print('Green light overscan trimmed, master bias and master green flat subtraced')
    v_reduced.write(os.path.join(v_reduced_dir, "v_light_reduced_"+file_name), overwrite=True)
    print('Reduced green light saved as: ',os.path.join(v_reduced_dir, "v_light_reduced_"+file_name))

print('Green light sucessfully reduced...')


# Reduce blue images and save
all_blues = []
blue_light_ccds = []
for light, file_name in b_collection.ccds(imagetyp='OBJECT', return_fname=True):
    blue_light_ccds.append(light)

    b_reduced = ccdproc.trim_image(light[:, :3078])
    b_reduced = ccdproc.subtract_bias(b_reduced, combined_bias)
    b_reduced = ccdproc.flat_correct(b_reduced, b_combined_flat)
    all_blues.append(b_reduced)
    print('Blue light overscan trimmed, master bias and master blue flat subtraced')
    b_reduced.write(os.path.join(b_reduced_dir, "b_light_reduced_"+file_name), overwrite=True)
    print('Reduced blue light saved as: ',os.path.join(b_reduced_dir, "b_light_reduced_"+file_name))

print('Blue light sucessfully reduced...')


# Stack red images and save
r_reduced_dir = os.path.dirname(__file__)+'/R/light/reduced'
r_reduced_collection = ccdproc.ImageFileCollection(r_reduced_dir)
print('Reduced red light directory and image file collection sucessfully built...')

red_light_to_combine = r_reduced_collection.files_filtered(imagetyp="OBJECT", include_path=True)
red_combined_light = ccdproc.combine(red_light_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

red_combined_light.meta['combined'] = True
red_combined_light.write(os.path.join(r_reduced_dir, "stacked_red_light_reduced.fits"), overwrite=True)

print('Red light stacked sucessfully,')
print('Red master light saved as: stacked_red_light_reduced.fits')



# Stack green / v images and save
v_reduced_dir = os.path.dirname(__file__)+'/V/light/reduced'
v_reduced_collection = ccdproc.ImageFileCollection(v_reduced_dir)
print('Reduced green light directory and image file collection sucessfully built...')

green_light_to_combine = v_reduced_collection.files_filtered(imagetyp="OBJECT", include_path=True)
green_combined_light = ccdproc.combine(green_light_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

green_combined_light.meta['combined'] = True
green_combined_light.write(os.path.join(v_reduced_dir, "stacked_green_light_reduced.fits"), overwrite=True)
print('Green light stacked sucessfully,')
print('Green master light saved as: stacked_green_light_reduced.fits')


# Stack blue images and save
b_reduced_dir = os.path.dirname(__file__)+'/B/light/reduced'
b_reduced_collection = ccdproc.ImageFileCollection(b_reduced_dir)
print('Reduced blue light directory and image file collection sucessfully built...')

blue_light_to_combine = b_reduced_collection.files_filtered(imagetyp="OBJECT", include_path=True)
blue_combined_light = ccdproc.combine(blue_light_to_combine,
                                method='average', 
                                scale=inv_median,
                                sigma_clip=True, 
                                sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5,
                                sigma_clip_func=np.ma.median, 
                                signma_clip_dev_func=mad_std,
                                mem_limit=350e6
                                )

blue_combined_light.meta['combined'] = True
blue_combined_light.write(os.path.join(b_reduced_dir, "stacked_blue_light_reduced.fits"), overwrite=True)
print('Blue light stacked sucessfully,')
print('Blue master light saved as: stacked_blue_light_reduced.fits')

# Confermation of science being done
for x in range(10):
   print('Science!')
