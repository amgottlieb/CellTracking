#Written by Amy Gottlieb 9/3/2020

# This python script takes in an nd2 file, gets the metadata, loops through all possible images and writes each image to a fits file.

#To run from python command line, syntax is:
# python extract_fits_from_nd2_general.py /blue/eiken/Tcell_Motility/test.nd2

import numpy as np
from nd2reader import ND2Reader
from astropy.io import fits
import os.path
import sys

#if the fits file does not already exist, write it to the file. Otherwise, skip it.
def write_fits(path,name,image):
    if os.path.isfile(path+name)==False:
        fits.writeto(path+name,image)
        print( 'Written: ',path+name)
        to_continue=False
    else:
        print('File already exists- ',path+name)
        to_continue=True
        
    return to_continue

# def rebin(a, shape):
#     sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
#     return a.reshape(sh).mean(-1).mean(1)

# def image_zsmooth2(imstack1):
#     imstack2 = np.zeros(imstack1.shape)
#     for i in range(imstack2.shape[2]-1):
#         imstack2[:,:,i] = (imstack1[:,:,i] + imstack1[:,:,i+1])/2.0
        
#     imstack2[:,:,imstack2.shape[2]-1] = imstack2[:,:,imstack2.shape[2]-2]
    
#     return imstack2

#full_path='/blue/eiken/Tcell_Motility/test.nd2'
full_path=sys.argv[1] 
# print('full path',full_path)

split_str=full_path.split('/')
# print('Split_str',split_str)

#path='/blue/eiken/Tcell_Motility'
path=full_path[:-1*len(split_str[-1])]
print('path',path)

#name = 'test.nd2'
name=split_str[-1]
print('name',name)

#pre = 'test'
pre=name[:-4]
print('pre',pre)

#folder='/blue/eiken/Tcell_Motility/test/'
folder=path+pre+'/'

#if there is not already a folder corresponding to this nd2 file, create one
if os.path.isdir(folder)==False:
    print('Creating directory ',folder)
    os.makedirs(folder)
else:
    print('Directory already exists:',folder)
    
fits_folder=folder+'fits_files/'

#if there is not already a folder for fits files in the folder corresponding to this nd2 file, create one
if os.path.isdir(fits_folder)==False:
    print('Creating directory ',fits_folder)
    os.makedirs(fits_folder)
else:
    print('Directory already exists:',fits_folder)

print('IMAGE',path+name)

with ND2Reader(path+name) as images:
    
    print( images.metadata )
    print( images.sizes )

    #Get metadata
    n_frames=images.metadata['num_frames'] #time slices
    if type(n_frames)==int:
        n_frames=range(n_frames)
    channels=images.metadata['channels']
    fov=images.metadata['fields_of_view']
    z_levels=images.metadata['z_levels']
    
    print(len(channels),len(fov),len(n_frames),len(z_levels))
    
    images.bundle_axes='yx'
    
    #iterate over time or z_level depending on how the data were taken; not sure if this actually does anything since I loop over both later anyways
    
    if len(z_levels)==0:# and n_frames!=0:
        print('Iterating over t')
        images.iter_axes='t'
    else:
        print('Iterating over z')
        images.iter_axes='z'

    #loop through all channels, fields of view, time slices and z_levels

    for c in range(len(channels)):
        #if there is only one channel, code breaks if you try to set default_coords
        if len(channels)>1:
            images.default_coords['c']=c
        for v in fov:
            #same for fields of view
            if len(fov)>1:
                images.default_coords['v']=v
            for t in n_frames:
                images.default_coords['t']=t
                
                # if there is more than one z level, loop through them here and write out each image. 
                if len(z_levels)!=0:
                    imf=np.zeros([np.int(images.frame_shape[0]/1.0),np.int(images.frame_shape[1]/1.0),images.sizes['z']])
                    im_bin = np.zeros([256,256,imf.shape[2]])
                    for z in z_levels:
                        images.default_coords['z']=z
                        
                        base=pre+'_'+channels[c]+'_fov'+str(v)+'_t'+str(t)+'_z'+str(z)
                        fits_name=base+'.fits'
                        to_continue=write_fits(fits_folder,fits_name,images[z]*1.0)
                        
#                         sep_name=base+'_SEP.fits'
#                         bin_name=base+'_binned_SEP.fits'
                        #sometimes background in an image is very low so add minimum noise for source extractor to latch onto
                        #so we can still use a reasonable threshold.
                        #NOTE: These images should ONLY be used for source extractor. Use the normal fits file for everything else.
#                         shape=images[z].shape
#                         noise=np.random.rand(shape[0], shape[1])*2.
#                         sep_arr=images[z]*1.0+noise
#                         to_continue=write_fits(sep_folder,sep_name,sep_arr)
                               
#                         imf[:,:,z]=images[z]*1.0
#                         im_bin[:,:,z] = rebin(imf[:,:,z],[256,256])

#                     #binning and smoothing images to increase S/N
#                     im_smooth = image_zsmooth2(im_bin)
#                     for zi in range(im_smooth.shape[2]):
#                         bin_name=pre+'_'+channels[c]+'_fov'+str(v)+'_t'+str(t)+'_z'+str(zi)+'_binned_SEP.fits'
#                         to_continue=write_fits(bin_folder,bin_name,im_smooth[:,:,zi])
    
                # Otherwise, just write out the image.            
                else:
                    
#                     sep_name=pre+'_'+channels[c]+'_fov'+str(v)+'_t'+str(t)+'_SEP.fits'
                    fits_name=pre+'_'+channels[c]+'_fov'+str(v)+'_t'+str(t)+'.fits'
                    
                    to_continue=write_fits(fits_folder,fits_name,images[t]*1.0)
                    if to_continue==True:
                        pass
                        
#                     shape=images[t].shape
#                     noise=np.random.rand(shape[0], shape[1])*2.
#                     sep_arr=images[t]*1.0+noise
                        
#                     to_continue=write_fits(sep_folder,sep_name,sep_arr)
#                     if to_continue==True:
#                         continue

