#Written by Amy Gottlieb 9/3/2020

# This python script takes in an nd2 file, gets the metadata, loops through all possible images and writes each image to a fits file.

#To run from the python command line, the syntax is:
# python extract_fits_from_nd2.py /path/to/test.nd2

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
        
    else:
        print('File already exists- ',path+name)
        
    return 


#full_path='/path/to/test.nd2'
full_path=sys.argv[1] 

split_str=full_path.split('/')

#path='/path/to'
path=full_path[:-1*len(split_str[-1])]
print('path',path)

#name = 'test.nd2'
name=split_str[-1]
print('name',name)

#pre = 'test'
pre=name[:-4]
print('pre',pre)

#folder='/path/to/test/'
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
    
    if len(z_levels)==0:#
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
            
            #if there is only one field of view, code breaks if you try to set default_coords
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
                        write_fits(fits_folder,fits_name,images[z]*1.0)
                        
                # Otherwise, just write out the image.            
                else:
                    
                    fits_name=pre+'_'+channels[c]+'_fov'+str(v)+'_t'+str(t)+'.fits'                    
                    write_fits(fits_folder,fits_name,images[t]*1.0)


