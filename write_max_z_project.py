#Written by Amy Gottlieb 11/2020

# This python script takes in an nd2 file, gets all of the z slices for a particular time/channel/fov, 
# takes the maximum of all the z slices, and writes that image to a file

import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
import pickle
import pandas as pd
import os
import sys
import csv
from functools import reduce
import time
import shutil

def write_fits(path,name,image,overwrite):
    if os.path.isfile(path+name)==False or overwrite==True:
        fits.writeto(path+name,image,overwrite=overwrite)
        print( 'Written: ',path+name)
        to_continue=False
    else:
        print('File already exists- ',path+name)
        to_continue=True
        
    return to_continue


#full_path='/blue/eiken/Tcell_motility/TCellMot_LLSsizeTest_Full-Small-Med-Large_FITC_10x_2020SEP06.nd2'
full_path=sys.argv[1] 

# path='/blue/eiken/Tcell_motility/test'+'/'
path=full_path[:-4]+'/'
print(path)

split_str=full_path.split('/')

#name = 'test.nd2'
name=split_str[-1]

# pre='test_'
pre=name[:-4]+'_'
print(pre)

with ND2Reader(full_path) as images:
    print(images.sizes)
#     print(images.metadata)
    channels=images.metadata['channels']
    z_levels=range(images.sizes['z'])
    n_frames=range(images.sizes['t'])
    fov=range(len(images.metadata['fields_of_view']))
    print(channels,'z',z_levels,'t',n_frames,'fov',fov)

check_folder=path+'max_project/'

if os.path.isdir(check_folder)==False:
    print('Creating directory ',check_folder)
    os.makedirs(check_folder)
else:
    print('Directory already exists:',check_folder)

    overwrite=False

for use_channel in channels:
    for use_fov in fov:
        for use_t in n_frames:

            old_fits_files=np.array(sorted([file for file in os.listdir(path+'fits_files/') if 
                                        (use_channel in file[-40:]) and 
                                        ('fov'+str(use_fov) in file[-30:]) and 
                                        ('t'+str(use_t)+'_' in file[-30:])]))

            nums=np.array([int(f.split('_z')[-1].split('.')[0]) for f in old_fits_files]).argsort()
            fits_files=old_fits_files[nums] #t0,t1,t2....

            fname=fits_files[0].split('_z')[0]+'_zstack.fits'

            all_z_slices=[]
            for use_z in z_levels:
                hdu=fits.open(path+'fits_files/'+fits_files[use_z])
                data=hdu[0].data
                data=data.byteswap(inplace=True).newbyteorder() 
                hdu.close()

                all_z_slices.append(data)
            all_z_slices=np.array(all_z_slices)

            use_data=np.max(all_z_slices,axis=0)
           
            to_continue=write_fits(check_folder,fname,use_data,overwrite)
