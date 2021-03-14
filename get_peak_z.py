#Written by Amy Gottlieb 11/2020

# This python script goes through each spot detected by source extractor and performs aperture photometry on its location in all images in a z stack.
# If plotted with good enough sampling, a plot of flux vs z slice should show a gaussian shape. 
# This script w


import numpy as np
from astropy.io import fits
# %matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
import pickle
import pandas as pd
import os
import sys
import time
import sep
from functools import reduce
from photutils import CircularAperture, CircularAnnulus

import track_cells_params as params
from track_cells_fns import *

print('Running get_peak_z.py')
overwrite=False
# This python script takes in an nd2 file 

#To run from python command line, syntax is:
# python get_spots.py /blue/eiken/Tcell_Motility/test.nd2

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
print('')

sep_folder=path+'sep_files_max/'

fits_folder=path+'fits_files/'
pk_folder=path+'peak_z_pk_files/'
if os.path.isdir(pk_folder)==False:
    print('Creating directory ',pk_folder)
    os.makedirs(pk_folder)
else:
    print('Directory already exists')#- removing current directory and creating a new one.')

    
with ND2Reader(full_path) as images:
    channels=images.metadata['channels']
    print('This nd2 file has the following channels:',channels)
    
    z_levels=range(images.sizes['z'])
    n_frames=range(images.sizes['t'])
    fov=images.metadata['fields_of_view']
    print('z',z_levels,'t',n_frames,'fov',fov)
    
flux_lim=params.flux_lim
sig_lim=params.sig_lim
n_points=params.n_points
r1=params.r1
r2=params.r2
r3=params.r3

time0=time.time()
for use_channel in channels:
    print('Channel',use_channel)
    time1=time.time()

    for use_fov in fov:
        print('   FOV',use_fov)
        time2=time.time()
        all_t_peak_z=[]
        all_zfilter=[]
        all_zinfo=[]
        old_txt=np.array(sorted([file for file in os.listdir(sep_folder) if file.endswith('.txt') 
                                 and ('fov'+str(use_fov) in file[-40:] )
                                 and (use_channel in file[-50:])]))
        txt=sort_files(old_txt,'_t','.')
        
        #need to load in all spots because they don't come with the attribute/column used flag 
        all_spots = np.array([pd.read_csv(sep_folder+f,sep='\t',lineterminator='\n') 
                                   for f in txt],dtype='object')
        for t in n_frames:
            time3=time.time()
            if t%20==0:
                print('      t',t)
                
            #From max project 
            x_cur=np.array(all_spots[t]['x']) #NOT sorted by flux ***** IS THIS A LIST????
            y_cur=np.array(all_spots[t]['y'])
            
            old_fnames=np.array(sorted([file for file in os.listdir(fits_folder) if 
                                        ('_t'+str(t)+'_' in file) 
                                        and ('fov'+str(use_fov) in file[-40:]) 
                                        and (use_channel in file[-40:])]))
            
            fits_files=sort_files(old_fnames,'_z','.') #z0,z1,z2...
            
            print('Getting fluxes of all spots in all z slices')
            #get fluxes of all spots in all z slices, so 1 spot has 41 z slices
            zstack=[]
            time4=time.time()
            for z in z_levels:

                #load in image
                hdu=fits.open(fits_folder+fits_files[z])
                data=hdu[0].data
                hdu.close()

                #subtract background
                data=data.byteswap(inplace=True).newbyteorder()    
                bkg = sep.Background(data)
                data_sub=data-bkg

                #do aperture photometry on current x,y position
                flux=do_ap_phot(data_sub,x_cur,y_cur,r1,r2,r3)
        
                zstack.append(flux)
            zstack=np.array(zstack)
            time5=time.time()
            print('This took ',(time5-time4)/60.,'min? for',z_levels[0],'z levels; t =',t)
#             print(len(zstack)) #41
#             print(len(zstack[0])) #195 all spots in first z slice
#             print(len(zstack[:,0])) #41 1st spot in all z slice

#             print('time for all z slices',time.time()-time4)

            print('Getting peak z')
            #find out where the spots peak in z space
            all_max=np.max(zstack,axis=0) #195
            peak_z=np.array([np.where(zstack[:,i]==all_max[i])[0][0] for i in range(len(zstack[0]))])
#             print(peak_z) #195
            zfilter=[]
            
            print('Calculating z filter for ',zstack.shape[1],'spots')
            for spot in range(zstack.shape[1]):
                lc= zstack[:,spot]
                sorted_flux=sorted(lc)[::-1]
                peak_flux_ind=np.where(np.array(lc)==np.max(lc))[0]
                flux50ind=sorted_flux[int(len(sorted_flux)/2)]  #throw out brightest half of points  ###############################
                bkg=lc[np.where(lc<flux50ind)[0]]
                bkg_std=np.std(bkg)
                above=np.where(lc>sig_lim*bkg_std)[0] ###############################
                print('above',above)
                n_conseq=findLongestConseqSubseq(above, len(above))
                conseq_ranges=np.array(consecutiveRanges(above, len(above)),dtype=object)
#                 print(len(conseq_ranges),conseq_ranges)
                flat_ranges= np.array([item for sublist in conseq_ranges for item in sublist])
                print('peak_flux_ind',peak_flux_ind[0],', conseq ranges (flattened)',flat_ranges)
#                 print('n_conseq',n_conseq,'n_points',n_points,'flux',flux,'flux_lim',flux_lim)
                if n_conseq>=n_points and (peak_flux_ind[0] in flat_ranges)==True and (lc>flux_lim).any()==True: #was orig flux>flux_lim
                    zfilter.append(1)
                    print('Passed zfiltering')
                else:
                    zfilter.append(0)
                    print('Failed zfiltering')
                    
#                 zinfo.append(list(lc))
            
#             all_zfilter.append(np.array(zfilter))
#             all_t_peak_z.append(peak_z)
#             all_zinfo.append(zstack)
            print('out of',zstack.shape[1],'spots',np.where(np.array(zfilter)==1)[0],'spots passed z filtering')
            print('Writing out info for t = ',t)
            col_names=["z"+str(i) for i in range(zstack.shape[0])]
            df_zinfo = pd.DataFrame(data=zstack.T, index=np.arange(zstack.shape[1]), columns=col_names)
            df_zinfo['peak_z']=peak_z
            df_zinfo['zfilter']=zfilter
            
            fname=pre+use_channel+'_fov'+str(use_fov)+'_t'+str(t)+'_peakz_zfilter_zinfo.csv'
            print(fname)
            df_zinfo.to_csv(path+'peak_z_pk_files/'+fname)
            print('-------------------------------------------------')
        #dump peak z info to pickle file
#         fname=pre+use_channel+'_fov'+str(use_fov)#+'_peakz_zfilter.pk'
#         pickle.dump( [all_t_peak_z,zfilter], open(path+'peak_z_pk_files/'+fname+'_peakz_zfilter.pk', "wb" ) )
#         pickle.dump( zfilter, open(path+'peak_z_pk_files/'+fname+'_zfilter.pk', "wb" ) )
#         pickle.dump( [zinfo], open(path+'peak_z_pk_files/'+fname+'_zinfo.pk', "wb" ) )

        print('Fov = ',use_fov,'; time for all t in 1 fov:',(time.time()-time2)/60.,'min')

    print('channel = ',use_channel,'; time for all t and all fov in one channel:',(time.time()-time1)/60.,'min')

print('total time',(time.time()-time0)/60.,'min')
print('Done running get_peak_z.py')