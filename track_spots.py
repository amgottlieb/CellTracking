#Written by Amy Gottlieb 11/2020

# This python script ****


import numpy as np
from astropy.io import fits
# %matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
from nd2reader import ND2Reader
import pickle
import pandas as pd
import os
import csv
from functools import reduce
import time
import sys
import shutil
import photutils
from photutils import CircularAperture, CircularAnnulus, aperture_photometry

import track_cells_params as params
from track_cells_fns import *

# flux_lim=400.#1221. #spots with fluxes below this limit will not be analyzed
# d_lim=20. #spots with distances greater than this limit will not be counted as a potential match/ 
# #brightest spot within this limit will be considered the match
# num_times=4 #if a spot appears less than this limit/# of time slices, it will not be written to a file
# num_empty=4 #if a spot doesn't appear in the next 4 time slices, stop looking for it
flux_lim=params.flux_lim
flux_drop_lim1=params.flux_drop_lim1
flux_drop_lim2=params.flux_drop_lim2
d_lim=params.d_lim
num_times=params.num_times
num_empty=params.num_empty
factor=params.factor

#r1 = (xm1 + ym1 + xym1)**0.5


print(sys.argv)
full_path=sys.argv[1]
# full_path='/blue/eiken/Tcell_motility/TCellMot_LLSsizeTest_Full-Small-Med-Large_FITC_10x_2020SEP06.nd2'

# d_lim=float(sys.argv[2])

# path='/blue/eiken/Tcell_motility/test'+'/'
path=full_path[:-4]+'/'
print(path)

split_str=full_path.split('/')

#name = 'test.nd2'
name=split_str[-1]

# pre='test_'
pre=name[:-4]+'_'
print(pre)

with ND2Reader(path[:-1]+'.nd2') as images:
    print(images.sizes)
#     print(images.metadata)
    channels=images.metadata['channels']
    z_levels=range(images.sizes['z'])
    n_frames=range(images.sizes['t'])
    fov=range(len(images.metadata['fields_of_view']))
    print(channels,'z',z_levels,'t',n_frames,'fov',fov)
    
track_folder=path+'track_spot_csv_files/'
if os.path.isdir(track_folder)==False:
    print('Creating directory ',track_folder)
    os.makedirs(track_folder)
else:
    print('Directory already exists')#- removing current directory and creating a new one.')
#     shutil.rmtree(track_folder)
#     os.makedirs(track_folder)


fits_folder=path+'max_project/'
sep_folder=path+'sep_files_max/'
to_plot_zstack=False
if len(n_frames)>1:
    #loop through all channels and fovs
    for c,use_channel in enumerate(channels):
        for use_fov in fov:
            #load in source extractor txt files of max project images and sort them by time
            old_txt=np.array(sorted([file for file in os.listdir(sep_folder) if file.endswith('.txt') 
                                     and ('fov'+str(use_fov) in file[-40:] )
                                     and (use_channel in file[-50:])]))
            txt=sort_files(old_txt,'_t','.')

            
            print('loading in all spots, adding columns')
            #need to load in all spots because they don't come with the attribute/column used flag 
            all_spots_init = np.array([pd.read_csv(sep_folder+f,sep='\t',lineterminator='\n') 
                                       for f in txt],dtype='object')

            #add used flag, time, distance, peak_z, zfilter, and compare_id columns to data frame
            all_spots,count = add_columns(all_spots_init)

#             #load in peak_z info
#             pk_name=pre+use_channel+'_fov'+str(use_fov)+'_peakz_zfilter.pk'
#             all_t_peak_z,all_zfilter,all_zinfo=pickle.load(open(path+'peak_z_pk_files/'+pk_name, "rb" ) )

#             for t in n_frames:
#                 #add peak z and zfilter info from file into dataframe for the current time
#                 all_spots[t]['peak_z']=all_t_peak_z[t]
#                 all_spots[t]['zfilter']=all_zfilter[t]
#                 all_spots[t]['zinfo']=list(all_zinfo[t])

            #load in peak z, zfilter, zinfo
            for t in n_frames:
                fname=pre+use_channel+'_fov'+str(use_fov)+'_t'+str(t)+'_peakz_zfilter_zinfo.csv'
                df_zinfo=pd.read_csv(path+'peak_z_pk_files/'+fname)
                all_spots[t]=pd.concat([all_spots[t], df_zinfo], axis=1)
            
            print('Looping through time slices....')
            new_spot_list=[]
            time0=time.time()
            for t in n_frames:

                if t%10 == 0:
                    print('Time: ',t)

                if t == 0:
                    #if on the first time slice, set all distances to 0
                    all_spots[t]['distance']=np.zeros(len(all_spots[t]['distance']))
                
                #get the x,y values of spots in the current time slice t
                x_cur=np.array(all_spots[t]['x']) #NOT sorted by flux
                y_cur=np.array(all_spots[t]['y'])

                #now sort spots by flux so that brightest spots are 1st
    #                 print('SORTING BY FLUX')
                sorted_spots=all_spots[t].sort_values(by='flux',inplace=False,ascending=False)
                index_order=np.array(sorted_spots.index)
                print('Index order',index_order)
                col_names=["z"+str(i) for i in z_levels]
                
                print('Looping through spots.............................')
                #start looping through all spots, brightest first
                for i,spot_num in enumerate(index_order):
                    #if the spot is brighter than the limit
                    lc=all_spots[t].loc[spot_num,col_names]
                    if max(lc)>flux_lim or all_spots[t].loc[spot_num,'flux']>flux_lim:
#                     if all_spots[t].loc[spot_num,'flux']>flux_lim: #CHANGE ME
                        
#                         if all_spots[t].loc[spot_num,'zfilter']==1:

                        print('t=',t)

                        compare_id=all_spots[t].loc[spot_num,'compare_id']
                        print('i,spot_num,compare_id',i,spot_num,compare_id)
                
                        #check if the spot has already been matched to another time slice/used
                        if all_spots[t].loc[spot_num,'used']==-1:
                            print('This spot already used- skipping')
                            continue

                        #set the current spot to used
                        all_spots[t].loc[spot_num,'used']=-1

                        #calculate the radius of the spot
                        all_spots[t].loc[spot_num,'radius']= np.sqrt(all_spots[t].loc[spot_num,'x2']+ all_spots[t].loc[spot_num,'y2']+ all_spots[t].loc[spot_num,'xy'])

                        #create a new data frame for this spot
                        spot_df=create_df(all_spots[t].columns,all_spots[t].loc[spot_num])#.copy()

                        #get the x,y coordinates of this spot
                        xcen=all_spots[t].loc[spot_num,'x']
                        ycen=all_spots[t].loc[spot_num,'y']
                        print('xcen,ycen',xcen,ycen)
                        print('flux',all_spots[t].loc[spot_num,'flux'])

                        if to_plot_zstack==True:
                            plot_zstack(xcen,ycen,lc,std,spot_num)
#                         print('-------')

                        empty_count=0

                        #loop over all the remaining time slices (ex. if initial t = 2, start here at t=3)
                        for new_t in range(t+1,len(n_frames)):
                            if new_t==t+1:
                                print('First new t:',new_t)
#                             if new_t%20==0:
#                                 print('new t',new_t)

                            #get the flags for all the spots in next time slice 
                            flags=np.array(all_spots[new_t]['used'])

                            # get all the spots in the next time slice that haven't been used already
                            inds=np.array(all_spots[new_t][flags==0].index)
        #                     print('number of spots that havent been used:',len(inds))
                            df_next=all_spots[new_t][flags==0]

                             # get list of positions (np.array(n,2)=np.array([np.array(x,y,),
                            pos1=np.array((df_next['x'],df_next['y']))
                            pos = pos1.T 

                            #calculate distances of all spots from current spots, find the minimum distance
                            loc=np.array([xcen,ycen])
                            all_dist=calc_all_dist(pos, loc)
                            min_dist,min_dist_ind=calc_min_dist(pos, loc)

                            #Check all spots within distance limit/get the indices of the close by spots in the next time slice
                            closest_ind=np.arange(len(df_next))[all_dist<d_lim] #[36]
                            print('Of ',len(inds),' spots that havent been used already in t =',new_t,',',len(closest_ind), 'are within',d_lim,'pixels' )
                            prev_flux=np.array(spot_df['flux'])[-1]
                            
                            #if there is a spot in the next time slice close enough to the current spot
                            if len(closest_ind)>0:

                                closest_fluxes,max_ind,max_ind2,f_diff,f_diff_percent=find_brightest_spot(df_next,all_dist,d_lim,prev_flux)
                                
                                #if ratio of fluxes<20 or >220, repeat with larger search radius to find spot
                                if f_diff_percent[max_ind2]<=flux_drop_lim1 or f_diff_percent[max_ind2]>=flux_drop_lim2:
                                    print('    Flux difference too large; need to find spot using larger search radius.')
                                    print('    fluxes so far:',np.array(spot_df['flux']))
                                    print('    prev_flux:',prev_flux,'closest fluxes:',closest_fluxes,
                                          'brightest',closest_fluxes[max_ind2])
                                    print('increasing d_lim by ',factor,'x= ',d_lim*factor)
                                
                                    closest_fluxes,max_ind,max_ind2,f_diff,f_diff_percent=find_brightest_spot(df_next,all_dist,d_lim*factor,prev_flux)    
                                    if f_diff_percent[max_ind2]<=flux_drop_lim1 or f_diff_percent[max_ind2]>=flux_drop_lim2:
                                        print('    couldnt find any spots after increasing search radius')
                                        empty_count+=1
                                        append_spot=False
                                    else:
                                        print('    f_diff_percent',f_diff_percent[max_ind2],'prev_flux',prev_flux,
                                          'new flux:',closest_fluxes[max_ind2])
                                        append_spot=True

                                else:
                                    print('    f_diff_percent',f_diff_percent[max_ind2],'prev_flux',prev_flux,
                                          'new flux:',closest_fluxes[max_ind2])
                                    append_spot=True
#                                     print('flux difference of ',f_diff_percent[max_ind2],'percent is okay')
                                
                            else: #if no spots were found within first radius
                                print('    no spots found within initial search radius of',d_lim,'pixels')
                                print('    fluxes so far:',np.array(spot_df['flux']))
                                print('    prev_flux',prev_flux)
                                print('    increasing d_lim to',factor,'*d_lim=',factor*d_lim)
                                
                                closest_fluxes,max_ind,max_ind2,f_diff,f_diff_percent=find_brightest_spot(df_next,all_dist,d_lim*factor,prev_flux)
                                
                                #if a spot is found by increasing the search radius, update it; otherwise, it wasn't found
                                if len(closest_fluxes)>0:
                                    if f_diff_percent[max_ind2]<=flux_drop_lim1 or f_diff_percent[max_ind2]>=flux_drop_lim2:
                                        print('    couldnt find any spots after increasing search radius')
                                        empty_count+=1
                                        append_spot=False
                                    else:
                                        print('    f_diff_percent',f_diff_percent[max_ind2],'prev_flux',prev_flux,
                                          'new flux:',closest_fluxes[max_ind2])
                                        append_spot=True
                                else:
                                    print('    couldnt find any spots after increasing search radius')
                                    empty_count+=1
                                    append_spot=False
                                #if there are no spots nearby in the next timeslice, right now do nothing
                                #check next time slice
#                                 print('Of ',len(inds),' spots that havent been used already, NO NEARBY SPOTS IN TIME SLICE ',new_t)
                                
                                #if the spot hasn't been found in 4 time slices in a row, stop looking for it
                                if empty_count==num_empty:
                                    print('spot not found in the last',num_empty,'slices; last new t',new_t)
                                    break
                                else:
                                    pass
                                
                            print('Append spot is',append_spot)
                            if append_spot==True:
                                radius=np.sqrt(all_spots[new_t].loc[inds[max_ind],'x2']+
                                all_spots[new_t].loc[inds[max_ind],'y2']+ 
                                all_spots[new_t].loc[inds[max_ind],'xy'])

                                #update distance, used_flag, compare_id,radius
                                all_spots[new_t].loc[inds[max_ind],'used']=-1
                                all_spots[new_t].loc[inds[max_ind],'compare_id']=compare_id
                                all_spots[new_t].loc[inds[max_ind],'distance']=min_dist
                                all_spots[new_t].loc[inds[max_ind],'radius']= radius 

                                #add matched spot to current spot dataframe; NEED PEAK_Z to be updated before this***
                                spot_df=spot_df.append(all_spots[new_t].loc[inds[max_ind]],ignore_index=False)

                                ###NEED TO UPDATE XCEN,YCEN so that you're comparing the distance from the previous
                                ##time slice, not always the first one
                                xcen=all_spots[new_t].loc[inds[max_ind],'x']
                                ycen=all_spots[new_t].loc[inds[max_ind],'y']
                                empty_count=0
                                
                        print('Appending new_spot list')
                        new_spot_list.append(spot_df)
                        print('-------------')
#                         else:
#                             print('Spot #',spot_num,'does not span >=3 zslices (zfilter=',all_spots[t].loc[spot_num,'zfilter'],
#                                   '; flux = ',all_spots[t].loc[spot_num,'flux'])

                    else:
                        print('Spot #',spot_num,' too faint',all_spots[t].loc[spot_num,'flux'])
                        #set the current spot to used
                        all_spots[t].loc[spot_num,'used']=-1
                time5=time.time()
                print('Done looping through all spots in t = ',t,)
                print('total for t=',t,':',time5-time0)

            write_spots(new_spot_list,num_times,pre,track_folder,use_channel,use_fov)

else:
    print('Only 1 time slice- cannot run track spots')

    

