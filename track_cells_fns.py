import numpy as np
from astropy.io import fits
# %matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
# from nd2reader import ND2Reader
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


def calc_all_dist(loc_list, prev_loc,zdist=False):
    
    #prev_loc=np.array([x,y,z])
    #loc_list=np.array([np.array([x,y,z]),np.array([x,y,z])...]
    
    dx = loc_list[:,0] - prev_loc[0]
    dy = loc_list[:,1] - prev_loc[1]
#     dz = loc_list[:,2] - prev_loc[2]
    dist0=dx*dx + dy*dy #+ dz*dz
    if zdist==True:
        dist0+=(dz*dz)
    dist=dist0**0.5
    return dist

def calc_min_dist(loc_list, prev_loc,zdist=False):
    
    #prev_loc=np.array([x,y,z])
    #loc_list=np.array([np.array([x,y,z]),np.array([x,y,z])...]
#     print('prev_loc',prev_loc)
#     print(loc_list[:,0]-prev_loc[0])
    
    dx = loc_list[:,0] - prev_loc[0]
#     print('min_dx',loc_list[np.argmin(abs(dx))])
    
    dy = loc_list[:,1] - prev_loc[1]
#     print('min_dy',loc_list[np.argmin(abs(dy))])
    
#     dz = loc_list[:,2] - prev_loc[2]
#     print('min_dz',loc_list[np.argmin(abs(dz))])
    
    dist0=dx*dx + dy*dy #+ dz*dz
    if zdist==True:
        dist0+=(dz*dz)
    dist=dist0**0.5
    
    try:
        min_dist_ind=np.argmin(dist)
        min_dist=dist[min_dist_ind]
    except ValueError:
        min_dist_ind=np.nan
        min_dist=np.nan
        
    return min_dist,min_dist_ind
    
def Repeat(x): 
    _size = len(x) 
    repeated = [] 
    for i in range(_size): 
        k = i + 1
        for j in range(k, _size): 
            if x[i] == x[j] and x[i] not in repeated: 
                repeated.append(x[i]) 
    return repeated 

def add_columns(df_list):
    print('Adding used flag, time, distance, radius, compare_id, peak z, and zfilter')
    count=0
    for time1,spot_list in enumerate(df_list):
        #set used_flag attribute to 0 (not used) for all spots
        spot_list['used']=np.zeros(len(spot_list),dtype=np.int8)
        spot_list['time']=np.ones(len(spot_list),dtype=np.int8)*int(time1)
        spot_list['distance']=np.zeros(len(spot_list),dtype=np.int8)
        spot_list['radius']=np.zeros(len(spot_list))
#         a = np.empty(len(spot_list))
#         a[:] = np.nan
#         spot_list['peak_z']=a
#         spot_list['zfilter']=a

        compare_id=[]
        for spot in range(len(spot_list)):
            compare_id.append(int(count))
            count+=1
        spot_list['compare_id']=compare_id
        
    return df_list,count

def do_ap_phot(data_sub,x,y,r1,r2,r3):
    #data should already be bkg subtracted

    positions=merge(x,y)

    target_apertures = CircularAperture(positions, r=r1)
    sky_annuli = CircularAnnulus(positions, r_in=r2, r_out=r3)

    target_masks = target_apertures.to_mask(method='subpixel', subpixels=5)
    flux = np.array([np.sum(el.multiply(data_sub)) for el in target_masks])
    
    return flux

def sort_files(old_txt,zt,punc):
    nums=np.array([int(f.split(zt)[-1].split(punc)[0]) for f in old_txt]).argsort()
    txt=old_txt[nums] #t0,t1,t2....
    return txt

def merge(list1, list2): 
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))] 
    return merged_list 

def create_df(col_names,df_old):
    data={}
    for col in col_names:
        data[col]=[df_old[col]]
    df_new = pd.DataFrame(data)
    return df_new

def write_spots(df_list,num_times,pre,track_folder,use_channel,v):
    #df_list is list of dataframes of each spot,
    #length of df_list is total number of spots found
    
#     num_times=4 #spot must be in at least 4 time slices
    lengths=[]
    for k in range(len(df_list)): # is total number of spots in all time slices
        length=len(df_list[k])
        lengths.append(length)
        compare_id=int(np.array(df_list[k]['compare_id'])[0])
        spot_name=pre+use_channel+'_fov'+str(v)+'_spot'+str(k)+'_id'+str(compare_id)+'_'+str(len(df_list[k]))+'slices.csv'

        #if the spot was found in more than 4 time slices, write all of the info to a csv file
        if length>num_times:
            print('Writing to ',track_folder+spot_name)
            df_list[k].to_csv(track_folder+spot_name)
        else:
            print('Less than',num_times,'time slices for spot #',k)
    print(lengths)

    
def write_fits(path,name,image):
    if os.path.isfile(path+name)==False:
        fits.writeto(path+name,image)
        print( 'Written: ',path+name)
        to_continue=False
    else:
        print('File already exists- ',path+name)
        to_continue=True
        
    return to_continue


def findLongestConseqSubseq(arr, n): 
    '''We insert all the array elements into unordered set.'''
  
    S = set(); 
    for i in range(n): 
        S.add(arr[i]); 
  
    # check each possible sequence from the start 
    # then update optimal length 
    ans = 0; 
    for i in range(n): 
          
        # if current element is the starting 
        # element of a sequence 
        if S.__contains__(arr[i]): 
              
            # Then check for next elements in the 
            # sequence 
            j = arr[i]; 
              
            # increment the value of array element 
            # and repeat search in the set 
            while(S.__contains__(j)): 
                j += 1; 
  
            # Update optimal length if this length 
            # is more. To get the length as it is 
            # incremented one by one 
            ans = max(ans, j - arr[i]); 
    return ans; 
 
    

def write_cat_to_txt(objects,fname):
    with open(fname, 'w') as file:
        for label in objects.dtype.names:
            file.write(label+'\t')
        file.write('\n')
        if len(objects)>0:
            for l in objects:
                for ele in l:
                    file.write(str(ele) + '\t')
                file.write('\n')
#         else:
#             for label in objects.dtype.names:
#                 file.write('Nan \t')
#             file.write('\n')
            
    
def consecutiveRanges(a, n):
 
    length = 1
    arr = []
     
    # If the array is empty,     # return the list 
    if (n == 0):
        return np.array(arr)
     
    # Traverse the array     # from first position 
    for i in range (1, n + 1):
     
        # Check the difference # between the current # and the previous elements 
        # If the difference doesn't # equal to 1 just increment  the length variable. 
        if (i == n or a[i] - a[i - 1] != 1):
        
            # If the range contains # only one element. # add it into the list. 
            if (length == 1):
#                 list.append(str(a[i - length]))
#                 arr.append(np.array([a[i-length]]))
                pass
            else:
     
                # Build the range between the first # element of the range and the 
                # current previous element as the # last range. 
#                 temp = (str(a[i - length]) + " -> " + str(a[i - 1]))          
                temp=np.arange(a[i-length],a[i-1]+1)
                arr.append(temp)
           
            # After finding the  # first range initialize # the length by 1 to # build the next range. 
            length = 1
        
        else:
            length += 1
    return np.array(arr,dtype=object)


def find_brightest_spot(df_next,all_dist,d_lim,prev_flux):
    closest_ind=np.arange(len(df_next))[all_dist<d_lim]
    if len(closest_ind)>0:
        #get the close by spots
        closest=df_next[all_dist<d_lim] #dataframe

        #get the fluxes of the close by spots
        closest_fluxes=np.array(df_next['flux'][all_dist<d_lim]) #[15000]

        #find the spot with the largest flux, use that as match
        max_ind=closest_ind[np.argmax(closest_fluxes)] #0
        max_ind2=np.argmax(closest_fluxes)

         #check the flux difference between current spot and matched spot in next time slice
        f_diff=abs(prev_flux-closest_fluxes)
    #                                 f_diff_percent=f_diff/prev_flux*100.
        f_diff_percent=closest_fluxes/prev_flux*100.

        # #                                 print('closest (dataframe)',closest)
    # #                                 print('closest fluxes',closest_fluxes,'max_ind',max_ind,'min_dist_ind',min_dist_ind)
    # #                                 print(df_next.loc[max_ind,'flux'], all_spots[new_t].loc[inds[min_dist_ind],'flux'])
    # #                                 print('previous flux',prev_flux)
        print('   Fluxes',closest_fluxes, 'f_diff_percent',f_diff_percent,'; max flux:',closest_fluxes[max_ind2])

        return closest_fluxes,max_ind,max_ind2,f_diff,f_diff_percent
    else: 
        print('No spots found within increased search radius of',d_lim,'pixels')
        return [],np.nan,np.nan,np.nan,np.nan
                                