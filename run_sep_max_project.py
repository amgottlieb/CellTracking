#Written by Amy Gottlieb 9/8/2020

# This python script takes in an nd2 file, all possible fits images in the corresponding folder, runs source extractor on it and writes the catalog to a txt file.

#To run from python command line, syntax is:
# python run_sep.py /blue/eiken/Tcell_Motility/test.nd2

import sys
import numpy as np
from astropy.io import fits
import sep
import os

overwrite=True
add_noise=False

print('In run_sep.py')
print('Overwrite = ',overwrite)
print('Add noise = ',add_noise)

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
            

#FULL PATH MUST HAVE .ND2 
#full_path='/blue/eiken/Tcell_Motility/test.nd2'            
full_path=sys.argv[1]
split_str=full_path.split('/')

#pre='/blue/eiken/Tcell_Motility/'
pre=full_path[:-1*len(split_str[-1])]

#name='test/'
name=split_str[-1][:-4]+'/'

#path='/blue/eiken/Tcell_Motility/test/'
path=pre+name


try_fits_folder1=path+'max_project/'

if os.path.isdir(try_fits_folder1):
    fits_folder=try_fits_folder1
    sep_folder=path+'sep_files_max/'
    use_ind=-12 ## z0_zstack.fits

print(fits_folder)

#if there is not already a folder corresponding to this nd2 file, create one
if os.path.isdir(sep_folder)==False:
    print('Creating directory ',sep_folder)
    os.makedirs(sep_folder)

fnames=np.array(sorted([file for file in os.listdir(fits_folder) if file.endswith(".fits")]))

for i in range(len(fnames)):
    txt=sep_folder+fnames[i][:use_ind]+'.txt'
    
    #if a txt file does not already exist:
    if os.path.isfile(txt)==False or overwrite==True:
        print('Working on ',fits_folder+fnames[i])
        
        ###################################################################
        #open the image and get the data
        hdu=fits.open(fits_folder+fnames[i])
        data=hdu[0].data
        data=data.byteswap(inplace=True).newbyteorder()    
    
        ###ADD NOISE HERE IF NECESSARY###
        if add_noise==True:
            shape=data.shape
            noise=np.random.rand(shape[0], shape[1])*2.
            data=data*1.0+noise
        
        #get the background and subtract it from the data
        bkg = sep.Background(data)
        data_sub=data-bkg
        
        print('Average of data:',np.mean(data))
        print('Average of background:',np.mean(bkg))
        print('Average of background subtracted data:',np.mean(data_sub))
        print('Global RMS',np.mean(bkg.globalrms))
        print('')
        
        #get a catalog of objects with that are 3sigma above the background from the background subtracted image using source extractor
        #NOTE: Sometimes 3 isn't high enough and it will throw an exception b/c it found one giant object. 
        #So, start with 3 and keep increasing until it doesn't throw an error; 
        #see https://github.com/kbarbary/sep/issues/15  for a better explanation
        thresh=3.
#         objects = sep.extract(data_sub, thresh, err=bkg.globalrms)

        while True:
            print('Threshold:',thresh)
            try:
                print('Threshold*GlobalRMS',thresh*np.mean(bkg.globalrms))
                objects = sep.extract(data_sub, thresh, err=bkg.globalrms)#,filter_kernel=None,clean=False) #,deblend_cont=0.1) #,deblend_nthresh=8)
            except Exception:
                thresh+=1
                continue
            break
        
        ####################################################################
        
        
        f=fnames[i][:use_ind]+'.txt'
        split=f.split('_')
        print('Split',split)
        #find index that has fov, stick threshold info in filename before that 
        for i,string in enumerate(split):
            if 'fov' in string:
#                 print('Found fov at ind',i)
                #If index is 3, the element is inserted after the 3rd element. Its position will be 4th.
                split.insert(i,'SEPthresh'+str(int(thresh)))
                print('new split',split)
                joined = '_'.join(split)
                new_txt=sep_folder+joined
#                 print('new txt file name',new_txt)
                break
#         new_txt=txt
        if os.path.isfile(new_txt)==False:
            #write the catalog to a txt file
            write_cat_to_txt(objects,new_txt)     
            print('Wrote',len(objects),'objects to ',new_txt)            
            hdu.close()
        else:
            print('File already exists- ',new_txt)
        print('')
        
    else:
        print('File already exists- ',txt)
        continue
        