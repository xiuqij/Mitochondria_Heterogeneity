#%%
import os
import pandas as pd

'''This script renames the image file acquired from Opera system, adding information on the sample for easier analysis.
Example of Opera image file name: r01c01f01p01-ch1sk1fk1fl1.tiff
Example of renamed image file: TUFM_1_A01_HPA000286_f01_ch1.tiff
'''
#%% Change directory
#directory = os.getcwd()
newdirectory = '' # Enter directory where image files are saved

#directory = directory + '/filerenametest'
#print(directory)
os.chdir(newdirectory)
directory = os.getcwd()
print(directory)
#%% 
plate_layout = '' # Enter the path for the file with Antibody/well information
plate_codes = '' # Enter the path for the file with plate barcode information

# Make dictionary for row names
rowdict = {'r01':'A','r02':'B', 'r03':'C', 'r04':'D', 'r05':'E', 'r06':'F', 'r07':'G', 'r08':'H'}   

# Make dictionary for plate name and its barcode
pl_code = pd.read_csv(plate_codes)
expdict = dict(zip(pl_code.code,pl_code.plate)) 

#fov_dict = {'f01':'f11','f02':'f01','f03':'f02','f04':'f03','f05':'f04','f06':'f05','f07':'f06','f08':'f07','f09':'f08','f10':'f09','f11':'f10'}

# Make dictionary for antibody in each well
pl_map = pd.read_csv(plate_layout)  
pl_map['plate_well'] = ''
for i in pl_map.index:
    pl_map.plate_well[i] = '{}_{}'.format(pl_map.plate[i],pl_map.well[i])
Abdict = dict(zip(pl_map.plate_well,pl_map.Ab))
#print(Abdict)
#%% Renaming
for subdir, dirs, files in os.walk(directory):
    for file in files:
        filepath = subdir + os.sep + file
        #print(filepath)
        if filepath.endswith(".tiff"):
            #print(filepath)
            filepathsplit = filepath.split(os.sep)
            oldfilename = filepathsplit[-1] #r01c01f01p01-ch1sk1fk1fl1.tiff
            
            #Get plate name
            folder = filepathsplit[-3] #P017151__2023-03-21T15_33_45-Measurement 1
            plate_id = folder.split("_")[0]
            plate_name = expdict.get(plate_id)
            plate_number = plate_name.split('_')[1]
            
            #Get well name
            filenamesplit = oldfilename.split("-")
            row = rowdict.get(filenamesplit[0][0:3])
            col = filenamesplit[0][4:6]
            well = row + col
            
            #Get the HPA antibody in this well
            plate_well = '{}_{}'.format(plate_number,well)
            Ab = Abdict.get(plate_well)

            #Get FOV
            fov = filenamesplit[0][6:9]
            #if fov in fov_dict:
                #fov = fov_dict[fov]
            #Get channel
            channel = filenamesplit[1][0:3]
            
            
            print('\nOld filename: ', oldfilename)

            newfilename =  "{}_{}_{}_{}_{}.tiff".format(plate_name,well,Ab,fov,channel)
            print('New filename: ', newfilename)
            newfilepath = subdir + os.sep + newfilename
            print("The new filepath plus file name is: ", newfilepath)

            os.rename(filepath, newfilepath)
            
