import os
import imageio.v2 as imageio
import cv2
from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
from collections import defaultdict
from cellpose import models, io, utils
from cellpose.io import imread
import time
from scipy.ndimage import find_objects
import argparse
import itertools


def remove_edges(image, percent =0.05, wipe_value=0):
    '''
    Change the values of pixels 5% from the edge of an image to 0./Crop the image 5% from the edges but keep the origianl image size. 
    Parameters
    image: single channel image
    percent: the percentage from the edge to crop
    wipe_value: the new pixel value.'''
    height, width = image.shape[:2]
    left = int(width * percent)
    top = int(height * percent)
    right = int(width * (1 - percent))
    bottom = int(height * (1 - percent))
    
    wiped_image = image.copy()
    wiped_image[:top, :] = wipe_value
    wiped_image[bottom:, :] = wipe_value
    wiped_image[:, :left] = wipe_value
    wiped_image[:, right:] = wipe_value
    
    return wiped_image

def remove_edge_masks(masks, change_index=True, percent = 0.05):
    """ remove masks with pixels on edge of image
    * copied from cellpose source code, and modified to remove masks with pixels that fall in the 5% margin from the edges.
    
    Parameters
    ----------------

    masks: int, 2D or 3D array 
        size [Ly x Lx] or [Lz x Ly x Lx], 0=NO masks; 1,2,...=mask labels

    change_index: bool (optional, default True)
        if True, after removing masks change indexing so no missing label numbers
    
    percent:
    Returns
    ----------------

    outlines: 2D or 3D array 
        size [Ly x Lx] or [Lz x Ly x Lx], 0=NO masks; 1,2,...=mask labels

    """
    slices = find_objects(masks.astype(int))
    height, width = masks.shape[0], masks.shape[1]
    border_pixels = int(percent * min(height, width))
    for i,si in enumerate(slices):
        remove = False
        if si is not None:
            for d,sid in enumerate(si):
                if sid.start < border_pixels or sid.stop > masks.shape[d] - border_pixels:
                    remove=True
                    break  
            if remove:
                masks[si][masks[si]==i+1] = 0
    shape = masks.shape
    if change_index:
        _,masks = np.unique(masks, return_inverse=True)
        masks = np.reshape(masks, shape).astype(np.int32)

    return masks

def load_image(image_paths: List[str],input_folder):
    '''Load the single channel images and stack them together.
    Parameters
    -----------------------
    image_paths: a list of files of the images that are to be stacked.
    input_folder: the directory that the images are in.

    Returns
    ------------------------
    the stacked image (nY x nX x channels).'''
    image_stack = []
    for path in image_paths: 
        image_part = imageio.imread(os.path.join(input_folder, path))
        #image_part = remove_edges(image_part) # not needed with remove_edge_masks
        image_stack.append(image_part)
    return np.stack(image_stack, axis=-1)

def load_folder(input_path):
    '''Read all image files in the input folder and load every image set.
    Parameters
    -----------------------
    input_path(str): path to the input folder
    Returns
    -----------------------
    stack(dict): a dictionary of stacked images with the file names of the corresponding image set as keys.'''
    filenames = []
    images = []
    file_list = os.listdir(input_path)
    file_list = [file for file in file_list if file.endswith('.tiff') and 'None' not in file] #get rid of index files and empty wells
    file_list.sort() #in alphabetical order ch1-5
    for iteration in range(0,len(file_list),5):
        current_image = file_list[iteration:(iteration+5)]  # get the images of 5 channels
        #print(current_image) 
        filename = file_list[iteration].split("_ch")[0] #common file name
        filenames.append(filename)
        current_image = load_image(current_image,input_path)
        images.append(current_image)
        #print(f"Image set added: {filename}")
    stack = dict(zip(filenames,images))
    print(f"Folder loaded: {input_path}.")
    return stack

def get_segmentation_masks(imgs, model_sepcs = {"model_type":"nuclei","channels":[[1, 0]],"diameter":None,"flow_threshold": 0.4}):
    '''Run the cellpose model on the images and get the object masks (edge masks removed).
    Parameters
    --------------------------
    imgs (lst): the input image sets 
    model_specs (dict): settings for model type, channels to segment, object diameter, flow threshold
    '''
    #[1,0] for nuclei
    #[3,1] for cells
    model = models.Cellpose(gpu=True, model_type = model_sepcs["model_type"])
    masks, flows, styles, diams = model.eval(
        imgs,
        diameter=model_sepcs["diameter"],
        channels=model_sepcs["channels"],
        flow_threshold=model_sepcs["flow_threshold"],)
    for mask in masks:
        mask = remove_edge_masks(mask)  # remove edge masks
    print(f"Segmentation completed for: {model_sepcs['model_type']} on channels {model_sepcs['channels']} with diameter {model_sepcs['diameter']} and flow threshold {model_sepcs['flow_threshold']}.")
    return masks

def get_object_size(mask, object_idx):
    mn = mask==object_idx
    return mn.sum()

def get_object_intensities(img,obj_mask,channel_idx):
    '''Get the pixel intensities of the selected channel within an object (normalized to the range [0,1])
    Parameters
    -------------------------------
    img: the stacked image
    obj_mask: mask of a single object. False=no mask, True=object
    channel_idx: the channel to extract intensities from
    Returns
    -------------------------------
    pixel_values_normalized: an array of pixel intensities within the object'''
    #channel_dict = {"DAPI":0,"Marker":1,"Tubulin":2,"HPA_H":3,"HPA_L":4}
    #np.set_printoptions(linewidth=np.inf)
    pixel_values = img[:,:,channel_idx][obj_mask]
    pixel_values_normalized = pixel_values/65535
    return pixel_values_normalized

def get_object_outline(mask, object_idx):
    '''copied from the source code for utils.outlines_list().
    Get the outline of a mask with a given object label.'''
    mn = mask==object_idx
    contours = cv2.findContours(mn.astype(np.uint8), mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE)
    contours = contours[-2]
    cmax = np.argmax([c.shape[0] for c in contours])
    pix = contours[cmax].astype(int).squeeze()
    return pix if len(pix) > 4 else np.zeros((0, 2))

def get_seg_data(filename,img,output_folder,model_sepcs = {"model_type":"nuclei","channels":[[1, 0]],"diameter":None,"flow_threshold": 0.4}):
    '''Save the segmentation masks for each image set, and save measurements to a dataframe.
    Parameters
    ----------------------------------
    filename:
    img: image set
    output_folder: the directory to save the segmentation masks in
    model_specs: settings for segmentation
    Returns
    ----------------------------------
    seg_data: a dataframe of all objects segmented in this run, with the filename, outline, location, size and mask for individual object.'''
    masks = get_segmentation_masks([img],model_sepcs=model_sepcs)
    # Save masks and measurements
    objects_data = []
    model_type = model_sepcs["model_type"]
    for mask in masks:# loop over every image set
        # Save the masks as numpy array and as png
        np.save(os.path.join(output_folder,f"{filename}_{model_type}_mask.npy"),mask)
        io.imsave(os.path.join(output_folder,f"{filename}_{model_type}_mask.png"),mask)
        # Make a dataframe for the measurements
        if model_type == "nuclei":
            for n in np.unique(mask)[1:]:# loop over every object within the image
                outline = get_object_outline(mask,n)
                centroid = outline.mean(axis=0)
                object_data = {
                    "Filename":filename,
                    #"Outline":outline,
                    "Location_Center_X":centroid[0],
                    "Location_Center_Y":centroid[1],
                    "Size": get_object_size(mask,n),
                    #"Object_Mask": mask == n,
                    "Object_Label": n,
                }
                objects_data.append(object_data)
        else:
            for n in np.unique(mask)[1:]:# loop over every object within the image
                outline = get_object_outline(mask,n)
                #centroid = outline.mean(axis=0)
                object_data = {
                    "Filename":filename,
                    "Outline":outline,
                    #"Location_Center_X":centroid[0],
                    #"Location_Center_Y":centroid[1],
                    "Size": get_object_size(mask,n),
                    #"Object_Mask": mask == n,
                    "Object_Label": n,
                }
                objects_data.append(object_data)
    seg_data = pd.DataFrame(objects_data)
    del objects_data
    print(f"Object data saved for {model_type}.")
    return seg_data

def get_object_mask(img_mask, obj_label):
    mask = np.load(img_mask)
    obj_mask = (mask == obj_label)
    return obj_mask

def get_cyto_data(nuc_locations: pd.DataFrame, cyto_outlines: pd.DataFrame, img_stack, output_folder):
    '''Match Cell and nuclei objects from nuclei locations and cell(cyto) outlines, and create final output dataframes of nuclei, cell, and cytosol with corresponding cell IDs.
    Parameters
    --------------------------------------
    nuc_locations: the dataframe output from segmenting nuclei objects
    cyto_outlines: the dataframe output from segmenting cyto(whole cell) objects
    img_stack: the dictionary of input image sets
    Returns
    --------------------------------------
    object_data: list of dict of the matched cell/nuclei/cytosol objects, with information on file name/cell ID/location and measurements of size and intensities.'''
    print("Matching cell and nuclei objects.")

    # Create a list to store the entries
    object_data = []

    # Create dictionary to store nuclei locations based on filename
    nuc_locations_dict  = defaultdict(list)
    for nuc_index, nuc_row in nuc_locations.iterrows():
        nuc_locations_dict[nuc_row["Filename"]].append(nuc_row)
    del nuc_locations
    # add Cell_ID column
    cyto_outlines = cyto_outlines.reset_index()
    cyto_outlines = cyto_outlines.rename(columns={"index":"Cell_ID"})

    # iterate over all cytoplasm outlines
    for cyto_index, cyto_row in cyto_outlines.iterrows():
        cell_id = cyto_row["Cell_ID"]

        #create path from outlines
        cytoplasm_path = mplPath.Path(cyto_row["Outline"])

        #check if nuclei are within the path
        filename = cyto_row["Filename"]
        nuc_locations_subset = nuc_locations_dict.get(filename,[])
        for nuc_row in nuc_locations_subset:
            nuc_center = (nuc_row["Location_Center_X"], nuc_row["Location_Center_Y"])
            # if nuclei is in cytoplasm path, give it the same cell ID as cytoplasm that it is within
            if cytoplasm_path.contains_point(nuc_center):
                #nuc_mask = nuc_row["Object_Mask"]
                nuc_mask = get_object_mask(os.path.join(output_folder,f"{filename}_nuclei_mask.npy"),nuc_row["Object_Label"])
                #cell_mask = cyto_row["Object_Mask"]
                cell_mask = get_object_mask(os.path.join(output_folder,f"{filename}_cyto_mask.npy"),cyto_row["Object_Label"])
                #cyto_mask = cell_mask ^ nuc_mask 
                cyto_mask = (cell_mask & ~nuc_mask)
                #centroid = cyto_row["Outline"].mean(axis=0)# convert cytoplasm outlines to center coords
                current_entry = {
                    "Filename":filename,
                    "Cell_ID":cell_id,
                    "Parent_Cell":cyto_row["Object_Label"],
                    "Parent_Nuc":nuc_row["Object_Label"],
                    "Nuc_Size":nuc_row["Size"],
                    "Cell_Size":cyto_row["Size"],
                    "Cyto_Size":cyto_row["Size"] - nuc_row["Size"],
                    "Cyto_Size_mask":cyto_mask.sum(),
                    "Nuc_Intensities_DAPI":get_object_intensities(img_stack[filename],nuc_mask,channel_idx=0),
                    "Nuc_Intensities_Marker":get_object_intensities(img_stack[filename],nuc_mask,channel_idx=1),
                    "Nuc_Intensities_Tubulin":get_object_intensities(img_stack[filename],nuc_mask,channel_idx=2),
                    "Nuc_Intensities_HPA_H":get_object_intensities(img_stack[filename],nuc_mask,channel_idx=3),
                    "Nuc_Intensities_HPA_L":get_object_intensities(img_stack[filename],nuc_mask,channel_idx=4),
                    "Cell_Intensities_DAPI":get_object_intensities(img_stack[filename],cell_mask,channel_idx=0),
                    "Cell_Intensities_Marker":get_object_intensities(img_stack[filename],cell_mask,channel_idx=1),
                    "Cell_Intensities_Tubulin":get_object_intensities(img_stack[filename],cell_mask,channel_idx=2),
                    "Cell_Intensities_HPA_H":get_object_intensities(img_stack[filename],cell_mask,channel_idx=3),
                    "Cell_Intensities_HPA_L":get_object_intensities(img_stack[filename],cell_mask,channel_idx=4),
                    "Cyto_Intensities_DAPI":get_object_intensities(img_stack[filename],cyto_mask,channel_idx=0),
                    "Cyto_Intensities_Marker":get_object_intensities(img_stack[filename],cyto_mask,channel_idx=1),
                    "Cyto_Intensities_Tubulin":get_object_intensities(img_stack[filename],cyto_mask,channel_idx=2),
                    "Cyto_Intensities_HPA_H":get_object_intensities(img_stack[filename],cyto_mask,channel_idx=3),
                    "Cyto_Intensities_HPA_L":get_object_intensities(img_stack[filename],cyto_mask,channel_idx=4),
                    #"Nuc_Outline"::nuc_row["Outline"].tolist(),
                    "Nuc_Location_Center_X":nuc_row["Location_Center_X"],
                    "Nuc_Location_Center_Y":nuc_row["Location_Center_Y"],
                    #"Cell_Outline":cyto_row["Outline"].tolist(),
                    "Cell_Location_Center_X":cyto_row["Outline"].mean(axis=0)[0],
                    "Cell_Location_Center_Y":cyto_row["Outline"].mean(axis=0)[1],
                }
                object_data.append(current_entry)
                del nuc_mask,cell_mask,cyto_mask
        cyto_outlines.drop(cyto_index, inplace=True)
    del nuc_locations_dict,cyto_outlines
    # create dataframe for object data using dicts created within the function
    #object_data = pd.DataFrame.from_dict(object_data)

    #print("Dataframes are created with matched cell IDs.")
    # cyto_data = cyto_data.drop_duplicates(subset = "Cell_ID")
    print("Matched object data saved.\n")
    return object_data

def chunked(it, size):
    it = iter(it)
    while True:
        p = tuple(itertools.islice(it,size))
        if not p:
            break
        yield p

def seg_folder(data_path, save_path, nuc_model_specs, cyto_model_specs):
    '''Segment the images from an input folder and save the outputs.'''
    start = time.time()
    
    os.makedirs(save_path, exist_ok=True)
    
    # Loading images
    image_stack = load_folder(data_path)
    end_load_folder = time.time()
    print(f"Loading images:{end_load_folder-start:.2f} seconds.\n")

    # Segment nuclei and whole cells and match to get object info
    start_segmentation = time.time()

    batch = 0
    #output_files = []
    for chunk in chunked(image_stack.items(), 100):
        batch += 1
        start_chunk = time.time()
        object_data = []
        print(f"Processing: batch {batch}")
        for filename, image in chunk:
            print(f"Processing image: {filename}")
            nuc_data_unlabeled = get_seg_data(filename,image,output_folder=save_path,model_sepcs=nuc_model_specs)
            cell_data_unlabeled = get_seg_data(filename,image,output_folder=save_path,model_sepcs=cyto_model_specs)
            object_data_labeled = get_cyto_data(nuc_data_unlabeled,cell_data_unlabeled,image_stack,output_folder=save_path)
            object_data.extend(object_data_labeled)
        end_chunk_segmentation = time.time()
        print(f"Processing of batch {batch} completed in: {end_chunk_segmentation-start_chunk:.2f} seconds.")
        
        start_save= time.time()
        object_data = pd.DataFrame.from_dict(object_data)

        plate_name = os.path.basename(save_path)
        output_file = f"{plate_name}_batch{batch}.feather"
        #output_files.append(output_file)
        object_data.to_feather(os.path.join(save_path,output_file))
        del object_data
        end_save = time.time()
        print(f"Output file saved as {output_file}")
        print(f"Saving output files of batch {batch}: {end_save-start_save:.2f} seconds.\n")
    '''
    # Concatenate all outputs
    print("Concatenating output files...")
    del image_stack
    start_concat = time.time()
    dfs = []
    for file in output_files:
        df = pd.read_feather(os.path.join(save_path,file))
        dfs.append(df)
    final_output = pd.concat(dfs,ignore_index = True)
    final_output.to_feather(os.path.join(save_path,f"{os.path.basename(save_path)}.feather"))
    '''
    end = time.time()
    print(f"Processing of {os.path.basename(save_path)} completed in {batch} batches in {end-start_segmentation:.2f} seconds.\n")
    #print(f"Processing output files:{end-start_concat:.2f} seconds.\n")
    print(f"Loading images: {end_load_folder-start:.2f} seconds.\n")
    print(f"Time in total: {end-start:.2f} seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Segment images using CellPose")
    parser.add_argument("data_path", type=str, help="Path to input data folder")
    parser.add_argument("save_path", type=str, help="Path to output data folder")
    args = parser.parse_args()

    nuc_sepcs = {"model_type": "nuclei", "channels": [[1, 0]], "diameter": None, "flow_threshold": 0.4}
    cell_specs = {"model_type": "cyto", "channels": [[3, 1]], "diameter": None, "flow_threshold": 0.4}

    seg_folder(data_path=args.data_path, save_path=args.save_path, nuc_model_specs=nuc_sepcs, cyto_model_specs=cell_specs)
