#%% Import packages
import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from functools import reduce
from statannotations.Annotator import Annotator
from itertools import combinations
import math 

############# Reading and processing output from CellProfiler #####################
#%% Functions for finding information from the CellProfiler output file name
def find_exp(filename):
    #20230404_IDH3A_5_Nuclei.csv
    filenamesplit = filename.split('_')
    exp = '{}_{}'.format(filenamesplit[1],filenamesplit[2])
    return exp

def find_prefix(filename):
    if 'Cells.csv' in filename:
        prefix = 'Cell_'
    elif 'FilteredNuclei.csv' in filename:
        prefix = 'Nuc_'
    elif 'Cytoplasm.csv' in filename:
        prefix = 'Cyto_'
    return prefix
#%% Functions for getting information from FileName_DAPI 
def find_plate(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    plate = '{}_{}'.format(filenamesplit[0],filenamesplit[1])
    return plate
def find_well(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    well = filenamesplit[2]
    return well 
def find_well_id(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    well_id = '{}_{}_{}'.format(filenamesplit[0],filenamesplit[1],filenamesplit[2])
    return well_id
def find_marker(filenameDAPI):
    #IDH3A_5_A01_HPA021995_f01_ch1.tiff
    filenamesplit = filenameDAPI.split("_")
    #plate = '{}_{}'.format(filenamesplit[0],filenamesplit[1])
    #well = filenamesplit[2]
    marker = filenamesplit[0]
    return marker
def find_hpa_ab(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    hpa_ab = filenamesplit[3]
    return hpa_ab
def find_hpa_gene(filenameDAPI):
    plate_layout = 'C:/Users/ji_xi/OneDrive/Studies/Degree Project MTLS/Image analysis/plate_layout_correct.csv'
    pl_map = pd.read_csv(plate_layout)
    dict_HPA_genes = dict(zip(pl_map.Ab,pl_map.gene))
    filenamesplit = filenameDAPI.split("_")
    hpa_ab = filenamesplit[3]
    hpa_gene = dict_HPA_genes.get(hpa_ab)
    return hpa_gene
def find_fixation_method(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    marker = filenamesplit[0]
    if marker == 'IDH3A' or marker == 'SOD2':
        fixation = 'Methanol'
    if marker == 'GLS' or marker == 'TUFM':
        fixation = 'PFA'  
    return fixation
def find_plate_well(filenameDAPI):
    filenamesplit = filenameDAPI.split("_")
    platewell = '{}_{}'.format(filenamesplit[1],filenamesplit[2])
    return platewell


#%%
############################################# Functions for reading output #####################################
def read_output(dir):
    '''
    Read in the measurement output from CellProfiler and merge into one dataframe
    '''
    plate_codes = ''# Enter the path for the file with plate barcode information
    pl_code = pd.read_csv(plate_codes)
    plate_list = pl_code['plate'].to_list() # Get a list of all plate names
    #all_columns = ['ImageNumber', 'ObjectNumber', 'FileName_DAPI', 'FileName_HPA_H', 'FileName_HPA_L', 'FileName_marker', 'FileName_tubulin', 'PathName_DAPI', 'PathName_HPA_H', 'PathName_HPA_L', 'PathName_marker', 'PathName_tubulin', 'AreaShape_Area', 'AreaShape_BoundingBoxArea', 'AreaShape_BoundingBoxMaximum_X', 'AreaShape_BoundingBoxMaximum_Y', 'AreaShape_BoundingBoxMinimum_X', 'AreaShape_BoundingBoxMinimum_Y', 'AreaShape_Center_X', 'AreaShape_Center_Y', 'AreaShape_Compactness', 'AreaShape_ConvexArea', 'AreaShape_Eccentricity', 'AreaShape_EquivalentDiameter', 'AreaShape_EulerNumber', 'AreaShape_Extent', 'AreaShape_FormFactor', 'AreaShape_MajorAxisLength', 'AreaShape_MaxFeretDiameter', 'AreaShape_MaximumRadius', 'AreaShape_MeanRadius', 'AreaShape_MedianRadius', 'AreaShape_MinFeretDiameter', 'AreaShape_MinorAxisLength', 'AreaShape_Orientation', 'AreaShape_Perimeter', 'AreaShape_Solidity', 'Correlation_Correlation_CropDAPI_CropHPAH', 'Correlation_Correlation_CropDAPI_CropHPAL', 'Correlation_Correlation_CropDAPI_CropMarker', 'Correlation_Correlation_CropDAPI_CropTubulin', 'Correlation_Correlation_CropHPAH_CropHPAL', 'Correlation_Correlation_CropHPAH_CropMarker', 'Correlation_Correlation_CropHPAH_CropTubulin', 'Correlation_Correlation_CropHPAL_CropMarker', 'Correlation_Correlation_CropHPAL_CropTubulin', 'Correlation_Correlation_CropMarker_CropTubulin', 'Correlation_Costes_CropDAPI_CropHPAH', 'Correlation_Costes_CropDAPI_CropHPAL', 'Correlation_Costes_CropDAPI_CropMarker', 'Correlation_Costes_CropDAPI_CropTubulin', 'Correlation_Costes_CropHPAH_CropDAPI', 'Correlation_Costes_CropHPAH_CropHPAL', 'Correlation_Costes_CropHPAH_CropMarker', 'Correlation_Costes_CropHPAH_CropTubulin', 'Correlation_Costes_CropHPAL_CropDAPI', 'Correlation_Costes_CropHPAL_CropHPAH', 'Correlation_Costes_CropHPAL_CropMarker', 'Correlation_Costes_CropHPAL_CropTubulin', 'Correlation_Costes_CropMarker_CropDAPI', 'Correlation_Costes_CropMarker_CropHPAH', 'Correlation_Costes_CropMarker_CropHPAL', 'Correlation_Costes_CropMarker_CropTubulin', 'Correlation_Costes_CropTubulin_CropDAPI', 'Correlation_Costes_CropTubulin_CropHPAH', 'Correlation_Costes_CropTubulin_CropHPAL', 'Correlation_Costes_CropTubulin_CropMarker', 'Correlation_K_CropDAPI_CropHPAH', 'Correlation_K_CropDAPI_CropHPAL', 'Correlation_K_CropDAPI_CropMarker', 'Correlation_K_CropDAPI_CropTubulin', 'Correlation_K_CropHPAH_CropDAPI', 'Correlation_K_CropHPAH_CropHPAL', 'Correlation_K_CropHPAH_CropMarker', 'Correlation_K_CropHPAH_CropTubulin', 'Correlation_K_CropHPAL_CropDAPI', 'Correlation_K_CropHPAL_CropHPAH', 'Correlation_K_CropHPAL_CropMarker', 'Correlation_K_CropHPAL_CropTubulin', 'Correlation_K_CropMarker_CropDAPI', 'Correlation_K_CropMarker_CropHPAH', 'Correlation_K_CropMarker_CropHPAL', 'Correlation_K_CropMarker_CropTubulin', 'Correlation_K_CropTubulin_CropDAPI', 'Correlation_K_CropTubulin_CropHPAH', 'Correlation_K_CropTubulin_CropHPAL', 'Correlation_K_CropTubulin_CropMarker', 'Correlation_Manders_CropDAPI_CropHPAH', 'Correlation_Manders_CropDAPI_CropHPAL', 'Correlation_Manders_CropDAPI_CropMarker', 'Correlation_Manders_CropDAPI_CropTubulin', 'Correlation_Manders_CropHPAH_CropDAPI', 'Correlation_Manders_CropHPAH_CropHPAL', 'Correlation_Manders_CropHPAH_CropMarker', 'Correlation_Manders_CropHPAH_CropTubulin', 'Correlation_Manders_CropHPAL_CropDAPI', 'Correlation_Manders_CropHPAL_CropHPAH', 'Correlation_Manders_CropHPAL_CropMarker', 'Correlation_Manders_CropHPAL_CropTubulin', 'Correlation_Manders_CropMarker_CropDAPI', 'Correlation_Manders_CropMarker_CropHPAH', 'Correlation_Manders_CropMarker_CropHPAL', 'Correlation_Manders_CropMarker_CropTubulin', 'Correlation_Manders_CropTubulin_CropDAPI', 'Correlation_Manders_CropTubulin_CropHPAH', 'Correlation_Manders_CropTubulin_CropHPAL', 'Correlation_Manders_CropTubulin_CropMarker', 'Correlation_Overlap_CropDAPI_CropHPAH', 'Correlation_Overlap_CropDAPI_CropHPAL', 'Correlation_Overlap_CropDAPI_CropMarker', 'Correlation_Overlap_CropDAPI_CropTubulin', 'Correlation_Overlap_CropHPAH_CropHPAL', 'Correlation_Overlap_CropHPAH_CropMarker', 'Correlation_Overlap_CropHPAH_CropTubulin', 'Correlation_Overlap_CropHPAL_CropMarker', 'Correlation_Overlap_CropHPAL_CropTubulin', 'Correlation_Overlap_CropMarker_CropTubulin', 'Correlation_RWC_CropDAPI_CropHPAH', 'Correlation_RWC_CropDAPI_CropHPAL', 'Correlation_RWC_CropDAPI_CropMarker', 'Correlation_RWC_CropDAPI_CropTubulin', 'Correlation_RWC_CropHPAH_CropDAPI', 'Correlation_RWC_CropHPAH_CropHPAL', 'Correlation_RWC_CropHPAH_CropMarker', 'Correlation_RWC_CropHPAH_CropTubulin', 'Correlation_RWC_CropHPAL_CropDAPI', 'Correlation_RWC_CropHPAL_CropHPAH', 'Correlation_RWC_CropHPAL_CropMarker', 'Correlation_RWC_CropHPAL_CropTubulin', 'Correlation_RWC_CropMarker_CropDAPI', 'Correlation_RWC_CropMarker_CropHPAH', 'Correlation_RWC_CropMarker_CropHPAL', 'Correlation_RWC_CropMarker_CropTubulin', 'Correlation_RWC_CropTubulin_CropDAPI', 'Correlation_RWC_CropTubulin_CropHPAH', 'Correlation_RWC_CropTubulin_CropHPAL', 'Correlation_RWC_CropTubulin_CropMarker', 'Intensity_IntegratedIntensityEdge_CropDAPI', 'Intensity_IntegratedIntensityEdge_CropHPAH', 'Intensity_IntegratedIntensityEdge_CropHPAL', 'Intensity_IntegratedIntensityEdge_CropMarker', 'Intensity_IntegratedIntensityEdge_CropTubulin', 'Intensity_IntegratedIntensity_CropDAPI', 'Intensity_IntegratedIntensity_CropHPAH', 'Intensity_IntegratedIntensity_CropHPAL', 'Intensity_IntegratedIntensity_CropMarker', 'Intensity_IntegratedIntensity_CropTubulin', 'Intensity_LowerQuartileIntensity_CropDAPI', 'Intensity_LowerQuartileIntensity_CropHPAH', 'Intensity_LowerQuartileIntensity_CropHPAL', 'Intensity_LowerQuartileIntensity_CropMarker', 'Intensity_LowerQuartileIntensity_CropTubulin', 'Intensity_MADIntensity_CropDAPI', 'Intensity_MADIntensity_CropHPAH', 'Intensity_MADIntensity_CropHPAL', 'Intensity_MADIntensity_CropMarker', 'Intensity_MADIntensity_CropTubulin', 'Intensity_MassDisplacement_CropDAPI', 'Intensity_MassDisplacement_CropHPAH', 'Intensity_MassDisplacement_CropHPAL', 'Intensity_MassDisplacement_CropMarker', 'Intensity_MassDisplacement_CropTubulin', 'Intensity_MaxIntensityEdge_CropDAPI', 'Intensity_MaxIntensityEdge_CropHPAH', 'Intensity_MaxIntensityEdge_CropHPAL', 'Intensity_MaxIntensityEdge_CropMarker', 'Intensity_MaxIntensityEdge_CropTubulin', 'Intensity_MaxIntensity_CropDAPI', 'Intensity_MaxIntensity_CropHPAH', 'Intensity_MaxIntensity_CropHPAL', 'Intensity_MaxIntensity_CropMarker', 'Intensity_MaxIntensity_CropTubulin', 'Intensity_MeanIntensityEdge_CropDAPI', 'Intensity_MeanIntensityEdge_CropHPAH', 'Intensity_MeanIntensityEdge_CropHPAL', 'Intensity_MeanIntensityEdge_CropMarker', 'Intensity_MeanIntensityEdge_CropTubulin', 'Intensity_MeanIntensity_CropDAPI', 'Intensity_MeanIntensity_CropHPAH', 'Intensity_MeanIntensity_CropHPAL', 'Intensity_MeanIntensity_CropMarker', 'Intensity_MeanIntensity_CropTubulin', 'Intensity_MedianIntensity_CropDAPI', 'Intensity_MedianIntensity_CropHPAH', 'Intensity_MedianIntensity_CropHPAL', 'Intensity_MedianIntensity_CropMarker', 'Intensity_MedianIntensity_CropTubulin', 'Intensity_MinIntensityEdge_CropDAPI', 'Intensity_MinIntensityEdge_CropHPAH', 'Intensity_MinIntensityEdge_CropHPAL', 'Intensity_MinIntensityEdge_CropMarker', 'Intensity_MinIntensityEdge_CropTubulin', 'Intensity_MinIntensity_CropDAPI', 'Intensity_MinIntensity_CropHPAH', 'Intensity_MinIntensity_CropHPAL', 'Intensity_MinIntensity_CropMarker', 'Intensity_MinIntensity_CropTubulin', 'Intensity_StdIntensityEdge_CropDAPI', 'Intensity_StdIntensityEdge_CropHPAH', 'Intensity_StdIntensityEdge_CropHPAL', 'Intensity_StdIntensityEdge_CropMarker', 'Intensity_StdIntensityEdge_CropTubulin', 'Intensity_StdIntensity_CropDAPI', 'Intensity_StdIntensity_CropHPAH', 'Intensity_StdIntensity_CropHPAL', 'Intensity_StdIntensity_CropMarker', 'Intensity_StdIntensity_CropTubulin', 'Intensity_UpperQuartileIntensity_CropDAPI', 'Intensity_UpperQuartileIntensity_CropHPAH', 'Intensity_UpperQuartileIntensity_CropHPAL', 'Intensity_UpperQuartileIntensity_CropMarker', 'Intensity_UpperQuartileIntensity_CropTubulin', 'Location_CenterMassIntensity_X_CropDAPI', 'Location_CenterMassIntensity_X_CropHPAH', 'Location_CenterMassIntensity_X_CropHPAL', 'Location_CenterMassIntensity_X_CropMarker', 'Location_CenterMassIntensity_X_CropTubulin', 'Location_CenterMassIntensity_Y_CropDAPI', 'Location_CenterMassIntensity_Y_CropHPAH', 'Location_CenterMassIntensity_Y_CropHPAL', 'Location_CenterMassIntensity_Y_CropMarker', 'Location_CenterMassIntensity_Y_CropTubulin', 'Location_CenterMassIntensity_Z_CropDAPI', 'Location_CenterMassIntensity_Z_CropHPAH', 'Location_CenterMassIntensity_Z_CropHPAL', 'Location_CenterMassIntensity_Z_CropMarker', 'Location_CenterMassIntensity_Z_CropTubulin', 'Location_Center_X', 'Location_Center_Y', 'Location_MaxIntensity_X_CropDAPI', 'Location_MaxIntensity_X_CropHPAH', 'Location_MaxIntensity_X_CropHPAL', 'Location_MaxIntensity_X_CropMarker', 'Location_MaxIntensity_X_CropTubulin', 'Location_MaxIntensity_Y_CropDAPI', 'Location_MaxIntensity_Y_CropHPAH', 'Location_MaxIntensity_Y_CropHPAL', 'Location_MaxIntensity_Y_CropMarker', 'Location_MaxIntensity_Y_CropTubulin', 'Location_MaxIntensity_Z_CropDAPI', 'Location_MaxIntensity_Z_CropHPAH', 'Location_MaxIntensity_Z_CropHPAL', 'Location_MaxIntensity_Z_CropMarker', 'Location_MaxIntensity_Z_CropTubulin', 'Number_Object_Number', 'Parent_Cells', 'Parent_FilteredNuclei']
    todrop = ['FileName_HPA_H','FileName_HPA_L','FileName_marker','FileName_tubulin','PathName_DAPI','PathName_HPA_H','PathName_HPA_L','PathName_marker','PathName_tubulin','AreaShape_BoundingBoxArea', 'AreaShape_BoundingBoxMaximum_X', 'AreaShape_BoundingBoxMaximum_Y', 'AreaShape_BoundingBoxMinimum_X', 'AreaShape_BoundingBoxMinimum_Y', 'AreaShape_Center_X', 'AreaShape_Center_Y', 'AreaShape_Compactness', 'AreaShape_ConvexArea', 'AreaShape_Eccentricity', 'AreaShape_EquivalentDiameter', 'AreaShape_EulerNumber', 'AreaShape_Extent', 'AreaShape_FormFactor', 'AreaShape_MajorAxisLength', 'AreaShape_MaxFeretDiameter', 'AreaShape_MaximumRadius', 'AreaShape_MeanRadius', 'AreaShape_MedianRadius', 'AreaShape_MinFeretDiameter', 'AreaShape_MinorAxisLength', 'AreaShape_Orientation', 'AreaShape_Perimeter', 'AreaShape_Solidity','Intensity_IntegratedIntensityEdge_CropDAPI', 'Intensity_IntegratedIntensityEdge_CropHPAH', 'Intensity_IntegratedIntensityEdge_CropHPAL', 'Intensity_IntegratedIntensityEdge_CropMarker', 'Intensity_IntegratedIntensityEdge_CropTubulin', 'Intensity_MADIntensity_CropDAPI', 'Intensity_MADIntensity_CropHPAH', 'Intensity_MADIntensity_CropHPAL', 'Intensity_MADIntensity_CropMarker', 'Intensity_MADIntensity_CropTubulin', 'Intensity_MassDisplacement_CropDAPI', 'Intensity_MassDisplacement_CropHPAH', 'Intensity_MassDisplacement_CropHPAL', 'Intensity_MassDisplacement_CropMarker', 'Intensity_MassDisplacement_CropTubulin', 'Intensity_MaxIntensityEdge_CropDAPI', 'Intensity_MaxIntensityEdge_CropHPAH', 'Intensity_MaxIntensityEdge_CropHPAL', 'Intensity_MaxIntensityEdge_CropMarker', 'Intensity_MaxIntensityEdge_CropTubulin','Intensity_MeanIntensityEdge_CropDAPI', 'Intensity_MeanIntensityEdge_CropHPAH', 'Intensity_MeanIntensityEdge_CropHPAL', 'Intensity_MeanIntensityEdge_CropMarker', 'Intensity_MeanIntensityEdge_CropTubulin', 'Intensity_MinIntensityEdge_CropDAPI', 'Intensity_MinIntensityEdge_CropHPAH', 'Intensity_MinIntensityEdge_CropHPAL', 'Intensity_MinIntensityEdge_CropMarker', 'Intensity_MinIntensityEdge_CropTubulin','Intensity_StdIntensityEdge_CropDAPI', 'Intensity_StdIntensityEdge_CropHPAH', 'Intensity_StdIntensityEdge_CropHPAL', 'Intensity_StdIntensityEdge_CropMarker', 'Intensity_StdIntensityEdge_CropTubulin','Location_CenterMassIntensity_X_CropDAPI', 'Location_CenterMassIntensity_X_CropHPAH', 'Location_CenterMassIntensity_X_CropHPAL', 'Location_CenterMassIntensity_X_CropMarker', 'Location_CenterMassIntensity_X_CropTubulin', 'Location_CenterMassIntensity_Y_CropDAPI', 'Location_CenterMassIntensity_Y_CropHPAH', 'Location_CenterMassIntensity_Y_CropHPAL', 'Location_CenterMassIntensity_Y_CropMarker', 'Location_CenterMassIntensity_Y_CropTubulin', 'Location_CenterMassIntensity_Z_CropDAPI', 'Location_CenterMassIntensity_Z_CropHPAH', 'Location_CenterMassIntensity_Z_CropHPAL', 'Location_CenterMassIntensity_Z_CropMarker', 'Location_CenterMassIntensity_Z_CropTubulin', 'Location_Center_X', 'Location_Center_Y', 'Location_MaxIntensity_X_CropDAPI', 'Location_MaxIntensity_X_CropHPAH', 'Location_MaxIntensity_X_CropHPAL', 'Location_MaxIntensity_X_CropMarker', 'Location_MaxIntensity_X_CropTubulin', 'Location_MaxIntensity_Y_CropDAPI', 'Location_MaxIntensity_Y_CropHPAH', 'Location_MaxIntensity_Y_CropHPAL', 'Location_MaxIntensity_Y_CropMarker', 'Location_MaxIntensity_Y_CropTubulin', 'Location_MaxIntensity_Z_CropDAPI', 'Location_MaxIntensity_Z_CropHPAH', 'Location_MaxIntensity_Z_CropHPAL', 'Location_MaxIntensity_Z_CropMarker', 'Location_MaxIntensity_Z_CropTubulin']
    #tokeep = []
    non_measurement = ['ImageNumber','ObjectNumber','FileName_DAPI','Plate','Well','Marker','HPA_Ab','HPA_gene','Fixation_method']
    merged_df_list =[]
    for plate in plate_list:
        df_list =[]
        for dirPath, dirNames, files in os.walk(dir):
            for file in files:
                if not plate in file:
                    continue
                if file.endswith('Image.csv') | file.endswith('Experiment.csv') | file.endswith('_Nuclei.csv'):
                    continue
                filepath = dirPath + os.sep + file
                print(filepath)
                #exp = find_exp(file)
                prefix = find_prefix(file)
                df = pd.read_csv(filepath)
                #df_full = pd.read_csv(filepath)
                #df = df_full[(df_full['ImageNumber'] % 16 == 1) | (df_full['ImageNumber'] % 16 == 2)]
                
                if prefix == 'Cyto_':
                    df = df.drop(columns=['AreaShape_Area'])
                df = df.drop(columns=todrop)
                df = df.add_prefix(prefix)
                df.rename(columns={prefix + 'ImageNumber': 'ImageNumber', prefix + 'ObjectNumber': 'ObjectNumber', prefix + 'FileName_DAPI': 'FileName_DAPI'}, inplace=True)
                df['Plate'] = plate
                df['Well'] = df['FileName_DAPI'].apply(find_well)
                df['Marker'] = df['FileName_DAPI'].apply(find_marker)
                df['HPA_Ab'] = df['FileName_DAPI'].apply(find_hpa_ab)
                df['HPA_gene'] = df['FileName_DAPI'].apply(find_hpa_gene)
                df['Fixation_method'] = df['FileName_DAPI'].apply(find_fixation_method)
        
                df_list.append(df)

        merged_df = reduce(lambda left,right: pd.merge(left,right,on=non_measurement), df_list)
        merged_df['Cyto_AreaShape_Area']=merged_df['Cell_AreaShape_Area'] - merged_df['Nuc_AreaShape_Area']
        merged_df_list.append(merged_df)
    merged_df_all = pd.concat(merged_df_list,ignore_index=True)
        
    return merged_df_all

#%%
def read_image_quality(dir):
    '''
    Read in the ImageQuality output from CellProfiler and merge into one dataframe
    '''
    tokeep = ['FileName_DAPI','ImageNumber', 'ImageQuality_Correlation_DAPI_30',
 'ImageQuality_Correlation_HPA_H_30',
 'ImageQuality_Correlation_HPA_L_30',
 'ImageQuality_Correlation_marker_30',
 'ImageQuality_Correlation_tubulin_30',
 'ImageQuality_FocusScore_DAPI',
 'ImageQuality_FocusScore_HPA_H',
 'ImageQuality_FocusScore_HPA_L',
 'ImageQuality_FocusScore_marker',
 'ImageQuality_FocusScore_tubulin',
 'ImageQuality_LocalFocusScore_DAPI_30',
 'ImageQuality_LocalFocusScore_HPA_H_30',
 'ImageQuality_LocalFocusScore_HPA_L_30',
 'ImageQuality_LocalFocusScore_marker_30',
 'ImageQuality_LocalFocusScore_tubulin_30', 'ImageQuality_MaxIntensity_DAPI',
 'ImageQuality_MaxIntensity_HPA_H',
 'ImageQuality_MaxIntensity_HPA_L',
 'ImageQuality_MaxIntensity_marker',
 'ImageQuality_MaxIntensity_tubulin', 'ImageQuality_PercentMaximal_DAPI',
 'ImageQuality_PercentMaximal_HPA_H',
 'ImageQuality_PercentMaximal_HPA_L',
 'ImageQuality_PercentMaximal_marker',
 'ImageQuality_PercentMaximal_tubulin','ModuleError_01Images',
 'ModuleError_02Metadata',
 'ModuleError_03NamesAndTypes',
 'ModuleError_04Groups',
 'ModuleError_05MeasureImageQuality']
    df_list =[]
    for dirPath, dirNames, files in os.walk(dir):
        for file in files:
            if not file.endswith('Image.csv'):
                 continue
            filepath = dirPath + os.sep + file
            print(filepath)
            df = pd.read_csv(filepath)
            df = df[tokeep]
            df['Plate'] = df['FileName_DAPI'].apply(find_plate)
            df['Well'] = df['FileName_DAPI'].apply(find_well)
            df['Marker'] = df['FileName_DAPI'].apply(find_marker)
            df['HPA_Ab'] = df['FileName_DAPI'].apply(find_hpa_ab)
            #df['HPA_gene'] = df['FileName_DAPI'].apply(find_hpa_gene)
            df['Fixation_method'] = df['FileName_DAPI'].apply(find_fixation_method)
            df_list.append(df)
    df_image_quality = pd.concat(df_list,ignore_index=True)
    return df_image_quality


############## Functions for plotting a summary of different measurements for quality control ####################
#%%
def plot_distribution(df,col1 = 'Plate',col2 = 'Cell_AreaShape_Area',annotate = False):
    df_sub = df[[col1,col2]]
    
    order = list(df[col1].unique())
    plt.figure(figsize=(16,9))
    boxplot = sns.boxplot(y=col2,x=col1,data=df_sub,palette="colorblind",order=order,showfliers=False)
    boxplot = sns.stripplot(y=col2,x=col1,data=df_sub,palette="colorblind",order=order,alpha=0.4)
    boxplot.set(xlabel = None)
    boxplot.set_ylim()
    
    if annotate:
        pairs = list(combinations(order, r=2))
        annotator = Annotator(ax=boxplot,pairs=pairs,data=df_sub,y=col2,x=col1,order=order)
        annotator.configure(test="Mann-Whitney",text_format = 'star')
        #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
        annotator.apply_and_annotate()
    
    plt.title("Distribution: {}".format(col2))
    plt.show()

#%%
def plot_violin(df,col1 = 'Plate',col2 = 'Cell_AreaShape_Area',annotate = False):
    df_sub = df[[col1,col2]]
    
    if col1 == 'Plate':
        order = ["GLS_1","GLS_2","GLS_3","GLS_4","GLS_5",
                "TUFM_1","TUFM_2","TUFM_3","TUFM_4","TUFM_5",
                "IDH3A_1","IDH3A_2","IDH3A_3","IDH3A_4","IDH3A_5",
                "SOD2_1","SOD2_2","SOD2_3","SOD2_5"]
    else:
        order = list(df[col1].unique())
    plt.figure(figsize=(16,9))
    boxplot = sns.violinplot(y=col2,x=col1,data=df_sub,palette="colorblind",order=order,showfliers=False)
    #boxplot = sns.stripplot(y=col2,x=col1,data=df_sub,palette="colorblind",order=order,alpha=0.4)
    boxplot.set(xlabel = None)
    boxplot.set_ylim()
    if annotate:
        pairs = list(combinations(order, r=2))
        annotator = Annotator(ax=boxplot,pairs=pairs,data=df_sub,y=col2,x=col1,order=order)
        annotator.configure(test="Mann-Whitney",text_format = 'star')
        #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
        annotator.apply_and_annotate()
    plt.title("Distribution: {}".format(col2))
    plt.show()
#%%
def plot_hist(df,cond='Plate',val ='IDH3A_1',col='Cyto_AreaShape_Area',add_line = False):
    df_sub = df[df[cond] == val]
    plt.hist(df_sub[col],bins=500)
    #plt.vlines((1, 3, 5,), 0, 10, colors = ("r", "g", "b"), linestyles = ("solid", "dashed", "dotted"))
    #plt.axvline(x = 2.5, ymin = 0.25, ymax = 0.75,linewidth = 4, linestyle ="--",color ='red')
    if add_line:
        mean_val = df_sub[col].mean()
        lower = df_sub[col].mean() - 3*df_sub[col].std()
        upper = df_sub[col].mean() + 3*df_sub[col].std()
        plt.axvline(mean_val,linestyle ='--',color='red')
        plt.axvline(lower,linestyle ='--',color='red')
        plt.axvline(upper,linestyle ='--',color='red')
    plt.xlim()
    plt.xlabel(col)
    plt.title("{} ({}:{})".format(col,cond,val))
    plt.show()
#%%
def plot_kde(df, col1='Cyto_Intensity_MeanIntensity_CropTubulin',col2='Plate'):
    plt.figure(figsize=(12,8))
    kdeplot = sns.kdeplot(data = df, x=col1,hue=col2)
    kdeplot.set(xlabel = col1) 
    kdeplot.set_xlim()
    plt.title('{} by {}'.format(col1,col2))


######################## Functions for filtering objects ##################
#%%
def filter_empty_obj(df):
    df_filtered = df[df['Cyto_AreaShape_Area']!=0]
    print('Before: {}; After: {}'.format(df.shape[0],df_filtered.shape[0]))
    return df_filtered
#%%
def filter_by_nuc_ratio(df,upper=0.7,lower=0.05):
    #mean_val = df['Nuc/Cell'].mean()
    #lower = df['Nuc/Cell'].mean() - 3*df['Nuc/Cell'].std()
    #upper = df['Nuc/Cell'].mean() + 3*df['Nuc/Cell'].std()
    df_filtered = df[(df['Nuc/Cell']<=upper)&(df['Nuc/Cell']>=lower)]

    print('Before: {}; After: {}'.format(df.shape[0],df_filtered.shape[0]))
    return df_filtered
#%%
def filter_by_max_intensity(df):
    df_filtered = df[(df['Cell_Intensity_MaxIntensity_CropDAPI']!=1)&(df['Cell_Intensity_MaxIntensity_CropMarker']!=1)&(df['Cell_Intensity_MaxIntensity_CropTubulin']!=1)]
    print('Before: {}; After: {}'.format(df.shape[0],df_filtered.shape[0]))
    return df_filtered
#%%
def filter_by_tubulin(df,lower=0.025,upper=0.3):
    df_filtered = df[(df['Cyto_Intensity_UpperQuartileIntensity_CropTubulin']<=upper)&(df['Cyto_Intensity_UpperQuartileIntensity_CropTubulin']>=lower)]
    print('Before: {}; After: {}'.format(df.shape[0],df_filtered.shape[0]))
    return df_filtered    

##Filter bad images
############# TBD

#%%
####################### Select measurement for target protein #########################
def HPA_channel_selection(df,cutoff = 0.05):
    dflist=[]
    for plate in list(df["Plate"].unique()):
        
        df_plate = df[df['Plate']==plate]
        for well in list(df_plate["Well"].unique()):
            df_well = df_plate[df_plate['Well']==well]
            total_obj = df_well.shape[0]
            #sat_obj = df_well[df_well['Cyto_Intensity_MaxIntensity_CropHPAH']==1].shape[0]
            sat_obj=df_well[df_well['Cell_Intensity_MaxIntensity_CropHPAH']==1].shape[0]
            ratio = sat_obj/total_obj
            df_well['Ratio'] = ratio
            if ratio >= cutoff:
                df_well['HPA'] = 'HPAL'
            else:
                df_well['HPA'] = 'HPAH'
            dflist.append(df_well)
    df_exp = pd.concat(dflist,ignore_index=True)
    df_exp['Well_id'] = df_exp['FileName_DAPI'].apply(find_well_id)
    df_exp['Plate_well'] = df_exp['Well_id'].apply(find_plate_well)
    print('{},{}'.format(df.shape[0],df_exp.shape[0]))
    return df_exp
#%%
def make_HPA_column(df_well, col_list):
    '''This function reads the selected HPA channel for a well
    and assign the values from the selected column to the new column.'''
    HPA_channel = df_well['HPA'].unique()[0]
    if HPA_channel=="HPAH":
        c = 'H'
    else:
        c = 'L'
    for col in col_list:
        df_well[col] = df_well[col+c]
    return df_well
def update_HPA_columns(df):
    '''This function makes a new column for the HPA channel (ends with "_CropHPA") for each of the listed measurements.
     Then the make_HPA_column function is applied to individual wells and the new columns are filled with values from either "_CropHPAH" or "_CropHPAL". '''
    col_list = []
    measurements = ['MeanIntensity','MaxIntensity','UpperQuartileIntensity','LowerQuartileIntensity','MinIntensity','MedianIntensity','StdIntensity','IntegratedIntensity']
    for measurement in measurements:
        for compartment in ['Nuc','Cyto','Cell']:
            col = '{}_Intensity_{}_CropHPA'.format(compartment,measurement)
            df[col]=0
            col_list.append(col)
    df = df.groupby(['Well_id'],as_index = False).apply(make_HPA_column,col_list)
    return df

#%%
####################### Funcions for statistical calculations #########################
# coefficient of variation
def cv(lst):
    '''This function takes a list of values as an input and returns the coefficient of variance of that list'''
    if len(lst)==0:
        return None
    return stats.variation(lst,nan_policy="omit", ddof=1)

#%%
#Gini index
def gini(lst):
    '''This function takes a list of values as an input and returns the Gini index of that list.
    n is the number of values in the list.
    s is the sum of all the values in the list.
    r is the sorted list of values.
    g is the cumulative sum of the Gini coefficients for each value in the sorted list.
    The loop iterates over each value in the sorted list, adding the Gini coefficient for that value to the running total g.
    The final Gini index is calculated by dividing the sum of the Gini coefficients by the product of n and s.
    Note that this implementation assumes that the input list contains non-negative values. If the input list contains negative values, the implementation needs to be adjusted accordingly.'''
    n = len(lst)
    s = sum(lst)
    r = sorted(lst)
    g = 0
    for i in range(n):
        g += (2 * (i + 1) - n - 1) * r[i]
    return g / (n * s)
#%%
def log_transform(lst):
    new_lst =[]
    for i in lst:
        new_lst.append(np.log10(i))
    return new_lst
def minmax_normalization(lst):
    new_lst =[]
    max_int = max(lst)
    min_int = min(lst)
    for i in lst:
        new_lst.append((i-min_int)/(max_int-min_int))
    return new_lst
#%%
def spearman_correlation(x, y):
    '''This function takes two lists of values as inputs and returns the Spearman correlation between them
    '''
    if len(x) != len(y):
        raise ValueError("Lists must have the same length.")
    return stats.spearmanr(x, y)[0]
def spearman_pval(x, y):
    '''This function takes two lists of values as inputs and returns the Spearman correlation p-value between them
    '''
    if len(x) != len(y):
        raise ValueError("Lists must have the same length.")
    return stats.spearmanr(x, y)[1]

#%%
####################### Funcions for analysis #########################

def generate_df_by_well(df, col = "MeanIntensity"): 
    '''This function makes calculations on the raw measurements, 
    and generate a new dataframe with statistics (CV,Gini,Spearman correlation) for each sample.'''

    #df['Plate_well'] = df['Well_id'].apply(find_plate_well)
    compartments = ['Nuc','Cyto','Cell']
    channels = ['DAPI','Marker','HPA','Tubulin']
    info_cols =['Well_id','Plate_well','Marker','HPA_Ab','HPA_gene','Plate','Fixation_method','HPA']
    columns =[]
    for compartment in compartments:
        for channel in channels:
            columns.append('{}_Intensity_{}_Crop{}'.format(compartment,col,channel))
            #columns.append('{}_Intensity_UpperQuartileIntensity_Crop{}'.format(compartment,channel))
    df_by_well = df.groupby(info_cols)[columns].agg(lambda x: list(x)).reset_index()
    cv_columns=[]
    
    for compartment in compartments:
        for target in channels:
            colname_cv = '{}_CV_{}'.format(compartment,target)
            colname_gini = '{}_Gini_{}'.format(compartment,target)
            cv_columns.append(colname_cv)
            df_by_well[colname_cv] = df_by_well['{}_Intensity_{}_Crop{}'.format(compartment,col,target)].apply(cv)
            df_by_well[colname_gini] = df_by_well['{}_Intensity_{}_Crop{}'.format(compartment,col,target)].apply(gini)
    
    df_by_well['Number_of_objects'] = df_by_well['Nuc_Intensity_{}_CropDAPI'.format(col)].apply(len)
    m = 'Cyto_Intensity_{}_CropMarker'.format(col)
    h = 'Cyto_Intensity_{}_CropHPA'.format(col)
    tub = 'Cyto_Intensity_{}_CropTubulin'.format(col)
    dapi = 'Nuc_Intensity_{}_CropDAPI'.format(col)
    m_log = 'Cyto_Intensity_{}_CropMarker_log'.format(col)
    h_log = 'Cyto_Intensity_{}_CropHPA_log'.format(col)
    m_norm = 'Cyto_Intensity_{}_CropMarker_norm'.format(col)
    h_norm = 'Cyto_Intensity_{}_CropHPA_norm'.format(col)
    df_by_well[m_log] = df_by_well[m].apply(log_transform)
    df_by_well[h_log] = df_by_well[h].apply(log_transform)
    df_by_well[m_norm] = df_by_well[m_log].apply(minmax_normalization)
    df_by_well[h_norm] = df_by_well[h_log].apply(minmax_normalization)

    df_by_well['Spearman_Cyto_MarkervsTub'] = df_by_well.apply(lambda row: spearman_correlation(row[m],row[tub]), axis=1)
    df_by_well['Spearman_Cyto_HPAvsTub'] = df_by_well.apply(lambda row: spearman_correlation(row[h],row[tub]), axis=1)
    
    df_by_well['Spearman_Cyto_MarkervsDAPI'] = df_by_well.apply(lambda row: spearman_correlation(row[m],row[dapi]), axis=1)
    df_by_well['Spearman_Cyto_HPAvsDAPI'] = df_by_well.apply(lambda row: spearman_correlation(row[h],row[dapi]), axis=1)
    
    df_by_well['Spearman_Cyto_HPAvsMarker'] = df_by_well.apply(lambda row: spearman_correlation(row[m],row[h]), axis=1)
    df_by_well['Spearman_Cyto_HPAvsMarker_log10'] = df_by_well.apply(lambda row: spearman_correlation(row[m_log],row[h_log]), axis=1)
    df_by_well['Spearman_Cyto_HPAvsMarker_norm'] = df_by_well.apply(lambda row: spearman_correlation(row[m_norm],row[h_norm]), axis=1)
    df_by_well['Spearman_pval'] = df_by_well.apply(lambda row: spearman_pval(row[m_norm],row[h_norm]), axis=1)
    df_by_well = df_by_well.explode(cv_columns)
    df_by_well = df_by_well.drop(columns=columns)
    return df_by_well


#%%
def plot_correlation(df_raw,well_id,col = "MeanIntensity"):
    '''This functions plots the linear correlation of the mean intensities of Marker and target protein in a sample.
    '''
    df_subset = df_raw[df_raw["Well_id"] == well_id]
    marker = df_subset['Marker'].unique()[0]
    HPA_gene = df_subset['HPA_gene'].unique()[0]

    col_x = 'Cyto_Intensity_{}_CropMarker'.format(col)
    col_y = 'Cyto_Intensity_{}_CropHPA'.format(col)
    col_x_log = 'Cyto_Intensity_{}_CropMarker_log'.format(col)
    col_y_log = 'Cyto_Intensity_{}_CropHPA_log'.format(col)
    col_x_norm = 'Cyto_Intensity_{}_CropMarker_norm'.format(col)
    col_y_norm = 'Cyto_Intensity_{}_CropHPA_norm'.format(col)
    df_subset[col_x_log] = df_subset[col_x].apply(np.log10)
    df_subset[col_y_log] = df_subset[col_y].apply(np.log10)
    df_subset[col_x_norm] = (df_subset[col_x_log]-min(df_subset[col_x_log]))/(max(df_subset[col_x_log])-min(df_subset[col_x_log]))
    df_subset[col_y_norm] = (df_subset[col_y_log]-min(df_subset[col_y_log]))/(max(df_subset[col_y_log])-min(df_subset[col_y_log]))
    #df_subset_norm_filtered = df_subset_norm[df_subset_norm[col_x]]
    #df_subset_norm = df_subset_norm[(df_subset_norm[col_x] <= 0.4) & (df_subset_norm[col_y] <= 0.4)]
    
    spearman_r = round(spearman_correlation(df_subset[col_x_norm],df_subset[col_y_norm]),2)
    spearman_p = '{:0.1e}'.format(spearman_pval(df_subset[col_x_norm],df_subset[col_y_norm]))
    
    plt.figure(figsize=(6,6))
    sns.set_theme()
    corr_plot = sns.regplot(data=df_subset,x=col_x_norm,y=col_y_norm,scatter_kws={'s':0.5},line_kws={'lw':1},label="R={} p={}".format(spearman_r,spearman_p))
    #corr_plot = sns.regplot(data=df_subset,x=col_x,y=col_y,scatter_kws={'s':0.5},line_kws={'lw':1},label="R={} p={}".format(spearman_r,spearman_p))
    #corr_plot = sns.regplot(data=df_subset_norm,x=col_x,y=col_y,hue="Nuc_AreaShape_Area",scatter_kws={'s':0.5},line_kws={'lw':1},label="R={} p={}".format(spearman_r[0],spearman_r[1]))
    corr_plot.legend(loc='upper right',prop={'size': 10})
    #corr_plot.set(xlim = (0,1))
    #corr_plot.set(ylim = (0,1))
    #corr_plot.set(xscale = "log2")
    #corr_plot.set(yscale = "log2")
    corr_plot.set(xlabel = 'Mean Intensity {}'.format(marker))
    corr_plot.set(ylabel = 'Mean Intensity {}'.format(HPA_gene))
    #plt.show()
    filename = "{}_{}.png".format(marker,HPA_gene)
    #plt.savefig(filename)
    return corr_plot

#%%
def plot_cv(df_cv):
    order = list(df_cv['Target'].unique())
    plt.figure(figsize=(9,10))
    cv_plot = sns.boxplot(y='CV',x='Target',data=df_cv,palette="colorblind",order=order,showfliers=False)
    cv_plot = sns.stripplot(y='CV',x='Target',data=df_cv,palette="colorblind",order=order,alpha=0.4,hue='Marker')
    cv_plot.set(xlabel = None)
    
    pairs = list(combinations(order, r=2))
    annotator = Annotator(ax=cv_plot,pairs=pairs,data=df_cv,y='CV',x='Target',order=order)
    annotator.configure(test="Kruskal",text_format = 'simple')
    #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
    annotator.apply_and_annotate()
    plt.title("Coefficient of variance")
    plt.show()

#%%
def plot_gini(df_gini):
    order = list(df_gini['Target'].unique())
    plt.figure(figsize=(9,10))
    gini_plot = sns.boxplot(y='Gini Index',x='Target',data=df_gini,palette="colorblind",order=order,showfliers=False)
    gini_plot = sns.stripplot(y='Gini Index',x='Target',data=df_gini,palette="colorblind",order=order,alpha=0.4)
    gini_plot.set(xlabel = None)
    
    pairs = list(combinations(order, r=2))
    annotator = Annotator(ax=gini_plot,pairs=pairs,data=df_gini,y='Gini Index',x='Target',order=order)
    annotator.configure(test="Kruskal",text_format = 'simple')
    #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
    annotator.apply_and_annotate()
    plt.title("Gini Index of mean intensity")
    plt.show()


#%%
def plot_boxplot(df_analysis,colx='Marker',coly='Spearman_Cyto_HPAvsDAPI'):
    order = list(df_analysis[colx].unique())
    plt.figure(figsize=(9,10))
    sns.set_theme()
    box_plot = sns.boxplot(y=coly,x=colx,data=df_analysis,palette="colorblind",order=order,showfliers=False)
    box_plot = sns.stripplot(y=coly,x=colx,data=df_analysis,palette="colorblind",order=order,alpha=0.4)
    box_plot.set(xlabel = None)
    
    pairs = list(combinations(order, r=2))
    annotator = Annotator(ax=box_plot,pairs=pairs,data=df_analysis,y=coly,x=colx,order=order)
    annotator.configure(test="Kruskal",text_format = 'simple')
    #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
    annotator.apply_and_annotate()
    plt.title(coly)
    plt.show()
#################################### Main ##############################
#%%
output_dir = 'D:\CellProfiler_output\output'
output_imagequality = 'E:\CellProfiler_output\ImageQuality'

#%%
df_merged = read_output(output_dir)

#print(df_merged.columns.tolist())
print(df_merged.shape)
# %%
df_merged.to_csv("df_merged.csv",index=False)


#%%
df = pd.read_csv('D:\CellProfiler_output\df_merged_subset.csv')
#%%
df_imagequality = read_image_quality(output_imagequality)

# %%
df['Nuc/Cell'] = df['Nuc_AreaShape_Area']/df['Cell_AreaShape_Area']
# %%
df_filtered = filter_empty_obj(df)
# %%
df_filtered2 = filter_by_nuc_ratio(df_filtered)
# %%
df_filtered3 = filter_by_max_intensity(df_filtered2)
# %%
df_filtered4 = filter_by_tubulin(df_filtered3)

#%%
df_exp = HPA_channel_selection(df_filtered4)
df_exp = update_HPA_columns(df_exp)

#%%
df_analysis = generate_df_by_well(df_exp)

#%%
thres = df_analysis[df_analysis['Spearman_pval']<0.05]['Spearman_Cyto_HPAvsMarker'].mean()+df_analysis[df_analysis['Spearman_pval']<0.05]['Spearman_Cyto_HPAvsMarker'].std()
df_correlating = df_analysis[(df_analysis['Spearman_pval']<0.05)&(df_analysis['Spearman_Cyto_HPAvsMarker']>thres)]


# %%
df_gini= df_analysis.melt(id_vars=['Well_id','Marker','HPA_gene'], value_vars=['Nuc_Gini_DAPI','Cyto_Gini_Tubulin','Cyto_Gini_HPA','Cyto_Gini_Marker'],var_name='Target',value_name='Gini Index')
plot_gini(df_gini)

# %%
df_heatmap_all = df_correlating.pivot_table(index=['HPA_gene'],columns=['Marker'],values=['Spearman_Cyto_HPAvsMarker'])
df_heatmap_all = df_heatmap_all.fillna(0)
#%%
sns.set(font_scale=1.5)
heatmap_all = sns.clustermap(figsize = (30,35),data=df_heatmap_all,cmap='YlGnBu',xticklabels = True,yticklabels = True)
heatmap_all.set(xlabel=None)
