import cellpose_analysis_utils as utils
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# Read in files
def read_df(directory,filelist):
    '''This function reads in the relevant columns of all the image quality output (.csv) files and combine them into one dataframe.'''
    dfs=[]
    col_list =['FileName_DAPI',
    'ImageQuality_Correlation_DAPI_30', 'ImageQuality_Correlation_HPA_H_30', 'ImageQuality_Correlation_HPA_L_30', 'ImageQuality_Correlation_Marker_30', 'ImageQuality_Correlation_Tubulin_30', 
    'ImageQuality_FocusScore_DAPI', 'ImageQuality_FocusScore_HPA_H', 'ImageQuality_FocusScore_HPA_L', 'ImageQuality_FocusScore_Marker', 'ImageQuality_FocusScore_Tubulin', 
    'ImageQuality_LocalFocusScore_DAPI_30', 'ImageQuality_LocalFocusScore_HPA_H_30', 'ImageQuality_LocalFocusScore_HPA_L_30', 'ImageQuality_LocalFocusScore_Marker_30', 'ImageQuality_LocalFocusScore_Tubulin_30',
    'ImageQuality_PercentMaximal_DAPI', 'ImageQuality_PercentMaximal_HPA_H', 'ImageQuality_PercentMaximal_HPA_L', 'ImageQuality_PercentMaximal_Marker', 'ImageQuality_PercentMaximal_Tubulin', 
    'ImageQuality_PercentMinimal_DAPI', 'ImageQuality_PercentMinimal_HPA_H', 'ImageQuality_PercentMinimal_HPA_L', 'ImageQuality_PercentMinimal_Marker', 'ImageQuality_PercentMinimal_Tubulin', 
    'ImageQuality_PowerLogLogSlope_DAPI', 'ImageQuality_PowerLogLogSlope_HPA_H', 'ImageQuality_PowerLogLogSlope_HPA_L', 'ImageQuality_PowerLogLogSlope_Marker', 'ImageQuality_PowerLogLogSlope_Tubulin', 
    ]
    for file in filelist:
        df = pd.read_csv(os.path.join(directory,file),usecols=col_list)
        dfs.append(df)
        #print(f'{file} read in')
    finaldf= pd.concat(dfs,ignore_index=True)
    del dfs
    return finaldf

def plot_hist(df,m,plate=None):
    '''Plot the distribution of the selected metrics.
    Parameters:
    df: the dataframe
    m: the selected metrics
    plate:'''
    arr = df[m]
    plt.hist(arr,bins=50)
    plt.axvline(np.mean(arr),linestyle='--',color='red')
    plt.axvline(np.percentile(arr,1),linestyle='--',color='blue')
    plt.axvline(np.percentile(arr,25),linestyle='--',color='blue')
    plt.axvline(np.percentile(arr,75),linestyle='--',color='blue')
    plt.axvline(np.percentile(arr,99),linestyle='--',color='blue')
    plt.xlabel(m)
    plt.savefig(f'hist_{plate}_{m}.png')
    plt.close()
    print(f'Max: {np.max(arr)}')
    print(f'Min: {np.min(arr)}')
    print(f'99%: {np.percentile(arr,99)}')
    print(f'1%: {np.percentile(arr,1)}')

# Mark each image
def mark_correlation(row,thres1 =0.2,thres2=0.4):
    '''count: 0-5'''
    count = 0
    for target in ['DAPI','Marker']:
        if row[f'ImageQuality_Correlation_{target}_30'] >= thres1:
            count += 1
    for target in ['Tubulin','HPA_H','HPA_L']:
        if row[f'ImageQuality_Correlation_{target}_30'] >= thres2:
            count += 1
    return count 

def mark_focus_score(row,thres1=0.02,thres2=0.005):
    '''count: 0-2'''
    count = 0
    if row[f'ImageQuality_FocusScore_DAPI'] <= thres1:
        count += 1
    if row[f'ImageQuality_FocusScore_Tubulin'] <= thres2:
        count += 1
    return count

def mark_slope(row,thres = -1):
    '''count: 0-3'''
    count =0
    for target in ['DAPI','Tubulin','Marker']:
        if row[f'ImageQuality_PowerLogLogSlope_{target}'] >= thres:
            count += 1
    return count 

def mark_all(row,thres=3):
    '''count the score for the image (0-10), and mark the image with the score above threshold (thres).
    Returns: 
    True: image with bad quality
    False: normal image'''
    count =0
    for i in ['Correlation','FocusScore','Slope']:
        count += row[i]
    return count >= thres

def find_filename(filename):
    return '_'.join(filename.split('_')[:-1])
def find_wellid(filename):
    return '_'.join(filename.split('_')[:-3])
def find_plate(filename):
    return '_'.join(filename.split('_')[:-4])
def mark_image(df,thres=3):
    '''mark the dataframe with: filename, well id, plate, score for each relevant metrics
    The image is marked to be removed if the total score if above the given threshold'''
    df['Filename'] = df['FileName_DAPI'].apply(find_filename)
    df['Well_id'] = df['FileName_DAPI'].apply(find_wellid)
    df['Plate'] = df['FileName_DAPI'].apply(find_plate)
    df['Correlation'] = df.apply(lambda row:mark_correlation(row),axis=1)
    df['FocusScore'] = df.apply(lambda row:mark_focus_score(row,thres2=0.01),axis=1)
    df['Slope'] = df.apply(lambda row:mark_slope(row),axis=1)
    df['Remove'] = df.apply(lambda row:mark_all(row,thres=thres),axis=1)
    return df


if __name__ == '__main__':
    directory = "/home/xiuqi/ImageQuality"
    dfs=[]
    for marker in ['IDH3A','SOD2','GLS','TUFM']:
        files_list = [f'{marker}_1_Image.csv',f'{marker}_2_Image.csv',f'{marker}_3_Image.csv',f'{marker}_4_Image.csv',f'{marker}_5_Image.csv']
        df = read_df(directory,files_list)
        df = mark_image(df,thres=1)
        '''
        for target in ['DAPI','Tubulin','Marker','HPA_H','HPA_L']:
            plot_hist(df,f'ImageQuality_PowerLogLogSlope_{target}')
            plot_hist(df,f'ImageQuality_FocusScore_{target}')
            plot_hist(df,f'ImageQuality_Correlation_{target}_30')
            plot_hist(df,f'ImageQuality_LocalFocusScore_{target}_30')
        for plate in np.unique(df['Plate']):
            for target in ['DAPI','Tubulin','Marker','HPA_H','HPA_L']:
                plot_hist(df[df['Plate']==plate],f'ImageQuality_PowerLogLogSlope_{target}',plate=plate)
                plot_hist(df[df['Plate']==plate],f'ImageQuality_FocusScore_{target}',plate=plate)
                plot_hist(df[df['Plate']==plate],f'ImageQuality_Correlation_{target}_30',plate=plate)
                plot_hist(df[df['Plate']==plate],f'ImageQuality_LocalFocusScore_{target}_30',plate=plate)
        '''
        print(marker)
        print(df['Correlation'].value_counts())
        print(df['FocusScore'].value_counts())
        print(df['Slope'].value_counts())
        print(df['Remove'].value_counts())
        #df[df['Remove']==True].to_csv(f'{marker}_removed.csv')
        dfs.append(df[df['Remove']==True])
    remove_df = pd.concat(dfs,ignore_index=True)
    del dfs
    remove_df.to_csv('imagetoremove_thres1.csv')