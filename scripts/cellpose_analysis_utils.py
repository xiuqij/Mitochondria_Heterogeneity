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
import time

#%% Processing filename
#SOD2_1_A01_HPA000286_f01
def find_plate(filename):
    filenamesplit = filename.split("_")
    plate = '{}_{}'.format(filenamesplit[0],filenamesplit[1])
    return plate
def find_well(filename):
    filenamesplit = filename.split("_")
    well = filenamesplit[2]
    return well 

def find_marker(filename):
    #IDH3A_5_A01_HPA021995_f01
    filenamesplit = filename.split("_")
    #plate = '{}_{}'.format(filenamesplit[0],filenamesplit[1])
    #well = filenamesplit[2]
    marker = filenamesplit[0]
    return marker
def find_hpa_ab(filename):
    filenamesplit = filename.split("_")
    hpa_ab = filenamesplit[3]
    return hpa_ab
def find_hpa_gene(filename,map="/home/xiuqi/plate_layout_correct.csv"):
    pl_map = pd.read_csv(map)
    dict_HPA_genes = dict(zip(pl_map.Ab,pl_map.gene))
    hpa_ab = find_hpa_ab(filename)
    hpa_gene = dict_HPA_genes.get(hpa_ab)
    return hpa_gene
def find_fixation_method(filename):
    marker = find_marker(filename)
    if marker == 'IDH3A' or marker == 'SOD2':
        fixation = 'Methanol'
    if marker == 'GLS' or marker == 'TUFM':
        fixation = 'PFA'  
    return fixation
def find_plate_well(filename):
    filenamesplit = filename.split("_")
    platewell = '{}_{}'.format(filenamesplit[1],filenamesplit[2])
    return platewell

def find_well_id(filename):
    filenamesplit = filename.split("_")
    well_id = '{}_{}_{}'.format(filenamesplit[0],filenamesplit[1],filenamesplit[2])
    return well_id



# General plotting of distributions
def plot_distribution(df,col_x = "", col_y = "",annotate=False,savename='temp'):

    order = list(df[col_x].unique())
    plt.figure(figsize=(16,9))
    boxplot = sns.boxplot(y=col_y,x=col_x,data=df,palette="colorblind",order=order,showfliers=False)
    boxplot = sns.stripplot(y=col_y,x=col_x,data=df,palette="colorblind",hue='Marker',order=order,alpha=0.4)
    boxplot.set(xlabel = None)
    boxplot.set_ylim()
    
    if annotate:
        pairs = list(combinations(order, r=2))
        annotator = Annotator(ax=boxplot,pairs=pairs,data=df,y=col_y,x=col_x,order=order)
        annotator.configure(test="Mann-Whitney",text_format = 'star')
        #test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, Brunner-Munzel.
        annotator.apply_and_annotate()
    
    plt.title("Distribution: {}".format(col_y))
    plt.savefig(f"{savename}.png")
    plt.close()  # Close the plot window

def plot_kde(df, col1='Cyto_MeanInt_Marker',col2='Well_id'):
    plt.figure(figsize=(12,8))
    pd.options.mode.use_inf_as_na = False

    kdeplot = sns.kdeplot(data = df, x=col1,hue=col2)
    kdeplot.set(xlabel = col1) 
    kdeplot.set_xlim()
    plt.title('{} by {}'.format(col1,col2))
    plt.savefig('kde.png')
    plt.close()
def plot_histogram(df,col='Nuc_Size',savename='histogram',bins=100):
    plt.hist(df[col],bins=bins)
    plt.savefig(f'{savename}.png')
    plt.close()
def plot_histplot(df, col1='Cyto_MeanInt_Marker',col2='Well_id'):
    sns.histplot(df, x=col1, y=col2)
    plt.savefig('histplottufm.png')
    plt.close()
def plot_correlation(df,well_id,col='Mito_MeanInt'):
    df_sample = df[df['Well_id']==well_id]
    marker = df_sample['Marker'].unique()[0]
    HPA_gene = df_sample['HPA_gene'].unique()[0]
    col_x= f'{col}_Marker'
    col_y= f'{col}_HPA'

    #drop outliers
    mean1 = df['Mito_MeanInt_Marker'].mean()
    std1 = df['Mito_MeanInt_Marker'].std()

    mean2 = df['Mito_MeanInt_HPA'].mean()
    std2 = df['Mito_MeanInt_HPA'].std()

    df_sample=df_sample[(df_sample['Mito_MeanInt_Marker'] >= mean1-2*std1)&(df['Mito_MeanInt_Marker']<=mean1+2*std1)]
    df_sample=df_sample[(df_sample['Mito_MeanInt_HPA'] >= mean2-2*std2)&(df['Mito_MeanInt_HPA']<=mean2+2*std2)]
    #df_sample.drop(df_sample[(df_sample['Mito_MeanInt_Marker'] > 0.6)|(df_sample['Mito_MeanInt_HPA'] > 0.6)].index,inplace=True)
 
    spearman_r = round(spearman_correlation(df_sample[col_x],df_sample[col_y]),2)
    spearman_p = '{:0.1e}'.format(spearman_pval(df_sample[col_x],df_sample[col_y]))
    
    plt.figure(figsize=(6,6))
    sns.set_theme()
    corr_plot = sns.regplot(data=df_sample,x=col_x,y=col_y,scatter_kws={'s':0.5},line_kws={'lw':1},label="R={} p={}".format(spearman_r,spearman_p))
    corr_plot.legend(loc='upper right',prop={'size': 10})
    #corr_plot.set(xlim = (0,1))
    #corr_plot.set(ylim = (0,1))
    #corr_plot.set(xscale = "log2")
    #corr_plot.set(yscale = "log2")
    corr_plot.set(xlabel = f'{col}_{marker}')
    corr_plot.set(ylabel = f'{col}_{HPA_gene}')
    filename = "{}_{}.png".format(marker,HPA_gene)
    plt.savefig(filename)
    plt.close()
# Initial QC for obvious segmentation artefacts
# remove cells with more than one nuleus; remove cells with nuclei not completely in the mask
def drop_duplicate_cellIDs(df,savepath=None,name=None,savedrop=True,col_list = ["Filename","Cell_ID","Parent_Cell","Parent_Nuc","Nuc_Size","Cell_Size","Cyto_Size","Cyto_Size_mask"]):
    '''find the objects that have more than one nucleus in a cell and remove them from the raw dataframe.
    a copy of the removed objects is also saved.'''
    duplicated_rows = df[df.duplicated(subset = ['Filename','Cell_ID'],keep=False)]
    #print(f"number of rows with duplicated cell IDs: {duplicated_rows.shape[0]}")
    df.drop_duplicates(subset = ['Filename','Cell_ID'],keep=False,inplace=True)
    if savedrop:
        duplicated_rows[col_list].to_csv(os.path.join(savepath,f"{name}_duplicated_CellID.csv"),sep='\t')
    return df
def drop_out_nucleus(df,savepath=None,name=None,savedrop=True,col_list = ["Filename","Cell_ID","Parent_Cell","Parent_Nuc","Nuc_Size","Cell_Size","Cyto_Size","Cyto_Size_mask"]):
    '''find in the dataframe the objects if the nucleus is not completely within the cell, and remove them from the dataframe.
    a copy of the removed objects is also saved.
    Note that drop_duplicated_cellIDs() should be performed first.'''
    rows_to_remove = df[df['Cyto_Size_mask']-df['Cyto_Size'] > 10]
    #print(f"number of rows to remove: {rows_to_remove.shape[0]}")
    df.drop(rows_to_remove.index, inplace=True)
    if savedrop:
        rows_to_remove[col_list].to_csv(os.path.join(savepath,f"{name}_nucleus_out_of_cell.csv"),sep='\t')
    return df
# Plotting and other QC
# plot nuc ratio
def add_nuc_ratio(df):
    df['Nuc/Cell'] = df['Nuc_Size'] / df['Cell_Size']
    return df
def plot_nuc_ratio(df):
    df = add_nuc_ratio(df)
    plot_histogram(df,col='Nuc/Cell',savename='nuc_ratio_distribution')
# filter by nuc/cell
def filter_by_nuc_ratio(df, upper=0.8, lower=0.05):
    add_nuc_ratio(df)
    df.drop(df[(df['Nuc/Cell'] > upper) | (df['Nuc/Cell'] < lower)].index, inplace=True)
    return df

# filter by max intensity
def add_max_int(df,cols = ['Cell_Intensities_DAPI','Cell_Intensities_Marker','Cell_Intensities_Tubulin']):
    for col in cols:
        df[f'{col}_Max'] = df[col].apply(np.max)
    return df
def filter_by_max_intensity(df,cols = ['Cell_Intensities_DAPI','Cell_Intensities_Marker','Cell_Intensities_Tubulin']):
    add_max_int(df,cols)
    for col in cols:
        max_int_filter = (df[f'{col}_Max'] != 1)
        df.drop(df[~max_int_filter].index,inplace=True)
    df.drop(columns=cols, axis=1, inplace=True)
    return df

# Select HPA channel
def calculate_sat_ratio(df):
    '''Calculate within a well the ratio of objects containing saturated pixel(s).'''
    return (df['Cell_MaxIntensity_HPA_H'] == 1).sum() / df.shape[0]
def select_for_well(well_df,cutoff=0.05):
    '''Create a new column for measurements of the HPA channel .
    If the ratio of objects with saturated pixels in HPA_H channel is higher than the cutoff, the measurements in the HPA_L channel will be used.
    '''
    ratio = calculate_sat_ratio(well_df)
    if ratio >= cutoff:
        well_df['Nuc_Intensities_HPA'] = well_df['Nuc_Intensities_HPA_L']
        well_df['Cell_Intensities_HPA'] = well_df['Cell_Intensities_HPA_L']
        well_df['Cyto_Intensities_HPA'] = well_df['Cyto_Intensities_HPA_L']
    else:
        well_df['Nuc_Intensities_HPA'] = well_df['Nuc_Intensities_HPA_H']
        well_df['Cell_Intensities_HPA'] = well_df['Cell_Intensities_HPA_H']
        well_df['Cyto_Intensities_HPA'] = well_df['Cyto_Intensities_HPA_H']
    return well_df
def HPA_column_selection(df,cutoff = 0.05):
    '''apply the selection well by well, and remove the original columns.'''
    df['Well_id'] = df['Filename'].apply(find_well_id)
    df['Cell_MaxIntensity_HPA_H'] = df['Cell_Intensities_HPA_H'].apply(np.max)
    df = df.groupby('Well_id').apply(lambda x: select_for_well(x, cutoff)).reset_index(drop=True)
    todrop = ['Nuc_Intensities_HPA_H','Cell_Intensities_HPA_H','Cyto_Intensities_HPA_H',
    'Nuc_Intensities_HPA_L','Cell_Intensities_HPA_L','Cyto_Intensities_HPA_L','Cell_MaxIntensity_HPA_H']
    df.drop(columns=todrop, inplace=True)
    return df

# Secondary QC
# plot tubulin intensities distribution
def add_upper_quartile_tub(df):
    df['Cyto_UpperQuartileIntensity_Tubulin'] = df['Cyto_Intensities_Tubulin'].apply(lambda x:np.percentile(x,75))
    #df.drop(columns=['Cyto_Intensities_Tubulin'], inplace=True)
    return df 
def plot_tubulin(df):
    df = add_upper_quartile_tub(df)
    plot_histogram(df,col='Cyto_UpperQuartileIntensity_Tubulin',savename='tubulin_distribution')
# filter by tubulin intensity
def filter_by_tubulin(df,upper,lower):
    '''Filter out objects with lower and upper 1% tubulin intensity(Cyto_UpperQuartile)'''
    add_upper_quartile_tub(df)
    tubulin_filter = (df['Cyto_UpperQuartileIntensity_Tubulin'] < upper) & (df['Cyto_UpperQuartileIntensity_Tubulin'] > lower)
    df.drop(df[~tubulin_filter].index, inplace=True)
    df.drop(columns=['Cyto_UpperQuartileIntensity_Tubulin'], axis=1, inplace=True)
    return df
# filter images
import cellpose_analysis_image_quality as imagequality
def filter_by_imagequality(df,iq_filepath='/home/xiuqi/image_quality/imagetoremove_1.csv'):
    df_iq = pd.read_csv(iq_filepath)
    toremove = df_iq[df_iq['Remove']==True]['Filename'].tolist()
    df = df[~df['Filename'].isin(toremove)]
    return df

# Add infomation for subsetting: marker, HPA, gene, fixation
def add_exp_info(df):
    df['Plate'] = df['Filename'].apply(find_plate)
    #df['Well'] = df['Filename'].apply(find_well)
    df['Marker'] = df['Filename'].apply(find_marker)
    df['HPA_Ab'] = df['Filename'].apply(find_hpa_ab)
    df['HPA_gene'] = df['Filename'].apply(find_hpa_gene)
    return df

# Statistics
# distribution and variation
def cv(arr):
    '''This function takes an array of values as an input and returns the coefficient of variance of that array'''
    if math.isnan(stats.variation(arr,nan_policy="omit", ddof=1)):
        return -1
    return stats.variation(arr,nan_policy="omit", ddof=1)
def gini(arr):
    '''This function takes an array of values as an input and returns the Gini index of that array.
    n is the number of values in the array.
    s is the sum of all the values in the array.
    r is the sorted array of values.
    g is the cumulative sum of the Gini coefficients for each value in the sorted array.
    The loop iterates over each value in the sorted array, adding the Gini coefficient for that value to the running total g.
    The final Gini index is calculated by dividing the sum of the Gini coefficients by the product of n and s.
    Note that this implementation assumes that the input array contains non-negative values. If the input array contains negative values, the implementation needs to be adjusted accordingly.'''
    if len(arr) <=1:
        return -1
    n = len(arr)
    s = np.sum(arr)
    r = np.sort(arr)
    g = 0
    for i in range(n):
        g += (2 * (i + 1) - n - 1) * r[i]
    return g / (n * s)

def calculate_skewness(arr):
    if math.isnan(stats.skew(arr,nan_policy='omit')):
        return -100
    return stats.skew(arr,nan_policy='omit')
def calculate_kurtosis(arr):
    if math.isnan(stats.kurtosis(arr,nan_policy='omit')):
        return -100
    return stats.kurtosis(arr,nan_policy='omit')

# Normalization
def log_transform(arr):
    return np.log10(arr)
def minmax_normalization(arr):
    max_int = np.max(arr)
    min_int = np.min(arr)
    arr_normalized = (arr - min_int)/(max_int - min_int)
    return arr_normalized

#Correlation
def spearman_correlation(arr_x, arr_y):
    '''This function takes two arrays of values as inputs and returns the Spearman correlation between them
    '''
    if len(arr_x) != len(arr_y):
        raise ValueError("Lists must have the same length.")
    if math.isnan(stats.spearmanr(arr_x, arr_y,nan_policy='omit')[0]):
        return -100
    return stats.spearmanr(arr_x, arr_y,nan_policy='omit')[0]
def spearman_pval(arr_x, arr_y):
    '''This function takes two lists of values as inputs and returns the Spearman correlation p-value between them
    '''
    if len(arr_x) != len(arr_y):
        raise ValueError("Lists must have the same length.")
    if math.isnan(stats.spearmanr(arr_x, arr_y,nan_policy='omit')[1]):
        return 100
    return stats.spearmanr(arr_x, arr_y,nan_policy='omit')[1]

def kendalls_tau(arr_x, arr_y):
    '''This function takes two arrays of values as inputs and returns the Spearman correlation between them
    '''
    if len(arr_x) != len(arr_y):
        raise ValueError("Lists must have the same length.")
    return stats.kendalltau(arr_x, arr_y,nan_policy='omit')[0]
def kendalls_tau_pval(arr_x, arr_y):
    '''This function takes two lists of values as inputs and returns the Spearman correlation p-value between them
    '''
    if len(arr_x) != len(arr_y):
        raise ValueError("Lists must have the same length.")
    return stats.kendalltau(arr_x, arr_y,nan_policy='omit')[1]

def pearson_correlation(arr_x,arr_y):
    return stats.pearsonr(arr_x,arr_y)[0]
def pearson_pval(arr_x,arr_y):
    return stats.pearsonr(arr_x,arr_y)[1]



# Add secondary measurements 
def calculate_top_quartile_mean(pixel_values):
    return np.mean(pixel_values[(pixel_values>= np.percentile(pixel_values,75)) & (pixel_values <= np.percentile(pixel_values,99))])
def calculate_mito_mean(marker,hpa):
    mito_mask = ((marker>=np.percentile(marker,75)) & (marker<=np.percentile(marker,99)))
    mito_intensities_marker = marker[mito_mask]
    mito_intensities_hpa = hpa[mito_mask]
    return np.mean(mito_intensities_marker),np.mean(mito_intensities_hpa)

def apply_measurement(df,m='Cyto_MeanInt'):
    '''calculate values for the selected measurement within a cell object.
    Cyto_MeanInt: mean intensity within the compartment.
    Cyto_Q3: upper quartile (75%highest datapoint) intensity within the compartment/
    Cyto_TopQuartileMean: mean intensity of the 75-99% highest pixels within the compartment.
    Mito_MeanInt: mean intensity within the 'mitochondria mask' (75-99% highest pixels in the mitochondria channel).'''
    if m=='Cyto_MeanInt':
        df['Nuc_MeanInt_DAPI']=df['Nuc_Intensities_DAPI'].apply(np.mean)
        df['Cyto_MeanInt_Marker']=df['Cyto_Intensities_Marker'].apply(np.mean)
        df['Cyto_MeanInt_HPA']=df['Cyto_Intensities_HPA'].apply(np.mean)
        df['Cyto_MeanInt_Tubulin']=df['Cyto_Intensities_Tubulin'].apply(np.mean)
    if m=='Cyto_Q3':
        df['Nuc_Q3_DAPI']=df['Nuc_Intensities_DAPI'].apply(lambda x:np.percentile(x,75))
        df['Cyto_Q3_Marker']=df['Cyto_Intensities_Marker'].apply(lambda x:np.percentile(x,75))
        df['Cyto_Q3_HPA']=df['Cyto_Intensities_HPA'].apply(lambda x:np.percentile(x,75))
        df['Cyto_Q3_Tubulin']=df['Cyto_Intensities_Tubulin'].apply(lambda x:np.percentile(x,75))
    if m=='Cyto_TopQuartileMean':
        df['Nuc_TopQuartileMean_DAPI']=df['Nuc_Intensities_DAPI'].apply(calculate_top_quartile_mean)
        df['Cyto_TopQuartileMean_Marker']=df['Cyto_Intensities_Marker'].apply(calculate_top_quartile_mean)
        df['Cyto_TopQuartileMean_HPA']=df['Cyto_Intensities_HPA'].apply(calculate_top_quartile_mean)
        df['Cyto_TopQuartileMean_Tubulin']=df['Cyto_Intensities_Tubulin'].apply(calculate_top_quartile_mean)
    if m=='Mito_MeanInt':
        df['Nuc_MeanInt_DAPI']=df['Nuc_Intensities_DAPI'].apply(np.mean)
        df['Mito_MeanInt_Marker']=df.apply(lambda row:calculate_mito_mean(row['Cyto_Intensities_Marker'],row['Cyto_Intensities_HPA'])[0],axis=1)
        df['Mito_MeanInt_HPA']=df.apply(lambda row:calculate_mito_mean(row['Cyto_Intensities_Marker'],row['Cyto_Intensities_HPA'])[1],axis=1)
        df['Mito_MeanInt_Tubulin']=df.apply(lambda row:calculate_mito_mean(row['Cyto_Intensities_Marker'],row['Cyto_Intensities_Tubulin'])[1],axis=1)
    #df['Cyto_UpperQuartileIntensity_Tubulin'] = df['Cyto_Intensities_Tubulin'].apply(lambda x:np.percentile(x,75))
    return df

def group_by_well(df, m='Cyto_MeanInt'):
    info_cols = ['Well_id', 'Marker', 'HPA_Ab', 'HPA_gene', 'Plate']

    # Define dictionaries to map measurement types to their corresponding columns
    measurement_columns = {
        'Cyto_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker':'Cyto_MeanInt_Marker', 'HPA':'Cyto_MeanInt_HPA', 'Tubulin':'Cyto_MeanInt_Tubulin'
        },
        'Cyto_Q3': {
            'DAPI':'Nuc_Q3_DAPI','Marker':'Cyto_Q3_Marker', 'HPA':'Cyto_Q3_HPA', 'Tubulin':'Cyto_Q3_Tubulin'
        }, 
        'Cyto_TopQuartileMean': {
            'DAPI':'Nuc_TopQuartileMean_DAPI','Marker': 'Cyto_TopQuartileMean_Marker', 'HPA':'Cyto_TopQuartileMean_HPA', 'Tubulin':'Cyto_TopQuartileMean_Tubulin'
        },
        'Mito_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker': 'Mito_MeanInt_Marker', 'HPA':'Mito_MeanInt_HPA', 'Tubulin':'Mito_MeanInt_Tubulin'
        }
    }
    # Get the measurement columns based on the selected measurement type
    mcols = measurement_columns.get(m)
    if not mcols:
        raise ValueError("Invalid measurement type")

    # Group by info_cols and aggregate the selected measurement columns as lists
    df_by_well = df.groupby(info_cols)[list(mcols.values())].agg(list).reset_index()
    # release some memory
    del df
    # Convert the lists in mcols to arrays
    df_by_well[list(mcols.values())].apply(np.array)

    return df_by_well

def get_statistics_by_well(df_sample, m='Cyto_MeanInt'):
    '''calculate statistics on each sample'''
    # Define dictionaries to map measurement types to their corresponding columns
    measurement_columns = {
        'Cyto_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker':'Cyto_MeanInt_Marker', 'HPA':'Cyto_MeanInt_HPA', 'Tubulin':'Cyto_MeanInt_Tubulin'
        },
        'Cyto_Q3': {
            'DAPI':'Nuc_Q3_DAPI','Marker':'Cyto_Q3_Marker', 'HPA':'Cyto_Q3_HPA', 'Tubulin':'Cyto_Q3_Tubulin'
        }, 
        'Cyto_TopQuartileMean': {
            'DAPI':'Nuc_TopQuartileMean_DAPI','Marker': 'Cyto_TopQuartileMean_Marker', 'HPA':'Cyto_TopQuartileMean_HPA', 'Tubulin':'Cyto_TopQuartileMean_Tubulin'
        },
        'Mito_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker': 'Mito_MeanInt_Marker', 'HPA':'Mito_MeanInt_HPA', 'Tubulin':'Mito_MeanInt_Tubulin'
        }
    }
    # Get the measurement columns based on the selected measurement type
    mcols = measurement_columns.get(m)
    if not mcols:
        raise ValueError("Invalid measurement type")

    # Calculate Spearman correlations, kurtosis, skewness, CV, and Gini for each measurement column
    df_sample['Total_objects'] = df_sample[mcols['DAPI']].apply(len)
    df_sample['Spearman_cor'] = df_sample.apply(lambda row: spearman_correlation(row[mcols['Marker']], row[mcols['HPA']]), axis=1)
    df_sample['Spearman_pval'] = df_sample.apply(lambda row: spearman_pval(row[mcols['Marker']], row[mcols['HPA']]), axis=1)
    for chan,col in mcols.items():
        df_sample[f'Kurtosis_{chan}'] = df_sample[col].apply(calculate_kurtosis)
        df_sample[f'Skewness_{chan}'] = df_sample[col].apply(calculate_skewness)
        df_sample[f'CV_{chan}'] = df_sample[col].apply(cv)
        df_sample[f'Gini_{chan}'] = df_sample[col].apply(gini)

    return df_sample

# File I/O
def get_file_paths(directory):
    feather_files=[]
    for root,dirs,files in os.walk(directory):
        for filename in files:
            if filename.endswith('.feather'):
                filepath=os.path.join(root,filename)
                feather_files.append(filepath)
    return feather_files
def get_file_paths_marker(directory,marker='TUFM'):
    feather_files=[]
    for root,dirs,files in os.walk(directory):
        for filename in files:
            if filename.endswith('.feather') and marker in filename:
                filepath=os.path.join(root,filename)
                feather_files.append(filepath)
    return feather_files
def create_measurement_df(filepath,measurement='Cyto_MeanInt'):
    col_list=['Filename','Cell_ID','Cyto_Size', 'Cyto_Size_mask', 'Nuc_Intensities_DAPI',
       'Nuc_Intensities_HPA_H', 'Nuc_Intensities_HPA_L','Cell_Intensities_Marker',
       'Cell_Intensities_Tubulin', 'Cell_Intensities_HPA_H',
       'Cell_Intensities_HPA_L', 'Cell_Intensities_DAPI',
       'Cyto_Intensities_Marker', 'Cyto_Intensities_Tubulin',
       'Cyto_Intensities_HPA_H', 'Cyto_Intensities_HPA_L']
    #Filter
    df_raw = pd.read_feather(filepath,columns=col_list)
    df = drop_duplicate_cellIDs(df_raw,savedrop=False)
    df = drop_out_nucleus(df,savedrop=False)
    df = filter_by_max_intensity(df,cols = ['Cell_Intensities_DAPI','Cell_Intensities_Marker','Cell_Intensities_Tubulin'])
    
    #Add HPA columns and filter by max intensity
    df = HPA_column_selection(df)
    df = filter_by_max_intensity(df,cols=['Cell_Intensities_HPA'])
    # Add experiment info
    df = add_exp_info(df)
    #Filter by tubulin
    marker = df['Marker'].values[0]
    if marker == 'IDH3A':
        df = filter_by_tubulin(df,upper=0.118,lower=0.024)
    if marker == 'GLS':
        df = filter_by_tubulin(df,upper=0.109,lower=0.015)
    if marker == 'TUFM':
        df = filter_by_tubulin(df,upper=0.073,lower=0.010)
    if marker == 'SOD2':
        df = filter_by_tubulin(df,upper=0.135,lower=0.010)
    
    # filter image quality 
    df = filter_by_imagequality(df)
    # get the pixels and measurement
    df = apply_measurement(df,m=measurement)
    #Cyto_MeanInt, Cyto_Q3, Cyto_TopQuartileMean, Mito_MeanInt
    # df by well
    df_by_well= group_by_well(df,m=measurement)
    df_by_well= get_statistics_by_well(df_by_well,m=measurement)
    return df_by_well
'''TUFM
max 0.3625200274662394
min 0.004638742656595712
99% 0.07289234760051881
1% 0.010376134889753566
IDH3A
0.3256618600747692
0.002822919050888838
0.11830319676508735
0.024406805523765927
GLS
0.2751049057755398
0.0025787747005416952
0.10877229724574658
0.014731136034180211
SOD2
0.47226672770275424
0.002319371328297856
0.13521782253757533
0.00988784618905928'''
def create_final_df(directory,measurement='Cyto_MeanInt'):
    dfs=[]
    for path in get_file_paths(directory):
        start= time.time()
        print(f'Processing: {path}')
        dfs.append(create_measurement_df(path,measurement=measurement))
        end = time.time()
        print(f'Processing completed in: {end-start:.2f} seconds.')


    final_df = pd.concat(dfs,ignore_index = True)
    print("Final dataframe generated")
    del dfs
    return final_df


# Functions for handling the analysis dataframe (df by sample)
def correct_duplicates(df, m='Cyto_MeanInt'):
    info_cols = ['Well_id', 'Marker', 'HPA_Ab', 'HPA_gene', 'Plate']

    # Define dictionaries to map measurement types to their corresponding columns
    measurement_columns = {
        'Cyto_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker':'Cyto_MeanInt_Marker', 'HPA':'Cyto_MeanInt_HPA', 'Tubulin':'Cyto_MeanInt_Tubulin'
        },
        'Cyto_Q3': {
            'DAPI':'Nuc_Q3_DAPI','Marker':'Cyto_Q3_Marker', 'HPA':'Cyto_Q3_HPA', 'Tubulin':'Cyto_Q3_Tubulin'
        }, 
        'Cyto_TopQuartileMean': {
            'DAPI':'Nuc_TopQuartileMean_DAPI','Marker': 'Cyto_TopQuartileMean_Marker', 'HPA':'Cyto_TopQuartileMean_HPA', 'Tubulin':'Cyto_TopQuartileMean_Tubulin'
        },
        'Mito_MeanInt': {
            'DAPI':'Nuc_MeanInt_DAPI','Marker': 'Mito_MeanInt_Marker', 'HPA':'Mito_MeanInt_HPA', 'Tubulin':'Mito_MeanInt_Tubulin'
        }
    }
    # Get the measurement columns based on the selected measurement type
    mcols = measurement_columns.get(m)
    if not mcols:
        raise ValueError("Invalid measurement type")

    # Group by info_cols and aggregate the selected measurement columns as lists
    df_by_well = df.groupby(info_cols)[list(mcols.values())].agg(np.hstack).reset_index()
    # release some memory


    # Calculate Spearman correlations, kurtosis, skewness, CV, and Gini for each measurement column
    df_by_well[f'Spearman_cor'] = df_by_well.apply(lambda row: spearman_correlation(row[mcols['Marker']], row[mcols['HPA']]), axis=1)
    df_by_well[f'Spearman_pval'] = df_by_well.apply(lambda row: spearman_pval(row[mcols['Marker']], row[mcols['HPA']]), axis=1)
    for chan,col in mcols.items():
        df_by_well[f'Kurtosis_{chan}'] = df_by_well[col].apply(calculate_kurtosis)
        df_by_well[f'Skewness_{chan}'] = df_by_well[col].apply(calculate_skewness)
        df_by_well[f'CV_{chan}'] = df_by_well[col].apply(cv)
        df_by_well[f'Gini_{chan}'] = df_by_well[col].apply(gini)

    return df_by_well

## Plotting functions on well dataframe
def plot_heatmap_correlation(df,savename='heatmap'):
    df_heatmap = df.pivot_table(index=['HPA_gene'],columns=['Marker'],values=['Spearman_cor'])
    #print(df_heatmap)
    df_heatmap = df_heatmap.fillna(0)
    #print(df_heatmap.isnull().sum())
    sns.set(font_scale=1.5)
    heatmap = sns.clustermap(figsize = (30,80),data=df_heatmap,cmap='YlGnBu',xticklabels = True,yticklabels = True)
    #heatmap.set(xlabel=None)
    plt.savefig(f'{savename}.png')
    plt.close() 

def plot_distribution_by_target(df,m='Gini',annotate=False):
    '''Plot the distribution of (Gini, CV, Kurtosis, Skewness) in different channels '''
    channels=['DAPI','Marker','HPA','Tubulin']
    col_list=[f'{m}_{chan}' for chan in channels]
    data = df.melt(id_vars=['Well_id','Marker'],value_vars=col_list,var_name='Channel',value_name=m)
    plot_distribution(data,col_x='Channel',col_y=m,savename=f'{m}_by_channel',annotate=annotate)

def plot_scatter(df,x,y,hue='Marker',addline=True):

    plt.figure(figsize=(8,8))
    sns.set_theme()
    #corr_plot = sns.regplot(data=df,x=x,y=y,scatter_kws={'s':0.5},line_kws={'lw':1})
    corr_plot = sns.lmplot(data=df,x=x,y=y,hue=hue,scatter_kws={'s':0.7},line_kws={'lw':1},fit_reg=False)
    #corr_plot.legend(loc='upper right',prop={'size': 10})
    #corr_plot.set(xlim = (0,100))
    #corr_plot.set(ylim = (0,10))
    #corr_plot.set(xscale = "log2")
    #corr_plot.set(yscale = "log2")#Use lmplot to plot scatter points
    corr_plot.set(xlabel = x)
    corr_plot.set(ylabel = y)
    if addline:
        plt.axline((0,0), slope=1,linestyle='--')
    filename = f'{x}_{y}.png'
    plt.savefig(filename)
    plt.close()