import cellpose_analysis_utils as utils
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import iqr
import numpy as np

'''
Index(['Well_id', 'Marker', 'HPA_Ab', 'HPA_gene', 'Plate', 'Nuc_MeanInt_DAPI',
       'Cyto_MeanInt_Marker', 'Cyto_MeanInt_HPA', 'Cyto_MeanInt_Tubulin',
       'Spearman_cor', 'Spearman_pval', 'Kurtosis_DAPI', 'Skewness_DAPI',
       'CV_DAPI', 'Gini_DAPI', 'Kurtosis_Marker', 'Skewness_Marker',
       'CV_Marker', 'Gini_Marker', 'Kurtosis_HPA', 'Skewness_HPA', 'CV_HPA',
       'Gini_HPA', 'Kurtosis_Tubulin', 'Skewness_Tubulin', 'CV_Tubulin',
       'Gini_Tubulin'],
      dtype='object')'''

# Add additional measurements
def add_pearson(df,m='Cyto_MeanInt'):
      df['Pearson_cor'] = df.apply(lambda row: utils.pearson_correlation(row[f'{m}_HPA'],row[f'{m}_Marker']),axis=1)
      df['Pearson_pval'] = df.apply(lambda row: utils.pearson_pval(row[f'{m}_HPA'],row[f'{m}_Marker']),axis=1)
      return df
def add_number_object(df,m='Cyto_MeanInt'):
      df['Total_objects'] = df[f'{m}_Marker'].apply(len)
      return df

def get_extremes(marker, hpa,percentage = 10):
      top_mask = ((marker>=np.percentile(marker,100-percentage)) & (marker<=np.percentile(marker,99)))
      bottom_mask = ((marker>=np.percentile(marker,1)) & (marker<=np.percentile(marker,percentage)))
      top_cells_marker = marker[top_mask]
      top_cells_hpa = hpa[top_mask]
      bottom_cells_marker = marker[bottom_mask]
      bottom_cells_hpa = hpa[bottom_mask]
      return top_cells_marker,top_cells_hpa,bottom_cells_marker,bottom_cells_hpa
def add_spearman_extremes(df,m='Cyto_MeanInt',p = 10):
      df[f'Marker_top{p}'] =df.apply(lambda row:get_extremes(row[f'{m}_Marker'],row[f'{m}_HPA'],percentage=p)[0],axis=1)
      df[f'HPA_top{p}'] = df.apply(lambda row:get_extremes(row[f'{m}_Marker'],row[f'{m}_HPA'],percentage=p)[1],axis=1)
      df[f'Marker_bottom{p}'] = df.apply(lambda row:get_extremes(row[f'{m}_Marker'],row[f'{m}_HPA'],percentage=p)[2],axis=1)
      df[f'HPA_bottom{p}'] = df.apply(lambda row:get_extremes(row[f'{m}_Marker'],row[f'{m}_HPA'],percentage=p)[3],axis=1)
      df[f'Spearman_top{p}']= df.apply(lambda row: utils.spearman_correlation(row[f'Marker_top{p}'], row[f'HPA_top{p}']), axis=1)
      df[f'Spearman_top{p}_pval']= df.apply(lambda row: utils.spearman_pval(row[f'Marker_top{p}'], row[f'HPA_top{p}']), axis=1)
      df[f'Spearman_bottom{p}']= df.apply(lambda row: utils.spearman_correlation(row[f'Marker_bottom{p}'], row[f'HPA_bottom{p}']), axis=1)
      df[f'Spearman_bottom{p}_pval']= df.apply(lambda row: utils.spearman_pval(row[f'Marker_bottom{p}'], row[f'HPA_bottom{p}']), axis=1)
      return df
# Normalization
def z_score(arr):
      z_scores = (arr - np.mean(arr))/np.std(arr)
      return z_scores
def log_scale(arr):
      return np.log2(arr)
def clipping(arr):
      arr_clipped = np.copy(arr)
      iqr_value = iqr(arr)
      lower = np.percentile(arr,25) - 1.5* iqr_value
      upper = np.percentile(arr,75) + 1.5* iqr_value
      outliers_low = (arr < lower) 
      outliers_high = (arr > upper)
      #arr_clipped[outliers_low] = lower 
      #arr_clipped[outliers_high] = upper
      arr_clipped = arr[~outliers_low & ~outliers_high]
      return arr_clipped

def min_max(arr):
      min_max_scaled = (arr -np.min(arr))/(np.max(arr) - np.min(arr))
      return min_max_scaled
def normalize_intensities(df,m='MeanInt'):
      df[f'{m}_HPA_scaled'] = df[f'{m}_HPA'].apply(lambda x: min_max(log_scale(x)))
      df[f'{m}_Marker_scaled'] = df[f'{m}_Marker'].apply(lambda x: min_max(log_scale(x)))
      df['Spearman_cor_scaled'] = df.apply(lambda row: utils.spearman_correlation(row[f'{m}_HPA_scaled'],row[f'{m}_Marker_scaled']),axis=1)
      df['Pearson_cor_scaled'] = df.apply(lambda row: utils.pearson_correlation(row[f'{m}_HPA_scaled'],row[f'{m}_Marker_scaled']),axis=1)
      return df

# Plot distribution of statistics across sample
# measurement by compartment
def plot_by_compartment(df):
      utils.plot_distribution_by_target(df,m='Gini',annotate=True)
      utils.plot_distribution_by_target(df,m='CV',annotate=True)
      utils.plot_distribution_by_target(df,m='Kurtosis',annotate=True)
      utils.plot_distribution_by_target(df,m='Skewness',annotate=True)

# measurement by marker
def plot_by_marker(df):
      utils.plot_distribution(df,col_x='Marker',col_y='CV_Marker',annotate=True,savename='CV_by_marker')
      utils.plot_distribution(df,col_x='Marker',col_y='Gini_Marker',annotate=True,savename='Gini_by_marker')
      utils.plot_distribution(df,col_x='Marker',col_y='Skewness_Marker',annotate=True,savename='Skewness_by_marker')
      utils.plot_distribution(df,col_x='Marker',col_y='Kurtosis_Marker',annotate=True,savename='Kurtosis_by_marker')

# Filtering samples
def get_mask(df,col='Marker',value='IDH3A'):
      return df[col]==value

def apply_filters(df):
      obj_num_mask = df['Total_objects']>150
      pval_mask = df['Spearman_pval'] < 0.01
      corr_mask = df['Spearman_cor'] >= 0.5 
      variable_mask = (df['CV_HPA'] > df['CV_DAPI'].mean()) & (df['CV_Marker'] > df['CV_DAPI'].mean())
      df_filtered = df[obj_num_mask&pval_mask&variable_mask&corr_mask]
      return df_filtered


# Annotation
def apply_pathway_groups(df):
      pathway_map = pd.read_csv('/home/xiuqi/test.csv')
      df = df.merge(pathway_map,on=['HPA_Ab','HPA_gene'],how='left')
      df.drop_duplicates(subset = ['HPA_Ab','HPA_gene','PathwayGroup'],keep='first',inplace=True)
      return df

def plot_countplot(df,x,hue=None,savename='count'):
      sns.countplot(data=df,x=x,hue=hue)
      plt.savefig(f'{savename}.png')
      plt.close()

def plot_scatter_correlation(df_batch,gene):
      df_batch = utils.drop_duplicate_cellIDs(df_batch,savedrop=False)
      df_batch = utils.drop_out_nucleus(df_batch,savedrop=False)
      df_batch = utils.filter_by_max_intensity(df_batch,['Cell_Intensities_DAPI','Cell_Intensities_Marker','Cell_Intensities_Tubulin'])
      df_batch = utils.HPA_column_selection(df_batch)
      df_batch = utils.filter_by_max_intensity(df_batch,['Cell_Intensities_HPA'])
      df_batch = utils.add_exp_info(df_batch)
      df_batch['Mito_MeanInt_Marker']=df_batch.apply(lambda row:utils.calculate_mito_mean(row['Cyto_Intensities_Marker'],row['Cyto_Intensities_HPA'])[0],axis=1)
      df_batch['Mito_MeanInt_HPA']=df_batch.apply(lambda row:utils.calculate_mito_mean(row['Cyto_Intensities_Marker'],row['Cyto_Intensities_HPA'])[1],axis=1)
      utils.plot_scatter(df_batch[df_batch['HPA_gene']==gene],x='Mito_MeanInt_Marker',y='Mito_MeanInt_HPA',hue=None,addline=False)

def plot_scatter_int(df,wellid,m='Cyto_MeanInt'):
      plt.figure(figsize=(6,6))
      #plt.scatter(df.loc[df['Well_id'] ==wellid,f'{m}_Marker'].values[0],df.loc[df['Well_id'] ==wellid,f'{m}_HPA'].values[0],s=0.5)
      plt.scatter(df.loc[df['Well_id'] ==wellid,f'{m}_Marker_scaled'].values[0],df.loc[df['Well_id'] ==wellid,f'{m}_HPA_scaled'].values[0],s=0.5)
      plt.xlabel(df.loc[df['Well_id'] ==wellid,'Marker'].values[0])
      plt.ylabel(df.loc[df['Well_id'] ==wellid,'HPA_gene'].values[0])
      plt.title(f'{m} (log scaled)')
      plt.savefig(f'{wellid}_{m}.png')
      plt.close()

def plot_hist_tubulin(combined_array,name):
      #print(np.max(combined_array))
      #print(np.min(combined_array))
      #print(np.percentile(combined_array,1))
      #print(np.percentile(combined_array,99))
      plt.hist(combined_array,bins=50)
      #plt.xlim(right=0.5)
      plt.axvline(np.mean(combined_array),linestyle='--',color='red')
      plt.axvline(np.median(combined_array),linestyle='--',color='blue')
      plt.axvline(np.percentile(combined_array,1),linestyle='--',color='blue')
      plt.axvline(np.percentile(combined_array,25),linestyle='--',color='blue')
      plt.axvline(np.percentile(combined_array,75),linestyle='--',color='blue')
      plt.axvline(np.percentile(combined_array,99),linestyle='--',color='blue')
      plt.xlabel('Cyto_UpperQuartileIntensity_Tubulin')
      plt.title(name)
      plt.savefig(f'{name}.png')
      plt.close()
def plot_tubulin(df,marker):
      targets = [f'{marker}_1',f'{marker}_2',f'{marker}_3',f'{marker}_4',f'{marker}_5']
      arrays =df[get_mask(df,'Marker',marker)]['Cyto_Q3_Tubulin'].values
      combined_array = np.concatenate(arrays)
      plot_hist_tubulin(combined_array,marker)
      for target in targets:
            arrays =df[get_mask(df,'Plate',target)]['Cyto_Q3_Tubulin'].values
            combined_array = np.concatenate(arrays)
            plot_hist_tubulin(combined_array,target)
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
def plot_heatmap_extreme(df,savename='heatmap'):
      df_heatmap = df[df['Marker']=='TUFM'][['HPA_gene','Spearman_top15','Spearman_bottom15','Spearman_cor']]
      df_heatmap=df_heatmap.set_index('HPA_gene')
      print(df_heatmap.head())
      sns.set(font_scale=1.5)
      heatmap = sns.clustermap(figsize = (30,80),data=df_heatmap,cmap='YlGnBu',xticklabels = True,yticklabels = True)
      #heatmap.set(xlabel=None)
      plt.savefig(f'{savename}.png')
      plt.close()


if __name__ == '__main__':
      measurement='Mito_MeanInt'
      iq_filter = 'unfiltered'      #imagequality_filtered_thres1
      # 'Cyto_MeanInt', 'Cyto_Q3','Cyto_TopQuartileMean','Mito_MeanInt'
      directory = "/mnt/hddmount/cellpose_analysis"
      df = pd.read_feather(os.path.join(directory,f'{measurement}_{iq_filter}.feather'))
      #print(df.columns)
      df = add_spearman_extremes(df,m=measurement,p=15)
 
      df_filtered = df[(df['Spearman_top15_pval']<0.05)|(df['Spearman_bottom15_pval']<0.05)|(df['Spearman_cor']>0.5)]
      plot_heatmap_extreme(df_filtered,savename='heatmap_TUFM_test')
      '''
      # Read in files
      for m in ['Cyto_MeanInt', 'Cyto_Q3','Cyto_TopQuartileMean','Mito_MeanInt']:

            measurement=m
            iq_filter = 'unfiltered'      #imagequality_filtered_thres1
            # 'Cyto_MeanInt', 'Cyto_Q3','Cyto_TopQuartileMean','Mito_MeanInt'
            directory = "/mnt/hddmount/cellpose_analysis"
            df = pd.read_feather(os.path.join(directory,f'{measurement}_{iq_filter}.feather'))

            df = add_number_object(df,m=measurement)
            df = normalize_intensities(df,m=measurement)

            df_filtered= apply_filters(df)
            print(df_filtered.head())
            utils.plot_heatmap_correlation(df_filtered,savename=f'{measurement}_231123_4')
      #df_filtered = apply_pathway_groups(df_filtered)
      '''
      '''
      col_list = ['Well_id','Marker','HPA_Ab','HPA_gene','Spearman_cor','Spearman_pval',
      'Kurtosis_DAPI','Skewness_DAPI','CV_DAPI','Gini_DAPI',
      'Kurtosis_Marker','Skewness_Marker','CV_Marker','Gini_Marker',
      'Kurtosis_HPA','Skewness_HPA','CV_HPA','Gini_HPA',
      'Kurtosis_Tubulin','Skewness_Tubulin','CV_Tubulin','Gini_Tubulin',
      'Total_objects','PathwayGroup','PathwayGroup_number','Multiple'
      ]
      df_filtered.sort_values(by=['Marker','Spearman_cor'],ascending=False)[col_list].to_csv('Cyto_TopQuartileMean_thres2_sorted.csv')
      '''
      '''
      for marker in ['IDH3A','TUFM','GLS','SOD2']:
            df_marker = df_filtered[df_filtered['Marker']==marker]
            df_marker = apply_pathway_groups(df_marker)
            plot_countplot(df_marker,x='PathwayGroup_number',hue='Multiple',savename=f'{marker}_pathway_count')
            plot_countplot(df_marker,x='Multiple',savename=f'{marker}_multiplepathway_count')
      '''
      '''
      plot_by_compartment(df)
      plot_by_marker(df)
      
      for marker in ['IDH3A','SOD2','GLS','TUFM']:
            utils.plot_histogram(df[df['Marker']==marker],col='Total_objects',savename=f'{marker}_objects_per_sample_thres1.png',bins=50)
      
      
      print(df[df['HPA_gene']=='IDH3A'])
      print(df[df['HPA_gene']=='DLST'])
      '''
      '''
      plot_scatter_int(df,'IDH3A_3_E01',m=measurement)
      plot_scatter_int(df,'IDH3A_4_A04',m=measurement)
      plot_scatter_int(df,'TUFM_3_E01',m=measurement)
      plot_scatter_int(df,'SOD2_3_E01',m=measurement)
      plot_scatter_int(df,'IDH3A_1_B03',m=measurement)
      plot_scatter_int(df,'SOD2_1_B03',m=measurement)
      
      
      #plot_scatter_int(df,'IDH3A_4_A05',m=measurement)
      col_list = ['Well_id', 'Marker', 'HPA_Ab', 'HPA_gene', 
       'Spearman_cor', 'Spearman_pval', 
       #'Kurtosis_DAPI', 'Skewness_DAPI','CV_DAPI', 'Gini_DAPI', 
       'Kurtosis_Marker', 'Skewness_Marker','CV_Marker', 'Gini_Marker', 
       'Kurtosis_HPA', 'Skewness_HPA', 'CV_HPA','Gini_HPA', 
       'Kurtosis_Tubulin', 'Skewness_Tubulin', 'CV_Tubulin','Gini_Tubulin', 
       'Total_objects', 'Pearson_cor_scaled']
      df_filtered = apply_filters(df)
      df_heatmap = df_filtered.pivot_table(index=['HPA_gene'],columns=['Marker'],values=['Spearman_cor'])
      df_heatmap = df_heatmap.fillna(0)
      print(df_heatmap)
      print(df_heatmap.shape)
      print(df_heatmap.columns)
      '''
      utils.plot_heatmap_correlation(df_filtered[df_filtered['Marker'] != 'TUFM'])

      
