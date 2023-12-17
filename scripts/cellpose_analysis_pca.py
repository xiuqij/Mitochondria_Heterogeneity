import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
import plotly.express as px
import plotly.io as pio
import cellpose_analysis_image_quality as imagequality

directory = "/mnt/hddmount/cellpose_analysis"
m = 'Cyto_Q3'
df = pd.read_feather(os.path.join(directory,f'{m}_imagequality_filtered_thres1.feather'))

'''
col_list =['FileName_DAPI',
'ImageQuality_Correlation_DAPI_30', 'ImageQuality_Correlation_HPA_H_30', 'ImageQuality_Correlation_HPA_L_30', 'ImageQuality_Correlation_Marker_30', 'ImageQuality_Correlation_Tubulin_30', 
'ImageQuality_FocusScore_DAPI', 'ImageQuality_FocusScore_HPA_H', 'ImageQuality_FocusScore_HPA_L', 'ImageQuality_FocusScore_Marker', 'ImageQuality_FocusScore_Tubulin', 
'ImageQuality_LocalFocusScore_DAPI_30', 'ImageQuality_LocalFocusScore_HPA_H_30', 'ImageQuality_LocalFocusScore_HPA_L_30', 'ImageQuality_LocalFocusScore_Marker_30', 'ImageQuality_LocalFocusScore_Tubulin_30',
'ImageQuality_PercentMaximal_DAPI', 'ImageQuality_PercentMaximal_HPA_H', 'ImageQuality_PercentMaximal_HPA_L', 'ImageQuality_PercentMaximal_Marker', 'ImageQuality_PercentMaximal_Tubulin', 
'ImageQuality_PercentMinimal_DAPI', 'ImageQuality_PercentMinimal_HPA_H', 'ImageQuality_PercentMinimal_HPA_L', 'ImageQuality_PercentMinimal_Marker', 'ImageQuality_PercentMinimal_Tubulin', 
'ImageQuality_PowerLogLogSlope_DAPI', 'ImageQuality_PowerLogLogSlope_HPA_H', 'ImageQuality_PowerLogLogSlope_HPA_L', 'ImageQuality_PowerLogLogSlope_Marker', 'ImageQuality_PowerLogLogSlope_Tubulin', 
]
df_iq = imagequality.read_df(directory,['TUFM_1_Image.csv','TUFM_2_Image.csv','TUFM_3_Image.csv','TUFM_4_Image.csv','TUFM_5_Image.csv'])
df_iq = imagequality.mark_image(df_iq)
print(df.shape)
toremove = df_iq[df_iq['Remove']==True]['Well_id'].tolist()
df = df[~df['Well_id'].isin(toremove)]
print(df.shape)
'''

features =[
        'Int_DAPI', 'Int_Tubulin', 'Int_Marker',
        'CV_DAPI', 'CV_Tubulin', 'CV_Marker',
        'Gini_DAPI', 'Gini_Tubulin', 'Gini_Marker',
        'Skewness_DAPI', 'Skewness_Tubulin', 'Skewness_Marker',
        'Kurtosis_DAPI', 'Kurtosis_Tubulin','Kurtosis_Marker'   
        ]


def df_prep(df,m='Cyto_MeanInt'):
    if m=='Cyto_MeanInt' or m=='Mito_MeanInt':
        df['Int_DAPI'] = df['Nuc_MeanInt_DAPI'].apply(np.mean)
    else:
        df['Int_DAPI'] = df[f"Nuc_{m.split('_')[1]}_DAPI"].apply(np.mean)
    df['Int_Tubulin'] = df[f'{m}_Tubulin'].apply(np.mean)
    df['Int_Marker'] = df[f'{m}_Marker'].apply(np.mean)
    return df
df = df_prep(df,m=m)

def pca_marker(df,marker):   
    features =[
        'Int_DAPI', 'Int_Tubulin', 'Int_Marker',
        'CV_DAPI', 'CV_Tubulin', 'CV_Marker',
        'Gini_DAPI', 'Gini_Tubulin', 'Gini_Marker',  
        ]
    df_marker = df[df['Marker']==marker].reset_index()
    x = df_marker.loc[:,features].values
    #y = df_marker.loc[:,['Plate']].values
    #print(df_marker.head())
    
    scaled_df = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents,columns=['principal component 1','principal component 2'])
    finalDf = pd.concat([principalDf,df_marker[['Plate','Int_DAPI', 'Int_Tubulin', 'Int_Marker',
        'CV_DAPI', 'CV_Tubulin', 'CV_Marker',
        'Gini_DAPI', 'Gini_Tubulin', 'Gini_Marker']]],axis=1)
    

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    

    targets = [f'{marker}_1', f'{marker}_2', f'{marker}_3',f'{marker}_4',f'{marker}_5']
    colors = ['red', 'green', 'blue','yellow','purple']
    for target, color in zip(targets,colors):
        indicesToKeep = finalDf['Plate'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                , finalDf.loc[indicesToKeep, 'principal component 2']
                , c = color
                , s = 5)
    #ax.set_xlim(0,100)
    #ax.set_ylim(0,100)
    ax.legend(targets)
    ax.grid()
    plt.savefig(f'pca_{marker}.png')
    plt.close()

    for feature in features:
        plt.figure(figsize = (8,8))
        plt.xlabel('Principal Component 1', fontsize = 15)
        plt.ylabel('Principal Component 2', fontsize = 15)
        plt.title('2 component PCA', fontsize = 20)
        plt.scatter(finalDf['principal component 1']
                , finalDf['principal component 2']
                , c = finalDf[feature]
                , s = 5
                , cmap='viridis')
        #ax.grid()
        plt.colorbar()
        plt.savefig(f'pca_{marker}_{feature}.png')
        plt.close()


    print(pca.explained_variance_ratio_)

for marker in ['IDH3A','SOD2','GLS','TUFM']:
    pca_marker(df,marker)


'''
def get_mask(df,col='Marker',value='IDH3A'):
      return df[col]==value
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
      arrays =df[df['Marker']==marker]['Cyto_Q3_Tubulin'].values
      combined_array = np.concatenate(arrays)
      plot_hist_tubulin(combined_array,marker)
      for target in targets:
            arrays =df[get_mask(df,'Plate',target)]['Cyto_Q3_Tubulin'].values
            combined_array = np.concatenate(arrays)
            plot_hist_tubulin(combined_array,target)
plot_tubulin(df,'TUFM')


features =[
    #'Int_DAPI', 'Int_Tubulin', 'Int_Marker',
    'CV_DAPI', 'CV_Tubulin', 'CV_Marker',
    #'Gini_DAPI', 'Gini_Tubulin', 'Gini_Marker',
    #'Skewness_DAPI', 'Skewness_Tubulin', 'Skewness_Marker',
    #'Kurtosis_DAPI', 'Kurtosis_Tubulin','Kurtosis_Marker'
    ]
df_marker = df[df['Marker']=='SOD2'].reset_index()
pca = PCA()
components = pca.fit_transform(df_marker[features])
labels = {
    str(i): f"PC {i+1} ({var:.1f}%)"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
}

fig = px.scatter_matrix(
    components,
    labels=labels,
    dimensions=range(2),
    color=df_marker["Plate"]
)
fig.update_traces(diagonal_visible=False)
pio.write_image(fig,'pca.png')
'''