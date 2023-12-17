import cellpose_analysis_utils as utils
import cellpose_analysis_image_quality as imagequality
import pandas as pd
import os



main_directory = "/mnt/hddmount/cellpose_output"
output_directory = "/mnt/hddmount/cellpose_analysis"

for measurement in ['Cyto_MeanInt', 'Cyto_Q3','Cyto_TopQuartileMean','Mito_MeanInt']:
    df = utils.create_final_df(main_directory,measurement=measurement)
    df_corrected = utils.correct_duplicates(df,m=measurement)
    del df
    df_corrected.to_feather(os.path.join(output_directory,f'{measurement}_imagequality_filtered_thres1.feather'))
    del df_corrected
