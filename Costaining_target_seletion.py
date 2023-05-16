#%% Import packages
import pandas as pd
import numpy as np

#%% Loading data files to dataframes
subcellular_location = ""  # Subcelluar location datafile (.tsv). The file can be downloaded from https://www.proteinatlas.org/about/download ("3 Subellular location data")
rna_celline = "" # File containing ranscript expression levels summarized per gene in different cell lines (.tsv). The file can be downloaded from https://www.proteinatlas.org/about/download ("23	RNA HPA cell line gene data")

subcellular = pd.read_csv(subcellular_location,sep="\t")
expression = pd.read_csv(rna_celline,sep="\t")

#%% Get all mitochodrial proteins
location_mask = (subcellular["Main location"].str.contains("Mitochondria")) | (subcellular["Additional location"].str.contains("Mitochondria"))
mito_proteins_raw = subcellular[location_mask]

#%% Filter out the uncertain proteins
reliability_mask = mito_proteins_raw["Reliability"] != "Uncertain"
mito_proteins = mito_proteins_raw[reliability_mask]

#%% Get the expression level of the proteins in U2OS cell line
isexpressed = []
pTPM = []
for index, row in mito_proteins.iterrows():
    gene = row["Gene"]
    expression_u2os = expression[(expression["Gene"] == gene) & (expression["Cell line"] == "U-2 OS")]
    if expression_u2os.shape[0] == 1:
        if expression_u2os["pTPM"].item() < 2:
            isexpressed.append(False)
            pTPM.append(expression_u2os["pTPM"].item())
        else:
            isexpressed.append(True)
            pTPM.append(expression_u2os["pTPM"].item())
    else:
        isexpressed.append(np.nan)
        pTPM.append(np.nan)

#%% Get variable mitochondria proteins 
variation_mask = (mito_proteins["Single-cell variation intensity"].str.contains("Mitochondria")) | (mito_proteins["Single-cell variation spatial"].str.contains("Mitochondria"))
variable_mito_proteins = mito_proteins[variation_mask]
variable_mito_proteins.to_csv("Variable mitochondria proteins.csv")
#%% Get stable mitochondria proteins
stable_mask = ~variation_mask
stable_mito_proteins = mito_proteins[stable_mask]
stable_mito_proteins.to_csv("Stable mitochondria proteins.csv")
#%% Stable mitochondria proteins but variable in other compartments
other_variable_proteins = stable_mito_proteins[(stable_mito_proteins["Single-cell variation intensity"].notna()) | (stable_mito_proteins["Single-cell variation spatial"].notna())]
other_variable_proteins.to_csv("Other variable proteins.csv")

#%% Make the dataframe for target selection
'''The dataframe for target selection, with the following columns:
Gene: (str) Ensembl ID of the gene 
Gene name: (str) Gene symbol 
Reliability: (str) the reliability score form HPA
SCV: (Bool) whether the protein contains single-cell variation annotation on HPA
SCV_mitochondria: (Bool) whether the protein contains single-cell variation annotation in mitochondria on HPA
expression_U2OS: (Bool) whether the protein is expressed in U2OS cell line (pTPM >=2)
pTPM_U2OS: (float) the protein expression level in U2OS cell line'''

d = {"Gene":mito_proteins["Gene"],
    "Gene name":mito_proteins["Gene name"],
    "Reliability":mito_proteins["Reliability"],
    "SCV":(mito_proteins["Single-cell variation intensity"].notna()) | (mito_proteins["Single-cell variation spatial"].notna()),
    "SCV_mitochondria":variation_mask,
    "expression_U2OS":isexpressed,
    "pTPM_U2OS":pTPM
    }

df = pd.DataFrame(data=d)
df.to_csv("Target_selection.csv",sep="\t")
