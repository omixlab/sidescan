from sidescan.cli import preprocess_argument_parser
from sidescan.features import smiles_to_morgan_fingerprints
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import pandas as pd
import numpy as np
import os

def main():

    arguments = preprocess_argument_parser.parse_args()
    drug_names_and_smiles_filepath = os.path.join(arguments.directory, 'drug_names_and_smiles.csv')
    if not os.path.isfile(drug_names_and_smiles_filepath):
        print(f'file drug_names_and_smiles.csv not found in {arguments.directory}')
        exit()
    
    df_freq = pd.read_csv(os.path.join(arguments.directory, 'meddra_freq.tsv.gz'), header=None, sep='\t', compression='gzip')
    df_freq = df_freq[df_freq[3] != 'placebo']

    df_labels = pd.read_csv(os.path.join(arguments.directory, 'meddra_all_label_se.tsv.gz'), header=None, sep='\t', compression='gzip')

    selected_labels = list(df_labels[6].unique())

    df_drug_names_and_smiles = pd.read_csv(drug_names_and_smiles_filepath)

    features = []
    y = []

    for r, row in tqdm(df_drug_names_and_smiles.iterrows(), total=df_drug_names_and_smiles.shape[0]):
        
        labels = list(df_labels[(df_labels[1] == row.CID) | (df_labels[2] == row.CID)][6].unique())
        labels = [label for label in labels if df_freq[(df_freq[9] == label) & (df_freq[5] >= 0.05)].shape[0] > 0]
        
        fingerprints = smiles_to_morgan_fingerprints(row.SMILES)
        molecule_labels = [1 if label in labels else 0 for label in selected_labels]
        y.append(molecule_labels)
        features.append(fingerprints)
    
    pd.DataFrame(features).to_csv(os.path.join(arguments.directory, 'fingerprints.csv'), index=False)
    pd.DataFrame(y,columns=selected_labels).to_csv(os.path.join(arguments.directory, 'labels.csv'), index=False)
    

if __name__ == '__main__':
    main()