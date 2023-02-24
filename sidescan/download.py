from sidescan.cli import download_argument_parser
from tqdm import tqdm
import pandas as pd
import requests
import json
import os

TEXT_FILES = [
    ('drug_names.tsv', 'http://sideeffects.embl.de/media/download/drug_names.tsv'),
    ('drug_atc.tsv', 'http://sideeffects.embl.de/media/download/drug_atc.tsv')
]
    
BINARY_FILES = [
    ('meddra_all_indications.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra_all_indications.tsv.gz'),
    ('meddra_all_se.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz'),
    ('meddra_freq.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra_freq.tsv.gz'),
    ('meddra_all_label_indications.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra_all_label_indications.tsv.gz'),
    ('meddra_all_label_se.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra_all_label_se.tsv.gz'),
    ('meddra.tsv.gz', 'http://sideeffects.embl.de/media/download/meddra.tsv.gz')
]

def main():
    
    arguments = download_argument_parser.parse_args()

    print(f'Downloading SIDERS file to {os.path.abspath(arguments.directory)}/ (it may take a while ...)')

    for filename, url in TEXT_FILES:
        download_text_files(url, os.path.join(arguments.directory, filename))
    for filename, url in BINARY_FILES:
        download_binary_files(url, os.path.join(arguments.directory, filename))

    df_drug_names = pd.read_csv(f'{arguments.directory}/drug_names.tsv', header=None, sep='\t')
    drug_names_and_smiles = []

    print(f'Downloading PubChem SMILES {os.path.abspath(arguments.directory)}/drug_names_and_smiles.csv (it may take a while ...)')

    for r, row in tqdm(df_drug_names.iterrows(), total=df_drug_names.shape[0]):
        cid = row[0][4:].lstrip('0')
        data = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON').text
        drug_names_and_smiles.append({
            'CID': row[0],
            'CID_CLEAN': cid,
            'SMILES': json.loads(data)['PropertyTable']['Properties'][0]['CanonicalSMILES']
        })
    
    df_drug_names_and_smiles = pd.DataFrame(drug_names_and_smiles)
    df_drug_names_and_smiles.to_csv(f'{arguments.directory}/drug_names_and_smiles.csv', index=False)

def download_text_files(url, filepath):
    with open(filepath, 'w') as writer:
        r = requests.get(url)
        writer.write(r.text)
    print(filepath, 'done!')

def download_binary_files(url, filepath):
    with open(filepath, 'wb') as writer:
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                writer.write(chunk)
    print(filepath, 'done!')

if __name__ == '__main__':
    main()