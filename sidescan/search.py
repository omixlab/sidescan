from sidescan.cli import search_argument_parser
from sidescan.logo import logo
from sidescan.features import mol_to_morgan_fingerprints
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import warnings
import json
import pickle
import os

warnings.simplefilter('ignore')

def main():

    print(logo)

    models = {}

    arguments = search_argument_parser.parse_args()

    with open(arguments.model, 'rb') as reader:
        models_data = pickle.loads(reader.read())

        for side_effect in models_data:
            if models_data[side_effect]['metrics']['f1'] >= arguments.minimum_f1_score:
                models[side_effect] = pickle.loads(open(models_data[side_effect]['path'], 'rb').read())

    fileext = os.path.splitext(arguments.input)[1]

    if fileext == '.mol2':
        mols = [Chem.MolFromMolFile(arguments.input)]

    elif fileext == '.sdf':
        mols = Chem.SDMolSupplier(arguments.input)
    
    results = []

    for m, mol in enumerate(mols):

        mol_data = mol_to_morgan_fingerprints(mol)

        results.append(
            {   
                'index': m,
                'InchI': Chem.MolToInchi(mol),
                'InchI_key': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'side_effects': {
                    side_effect: float(models[side_effect].predict_proba([mol_data])[:,1]) for side_effect in tqdm(models)
                }
            }
        )


    with open(arguments.output, 'w') as writer:
        writer.write(json.dumps(results, indent=4))

if __name__ == '__main__':
    main()
