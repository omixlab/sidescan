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

    arguments = search_argument_parser.parse_args()

    with open(arguments.model, 'rb') as reader:
        models = pickle.loads(reader.read())
    
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
                    side_effect: float(models[side_effect].predict_proba([mol_data])[:,1]) for side_effect in tqdm(models) if models[side_effect].f1 >= arguments.minimum_f1_score
                }
            }
        )


    with open(arguments.output, 'w') as writer:
        writer.write(json.dumps(results, indent=4))

if __name__ == '__main__':
    main()
