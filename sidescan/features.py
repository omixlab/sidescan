from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def mol_to_morgan_fingerprints(mol):
    fingerprints = AllChem.GetMorganFingerprintAsBitVect(mol, useChirality=True, radius=2, nBits = 1024, bitInfo={})
    return np.array(fingerprints)

def smiles_to_morgan_fingerprints(smiles):
    mol = Chem.rdmolfiles.MolFromSmiles(smiles)
    return mol_to_morgan_fingerprints(mol)
