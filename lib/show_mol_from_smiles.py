"""draw chem pic by smiles code"""
import argparse
from rdkit import Chem
from rdkit.Chem import Draw



def show_mol_from_smiles(smiles:str, size:tuple, kekulize:bool):
    mol = Chem.MolFromSmiles(smiles)
    Draw.ShowMol(mol=mol, size=size, kekulize=kekulize)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='show mol from smiles')
    parser.add_argument('--smiles', "-s", type=str, help='mol file full path')
    parser.add_argument('--height', "-H", type=int, help='image height')
    parser.add_argument('--length', "-L", type=int, help='image length')
    parser.add_argument('--kekulize', "-K", type=bool, help='kekulize')
    args = parser.parse_args()
    show_mol_from_smiles(smiles=args.smiles,
                         size=(args.length, args.height,),
                         kekulize=args.kekulize)
