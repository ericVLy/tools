from rdkit import Chem
from rdkit.Chem import Draw
import argparse



def show_mol_from_file(mol_file:str, size:tuple, kekulize:bool):
    mol = Chem.MolFromMolFile(mol_file)
    return Draw.ShowMol(mol=mol, size=size, kekulize=kekulize)


def get_smiles_from_mol(mol_file:str):
    mol = Chem.MolFromMolFile(mol_file)
    return Chem.MolToSmiles(mol)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='show mol file')
    parser.add_argument('--mol_file', type=str, help='mol file full path')
    parser.add_argument('--height', type=int, help='image height')
    parser.add_argument('--length', type=int, help='image length')
    parser.add_argument('--kekulize', type=bool, help='kekulize')
    args = parser.parse_args()
    print(get_smiles_from_mol(mol_file=args.mol_file))
    pic = show_mol_from_file(mol_file=args.mol_file, size=(args.length, args.height,), kekulize=args.kekulize)
    pic.show()