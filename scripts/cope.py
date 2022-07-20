import argparse
import pandas as pd
import itertools 
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule
import numpy as np
import networkx as nx 
from random import choice
from random import shuffle

SINGLE = Chem.BondType.SINGLE
DOUBLE = Chem.BondType.DOUBLE
TRIPLE = Chem.BondType.TRIPLE
AROMATIC = Chem.BondType.AROMATIC



parser = argparse.ArgumentParser()
parser.add_argument('--smiles', help='Smiles to find cope rearrangments on', type=str)
parser.add_argument('--write_smiles', help='write SMILES to filename', type=str)
parser.add_argument('--write_xyz', help='Write xyz files for unique molecules', type=str)


def copy_atom(atom):
    new_atom = Chem.Atom(atom.GetSymbol())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetAtomMapNum(atom.GetAtomMapNum())
    return new_atom


def mol_to_edit(mol, indices, rxn='cope'):
    emol = Chem.EditableMol(mol)
    if rxn == 'cope':
        pairs =[(indices[k], indices[k+1]) for k in range(len(indices)-1)]
        pairs.append((indices[-1], indices[0])) # Pairs to edit bonds of 
        bonds = [mol.GetBondBetweenAtoms(pair[0], pair[1]) for pair in pairs] # Associated existing or to exist bonds 
        order = False # Go down 
        for (bond, pair) in zip(bonds, pairs):
            if not order:
                # Should never have blank bond be at False step 
                if bond is None:
                    print("Error: Cannot reduce bond order of non-existing bond ... exiting.")
                    exit()
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    emol.RemoveBond(pair[0], pair[1])
                elif bond.GetBondType() == Chem.BondType.DOUBLE:
                    emol.RemoveBond(pair[0], pair[1])
                    emol.AddBond(pair[0], pair[1], Chem.BondType.SINGLE)
                else:
                    print("Tried to remove a bond that wasn't none, single, or double, something is wrong ... exiting")
                    exit()
            if order:
                if bond is None:
                    emol.AddBond(pair[0], pair[1], Chem.BondType.SINGLE)
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    emol.RemoveBond(pair[0], pair[1])
                    emol.AddBond(pair[0], pair[1], Chem.BondType.DOUBLE)
                elif bond.GetBondType() == Chem.BondType.DOUBLE:
                    print("Tried to make an existing double bond a triple ... not valid for cope structures ... exiting")
                    exit()
                else:
                    print("Tried to add to a bond that wasn't none, single, or double, something is wrong ... exiting")
                    exit()
            order = not order # Flip 
    return emol.GetMol()





args = parser.parse_args()
if args.smiles is None:
    raise ValueError('SMILES string not specified by user ... need molecule to search on!')

molecule = Chem.MolFromSmiles(args.smiles)
molecule = Chem.AddHs(molecule)
file_smiles = args.write_smiles
file_xyz = args.write_xyz

# Three structures that enable cope rearrangments 
substruct_1 = Chem.MolFromSmarts('C=CCCC=C') #Structure one reacts to reform itself 
#substruct_3 = Chem.MolFromSmarts('C=CC=CC=C')


# Get all cope rearrangments
def get_copes(mol):
    # Three structures that enable cope rearrangments 
    substruct_1 = Chem.MolFromSmarts('C=CCCC=C') #Structure one reacts to reform itself 
    substruct_3 = Chem.MolFromSmarts('C=CC=CC=C')
    cope_2_pi = None
    if molecule.HasSubstructMatch(substruct_1):
        cope_2_pi = molecule.GetSubstructMatches(substruct_1)

    cope_3_pi = None
    if molecule.HasSubstructMatch(substruct_3):
        cope_3_pi = molecule.GetSubstructMatches(substruct_3)

    all_copes = []
    if cope_2_pi:
        all_copes.extend(cope_2_pi)
    if cope_3_pi:
        all_copes.extend(cope_3_pi)  
    return all_copes

unique_molecules, start, N = [Chem.CanonSmiles(Chem.MolToSmiles(molecule))], 0, 1
# With the reacting substructures in the molecule, modify molecule to correspond to rearrangment product
search = True
while search:
    search = False # Assume no new molecules will be found
    update = True # Update start 
    for n in range(start, N):
        cope_rxns = get_copes(molecule)
        if cope_rxns is not None:
            #print('Cope rearrangments found! ... Checking reactions generated for unique molecule: ', n + 1)
            for atoms in cope_rxns:
                molecule_mod = mol_to_edit(molecule, atoms)
                canonical_smiles = Chem.CanonSmiles(Chem.MolToSmiles(molecule_mod))
                if canonical_smiles not in unique_molecules:
                    N += 1
                    #print("New molecule found! ... Current unique molecules:", N)
                    # A new molecule found! Continue search from however many starting indexes are found  
                    unique_molecules.append(canonical_smiles)
                    if update:
                        start += 1
                        update = False
                    search = True 
            if not search:
                #print('No new molecules found for unique molecule:', n + 1)
                pass

print('Total molecules found:',len(unique_molecules))
if file_smiles is not None:
    print('Writing SMILES to csv')
    pd.DataFrame(unique_molecules, columns=['SMILES']).to_csv(file_smiles + '.csv')

if file_xyz is not None:
    print('Generating and writing xyz coordinates for molecules')
    for n, smiles in enumerate(unique_molecules):
        path = file_xyz + '_' + str(n) + '.xyz'
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        EmbedMolecule(mol)
        Chem.MolToXYZFile(mol,filename=path)
