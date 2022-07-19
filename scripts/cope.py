import argparse
import itertools 
from rdkit import Chem
import numpy as np
import networkx as nx 
from random import choice
from random import shuffle

SINGLE = Chem.BondType.SINGLE
DOUBLE = Chem.BondType.DOUBLE
TRIPLE = Chem.BondType.TRIPLE
AROMATIC = Chem.BondType.AROMATIC


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



# bvl_bvl smiles 

bvl_bvl = Chem.MolFromSmiles('[C@H]12C=C[C@H]3[C@@H](C(C4=C[C@H]5[C@H]6C=C[C@@H]4C=C[C@@H]56)=C2)[C@H]3C=C1')
bvl_bvl = Chem.AddHs(bvl_bvl)

#bvl = Chem.MolFromSmiles('C1=CC2C3C2C=CC1C=C3')
#bvl = Chem.AddHs(bvl)

# Three structures that enable cope rearrangments 
substruct_1 = Chem.MolFromSmarts('C=CCCC=C') #Structure one reacts to reform itself 
substruct_3 = Chem.MolFromSmarts('C=CC=CC=C')


# Get all cope rearrangments
def get_copes(mol):
    # Three structures that enable cope rearrangments 
    substruct_1 = Chem.MolFromSmarts('C=CCCC=C') #Structure one reacts to reform itself 
    substruct_3 = Chem.MolFromSmarts('C=CC=CC=C')
    cope_2_pi = None
    if bvl_bvl.HasSubstructMatch(substruct_1):
        cope_2_pi = bvl_bvl.GetSubstructMatches(substruct_1)

    cope_3_pi = None
    if bvl_bvl.HasSubstructMatch(substruct_3):
        cope_3_pi = bvl_bvl.GetSubstructMatches(substruct_3)

    all_copes = []
    if cope_2_pi:
        all_copes.extend(cope_2_pi)
    if cope_3_pi:
        all_copes.extend(cope_3_pi)  
    return all_copes

unique_molecules, start, N = [Chem.CanonSmiles(Chem.MolToSmiles(bvl_bvl))], 0, 1
# With the reacting substructures in the molecule, modify molecule to correspond to rearrangment product
search = True
while search:
    search = False # Assume no new molecules will be found
    update = True # Update start 
    for n in range(start, N):
        cope_rxns = get_copes(bvl_bvl)
        if cope_rxns is not None:
            print('Cope rearrangments found! ... Checking reactions generated for unique molecule: ', n + 1)
            for atoms in cope_rxns:
                bvl_bvl_modified = mol_to_edit(bvl_bvl, atoms)
                canonical_smiles = Chem.CanonSmiles(Chem.MolToSmiles(bvl_bvl_modified))
                if canonical_smiles not in unique_molecules:
                    N += 1
                    print("New molecule found! ... Current unique molecules:", N)
                    # A new molecule found! Continue search from however many starting indexes are found  
                    unique_molecules.append(canonical_smiles)
                    if update:
                        start += 1
                        update = False
                    search = True 
            if not search:
                print('No new molecules found for unique molecule:', n + 1)

print(len(unique_molecules))
    
