import argparse
import itertools 
from rdkit import Chem
import numpy as np
import networkx as nx 
from random import choice
from random import shuffle




# bvl_bvl smiles 

bvl_bvl = Chem.MolFromSmiles('[C@H]12C=C[C@H]3[C@@H](C(C4=C[C@H]5[C@H]6C=C[C@@H]4C=C[C@@H]56)=C2)[C@H]3C=C1')



for bond in bvl_bvl.GetBonds():
    print(bond.GetBondType())