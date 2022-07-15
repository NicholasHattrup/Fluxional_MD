import autode as ade
import numpy as np
import networkx as nx
import argparse
import itertools

''' 
	Workflow for molecular processing:
1. Generate molecular graph from reliable xyz data
2. Remove all hydrogen atoms to get condensed molecular graph
3. Generate all possible 6 heavy atom combinations from the remaining nodes in the molecular graph

	Filtering section
4. Check if the subgraphs formed from these combinations are connected, if not remove
5. Check if first pass subgraphs contain at least two pi bonds, if not remove
6. Of the remaining combinations, find carbons at end of "line" and ensure both are sp2 or sp at least 
7. Finally, ensure that at least two pi bonds are 
'''

parser = argparse.ArgumentParser(description='Detect possible cope arrangments in a given fluxional system')
parser.add_argument('--file', type=str,
                    help='file with associated xyz coordinates')
parser.add_argument('--smiles', type=str, 
					help='associated molecular smiles string represenation' )
args = parser.parse_args()


def subgraph(mol_graph, atoms):
	mol_sub = mol_graph.subgraph(atoms)
	if nx.is_connected(mol_sub):
		return mol_sub
	else:
		return None

# Assumes edge data format from molecular graph generated from autode 
def countpi(mol_graph):
	num_pi_bonds = 0
	for edge in mol_graph.edges:
		if mol_graph.get_edge_data(edge[0],edge[1])['pi']:
			num_pi_bonds += 1
	return num_pi_bonds

def is_pi(mol_graph, atoms):
	for atom in atoms: 
		bonds = mol_graph.edges(atom)
		for bond in bonds:
			if not mol_graph.get_edge_data(bond[0], bond[1])['pi']:
				return False  
		return True 

def get_pi_atoms(mol_graph):
	pi_atoms = [] 
	for edge in mol_graph.edges:
		if mol_graph.get_edge_data(edge[0], edge[1])['pi']:
			if not edge[0] in pi_atoms:
				pi_atoms.append(edge[0])
			if not edge[1] in pi_atoms:
				pi_atoms.append(edge[1])
	return pi_atoms



# In house check for carbon valency issues encountered with
# the moleculer graph generated from autodE
def check_pi_assignemnt(mol_graph):
	atoms_to_check = []
	pi_atoms = get_pi_atoms(mol_graph)
	for atom in pi_atoms:
		valency = 0
		for edge in mol_graph.edges(atom):
			if mol_graph.get_edge_data(edge[0], edge[1])['pi']:
				valency += 2
			else:
				valency += 1
		if valency > 4:
			atoms_to_check.append(atom)
	return atoms_to_check



# Code to reassign bond attributes for carbon atoms flagged
# with valency issues (i.e. bond order > 4)
def reassign_valency(mol_graph, atoms):
	


# For a given graph, find if there exists a simple path
# For now only developed for acylic graphs .... will need to augment later on  
# I will need to rewrite for method using used and unused nodes 




def simple_path(mol_graph, is_cyclic = False):
	# Get atoms and bonds in heavy atom molecular graph 
	bonds = list(mol_graph.edges)
	atoms = list(mol_graph.nodes)
	# Simple path to be constructed
	first_bond = bonds[0]
	bond_path = [first_bond]
	bonds.remove(first_bond)

	# Get start edges in sequence and remove connecting atom from "active" atom list
	for bond in bonds:
		atom1, atom2 = bond
		if atom1 in first_bond:
			second_bond = (atom1, atom2)
			bond_path.insert(len(bond_path), second_bond)
			bonds.remove(bond)
			atoms.remove(atom1)
			break
		elif atom2 in first_bond:
			second_bond = (atom1, atom2)
			bond_path.insert(len(bond_path), second_bond)
			bonds.remove(bond)
			atoms.remove(atom2)
			break
	
	# So now we have a first and second bond that we can use to continue to build up our sequence
	# Can only add if atom is in active atoms
	build = True 
	while build: 
		build = False 
		for bond in bonds:
			atom1, atom2 = bond 
			if atom1 in first_bond and atom1 in atoms:
				first_bond = (atom1, atom2)
				bond_path.insert(0, first_bond)
				bonds.remove(bond)
				atoms.remove(atom1)
				build = True 
				break
			elif atom2 in first_bond and atom2 in atoms:
				first_bond = (atom1, atom2)
				bond_path.insert(0, first_bond)
				bonds.remove(bond)
				atoms.remove(atom2)
				build = True
				break
			elif atom1 in second_bond and atom1 in atoms:
				second_bond = (atom1, atom2)
				bond_path.insert(len(bond_path), second_bond)
				bonds.remove(bond)
				atoms.remove(atom1)
				build = True
				break
			elif atom2 in second_bond and atom2 in atoms:
				second_bond = (atom1, atom2)
				bond_path.insert(len(bond_path), second_bond)
				bonds.remove(bond)
				atoms.remove(atom2)
				build = True
				break
			else:
				# Cannot keep generating sequence from any bonds in the sequence .... simple path does NOT exist
				return None 
		if not is_cyclic and len(atoms) == 2:
			return atoms 
		# Need to add condition for cylic graph

				
		

	
 
# Generate autode Molecule object and associated molecular graph
mol = ade.Molecule(args.file)
#mol = ade.Molecule(smiles=args.smiles)
mol_graph = mol.graph
nodes = mol_graph.nodes

edges = mol_graph.edges 


# Need to check initial graph construction before removing hydrogens 



check_pi_assignemnt(mol_graph)

exit()

# Remove all hydrogens
# Need to decide if I want to remove all non-carbon atoms 
for index in range(len(nodes)):
	if nodes[index]['atom_label'] == 'H':
		mol_graph.remove_node(index)


# Cope rearrangment takes place between five carbon atoms (possibly 5 heavy atoms?)
# Generate all 6 member combinations of carbon atoms in the molecule
carbon_combins = [combin for combin in itertools.combinations(mol_graph.nodes, 6)]
filtered_combins = []
# 4 step filtering of atom combinations for cope rearrangment 
for combin in carbon_combins:
	break
	mol_sub = subgraph(mol_graph, combin) 
	if mol_sub is None: # Is the subgraph connected ?
		continue
	if 2 > countpi(mol_sub): # Does the subgraph have *atleast* two pi bonds 
		continue 
	end_atoms = simple_path(mol_sub) # Does the group of atoms contain  a simple path
	if end_atoms is None:
		continue 
	if not is_pi(mol_sub , end_atoms):
		continue
	filtered_combins.append(combin)


edges = mol_graph.edges(9) 

for edge in edges:
	print(mol_graph.get_edge_data(edge[0], edge[1]))


for combin in filtered_combins:
	break
	#mol_sub = subgraph(mol_graph, combin) 
	#end_atoms = simple_path(mol_sub)
	print(combin)


