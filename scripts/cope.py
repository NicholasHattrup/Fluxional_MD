import xyz2graph 
import argparse
import itertools 
from rdkit import Chem
import numpy as np
import networkx as nx 



# Get all valid 6 membered carbons 
def getCombinations(mol_graph, atoms = None):

	# If no argument given, getCombinations default to only checking carbons 
	if atoms is None:
		atoms = 'C'
	atoms = mol_graph.nodes(data=True)
	nodes = [n for n in range(len(atoms)) if atoms[n]['symbol'] == 'C']
	sub_mol = nx.subgraph(mol_graph, nodes)
	combinations = [combination for combination in itertools.combinations(sub_mol.nodes, 6)]
	for combination in combinations:
		if not nx.is_connected(nx.subgraph(sub_mol, combination)):
			combinations.remove(combination)
	return combinations



def countpi(mol_graph, nodes):
	pi = 0 
	edges = mol_graph.edges(nbunch=nodes,data=True)
	for edge in edges:
		if edge[2]['bondtype'] == str(Chem.BondType.DOUBLE):
			pi +=1
	return pi 

def nodes2edges(mol_graph, nodes):
	sub_mol = nx.subgraph(mol_graph, nodes)
	return sub_mol.edges


def get_hamiltonian(graph):
	



if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage='%(prog)s [options] molecule.xyz')
    parser.add_argument('structure', metavar='structure', type=str)
    parser.add_argument('--no-charged-fragments',
        action="store_true",
        help="Allow radicals to be made")
    parser.add_argument('--no-graph',
        action="store_true",
        help="Run xyz2mol without networkx dependencies")

    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    parser.add_argument('--use-huckel',
        action="store_true",
        help="Use Huckel method for atom connectivity")
    parser.add_argument('-o', '--output-format',
        action="store",
        type=str,
        help="Output format [smiles,sdf] (default=sdf)")
    parser.add_argument('-c', '--charge',
        action="store",
        metavar="int",
        type=int,
        help="Total charge of the system")
    
    parser.add_argument('--aromaticity',
        action="store_true",
        help="Detect and set aromatic bonds")

    args = parser.parse_args()

    # read xyz file
    filename = args.structure

    # allow for charged fragments, alternatively radicals are made
    charged_fragments = not args.no_charged_fragments

    # quick is faster for large systems but requires networkx
    # if you don't want to install networkx set quick=False and
    # uncomment 'import networkx as nx' at the top of the file
    quick = not args.no_graph


    detect_aromaticity = not args.aromaticity

    # read atoms and coordinates. Try to find the charge
    atoms, charge, xyz_coordinates = xyz2graph.read_xyz_file(filename)
    # huckel uses extended Huckel bond orders to locate bonds (requires RDKit 2019.9.1 or later)
    # otherwise van der Waals radii are used
    use_huckel = args.use_huckel

    # if explicit charge from args, set it
    if args.charge is not None:
        charge = int(args.charge)


    # Get atom connectivity using Huckle Extended Connectivity Method 
    AC, _ = xyz2graph.xyz2AC(atoms, xyz_coordinates, charge, use_huckel=use_huckel)

    # convert AC matrix to bond order (BO) matrix
    BO, atomic_valence_electrons = xyz2graph.AC2BO(
        AC,
        atoms,
        charge,
        allow_charged_fragments=charged_fragments,
        use_graph=quick)

    # Generate molecular graph
    mol_graph = xyz2graph.BO2graph(BO, atoms, atomic_valence_electrons, charge, allow_charged_framents=charged_fragments, use_atoms_maps=False)
    connected_combins = getCombinations(mol_graph)
    if not connected_combins:
        print('No valid 6 membered structures that are fully connected ... stopping.')
        exit()

    pi_combins = [combin for combin in connected_combins if countpi(mol_graph, combin) > 1]
    if not pi_combins:
        print('No valid 6 membered structures with at least two pi bonds ... stopping.')
    simple_combins = [combin for combin in pi_combins if  xyz2graph.edges2path(nodes2edges(mol_graph, combin)) is not None]
    print(len(simple_combins))
    

    
    