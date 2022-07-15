import numpy as np
import argparse
from ase import Atoms
from ase.io import read
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator



# read in xyz files and create combined extended xyz file (ase format) 

parser = argparse.ArgumentParser(description='Get files to combine into .extxyz file.')
parser.add_argument("--forces_file", help="File for forces from MD trajectory",
                    type=str)
parser.add_argument("--positions_file", help="File for positions from MD trajectory",
                    type=str)
parser.add_argument("--out_file", help="Output file for combined .extxyz file", 
		   type=str)

args = parser.parse_args()




positions_data = args.positions_file
forces_data = args.forces_file

# Add to get path of positions/forces directory

if args.out_file is None and forces_data is None:
	out_filename = 'nequip_data_energy_only.extxyz'
elif args.out_file is None:	
	out_filename = 'nequip_data_energy_only.extxyz'
else:	
	out_filename = args.out_file + '.extxyz'

if positions_data is None:
	raise Exception('Error: Did not provide positions file, need atomic positions for model-training')


# Get the number of frames 
with open(positions_data, 'r') as f:
	
	num_atoms = int(f.readline())
	
	# I added the + 1 becuase the .readline() functions goes to the next line when it then counts to get num_lines
	num_lines = sum(1 for line in f) + 1
	# based on output from orca MD runs .... may generallize long term to other packages, but who needs anything besides orca + xTB :)
	# + 2 is based on the number of lines per xyz-frame, line for atom, line for energy/orca output, line for actual number of atoms 
	
	num_frames = int(num_lines/(num_atoms + 2))


Atoms =  read(positions_data, index=':')
energies_array = np.empty(shape = (num_frames,1))

#Generate energy array
index = 0
with open(positions_data,'r') as f:
        for line in f.readlines():
                if 'E_Pot' in line:
                        start = line.index('E_Pot')
                        line = line[start + 6:]
                        end = line.index(' ')
                        energies_array[index,0] = float(line[:end])
                        index += 1

if forces_data is not None:
	Forces = read(forces_data, index =':')
	forces_array = np.empty(shape = (num_frames, num_atoms, 3))
	# Generate force array
	index = 0
	for force in Forces:
        	forces_array[index, :, :] = force.get_positions()
        	index += 1


	for index in range(num_frames):
		Atoms[index].calc = SinglePointCalculator(Atoms[index], energy = energies_array[index,0], forces = forces_array[index,:,:]) 
else:
	print('Force file not specified, constructing energy only model for nequip-training')
	for index in range(num_frames):
                Atoms[index].calc = SinglePointCalculator(Atoms[index], energy = energies_array[index,0])


write(out_filename, Atoms, format='extxyz', append=True)







