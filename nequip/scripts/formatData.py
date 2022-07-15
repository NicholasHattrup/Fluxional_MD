import numpy as np

from ase import Atoms
from ase.io import read
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator



# read in npz file


in_filename = 'bvl_bvl_TS_trajectory_500K.xyz'
in2_filename = 'bvl_bvl_TS_forces_500K.xyz'
out_filename = 'nequip-data.extxyz'



num_frames = 413
num_atoms = 38



Atoms =  read(in_filename, index=':')
Forces = read(in2_filename, index =':')

forces_array = np.empty(shape = (num_frames, num_atoms, 3))
energies_array = np.empty(shape = (num_frames,1))




#Generate energy array
index = 0  
with open(in_filename,'r') as f:
	for line in f.readlines():
		if 'E_Pot' in line:
			start = line.index('E_Pot')
			line = line[start + 6:]
			end = line.index(' ')
			energies_array[index,0] = float(line[:end])
			index += 1


# Generate force array
index = 0
for force in Forces:
	forces_array[index, :, :] = force.get_positions()
	index += 1



for index in range(num_frames):
	Atoms[index].calc = SinglePointCalculator(Atoms[index], energy = energies_array[index,0], forces = forces_array[index,:,:]) 



write(out_filename, Atoms, format='extxyz', append=True)







