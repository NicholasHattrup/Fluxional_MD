import autode as ade
import numpy as np


pes = ade.pes.RelaxedPESnD.from_file('PES_surface.npz')

stationary_point_one = pes._coordinates[3,6]
stationary_point_two = pes._coordinates[4,5]



np.savetxt('stationary_point_one.xyz', stationary_point_one)



