import autode as ade
import numpy as np
ade.Config.n_cores = 4
import sys

molecule = sys.argv[1]

# Need to modify this to accept atom indices from command line 





rs_dict = {tuple([11,30]): np.linspace(1.5, 2.5,10),
	  tuple([9,28]): np.linspace(2.5,1.5,10)}
pes = ade.pes.RelaxedPESnD(species=ade.Molecule(molecule),
                           rs=rs_dict)
pes.calculate(method=ade.methods.XTB())


print("These are the saddle points")
for index in pes._saddle_points():
	print(index)

print("These are the stationary points")
for index in pes._stationary_points():
	print(index)



for i,TS in enumerate(pes.ts_guesses()):
	TS.print_xyz_file(f"{i}.xyz")

pes.save('PES_surface.npz')
