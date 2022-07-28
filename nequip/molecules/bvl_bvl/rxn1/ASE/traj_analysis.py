import os
import argparse
from numpy import arange
from numpy.random import choice
from numpy.linalg import norm
import matplotlib.pyplot as plt
plt.style.use('ipynb')
from ase.io.trajectory import TrajectoryReader
from ase.io.xyz import write_xyz
parser=argparse.ArgumentParser()
parser.add_argument('--samples', type=int, default=None, help='Number of trajectory samples to plot')
parser.add_argument('--base', type=str, help='Basename of trajectory files') 
# Please note that this expects 0 indexed pairs!i.e. if in avogrado/chimeraX the atoms are (17, 18) then specify the pair as (16, 17)
parser.add_argument('--pair_one', nargs="+", type=int, help='Atom pair for first reacting bond')
parser.add_argument('--pair_two', nargs="+", type=int, help='Atom pair for second reacting bond')
args=parser.parse_args()


# Will only check current directory
trajectories = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f)) and args.base in f]


# Randomly select k trajectories from n total trajectories for viewing/analysis 

indices=arange(0,len(trajectories))
if not args.samples is None:
	indices = choice(a=indices, size=args.samples, replace=False) 



atoms_one=args.pair_one
atoms_two=args.pair_two
pair_one_all=[]
pair_two_all=[]
for index in indices:
	trajectory=TrajectoryReader(filename=trajectories[index])
	pair_one=[]
	pair_two=[]
	for frame in trajectory:
		dist=norm(frame.positions[atoms_one[0]]-frame.positions[atoms_one[1]])
		pair_one.append(dist)
		dist=norm(frame.positions[atoms_two[0]]-frame.positions[atoms_two[1]])
		pair_two.append(dist)
	pair_one_all.append(pair_one)
	pair_two_all.append(pair_two)


def check_dist(traj):
	for dist in traj:
		if dist > 3:
			return True	
	return False

#print(pair_one_all)
#print(pair_two_all)

steps = arange(len(pair_one_all[0]))
for k, (pair_one, pair_two) in enumerate(zip(pair_one_all, pair_two_all)):
	if check_dist(pair_one):
		k_to_write = k
		print(k)
		continue
	plt.plot(steps, pair_two, label = 'Pair_Two')
	plt.ylim([1.25, 2.75])
	plt.title('C-C atom distance vs. simulation time')
	plt.xlabel('Step (1 step = 5 fs)')
	plt.ylabel('C-C atom distance A')
	#plt.plot(steps, pair_two, label = 'Pair_Two')
plt.show()


#trajectory = TrajectoryReader(filename=trajectories[k_to_write])
#for atoms in trajectory:
#	atoms.write('Bad_Trajectory.xyz',format='xyz',append=True)


			
	
	




