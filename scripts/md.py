from nequip.ase import NequIPCalculator
from ase import units
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import Stationary, ZeroRotation
from ase.io.trajectory import Trajectory as ASETrajectory
import argparse
import torch
import time
parser=argparse.ArgumentParser()
parser.add_argument('--xyz', type=str,help='Path to xyz to construct ASE Atoms class from')
parser.add_argument('--model',type=str,help='Path to nequip deployed model to run dynamics with')
parser.add_argument('--temperature', type=float,help='Temperature to run dynamics at')
parser.add_argument('--dt', type=float,help='Timestep for dynamics')
parser.add_argument('--energy_units_to_eV', type=float, default=27.2114, help='model units -> eV factor')
parser.add_argument('--length_units_to_A', type=float, default=1, help='model units -> A')
parser.add_argument('--num_steps', type=float, default=1000, help='Number of steps to run')
parser.add_argument('--interval', type=float, default=10, help='Write current trajectory frame every n steps')
parser.add_argument('--samples', type=int, default=1, help='Number of trajectories to sample/run')
parser.add_argument('--dir',type=str,default='./',help='Path to write trajectory files to')
#Define system
args = parser.parse_args()

atoms = read(args.xyz, index=0)
init_xyz = atoms.positions.copy()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

atoms.calc = NequIPCalculator.from_deployed_model(
    model_path=args.model,
    device=device,
    species_to_type_name = {
        "H": "H",
        "B": "B",
	"C": "C",
        "N": "N",
	"O": "O",
        "F": "F",
        "S": "S",
        "Cl": "Cl",
        "Br": "Br",
        "I": "I"
    },
    energy_units_to_eV=args.energy_units_to_eV,
    length_units_to_A=args.length_units_to_A,
)

'''
nvt_dyn = Langevin(
	atoms=atoms,
	temperature_K=args.temperature,
	timestep=args.dt * units.fs,
	friction=0.02)
'''
print('Running Trajectories using ' + str(device))
start_time = time.time()
for i in range(91,args.samples):
	nvt_dyn = Langevin(
        atoms=atoms,
        temperature_K=args.temperature,
        timestep=args.dt * units.fs,
        friction=0.02)
	traj_file = args.dir + '/' + 'Trajectory_' + str(i) + '.traj'
	print(i, traj_file)
	MaxwellBoltzmannDistribution(atoms=atoms, temp=args.temperature * units.kB)
	ZeroRotation(atoms) # Set center of mass momentum to zero
	Stationary(atoms) # Set rotation about center of mass zero
	traj = ASETrajectory(traj_file, 'w', atoms)
	traj.write(atoms)
	nvt_dyn.attach(traj.write, interval=args.interval)
	nvt_dyn.run(steps=args.num_steps)
	traj.close()
	# reset atom positions to initial sampling geometry
	atoms.set_positions(init_xyz.copy())

print("--- %s seconds ---" % (time.time() - start_time))
