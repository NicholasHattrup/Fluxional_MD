import argparse
import os
from ase import io
import numpy as np
import matplotlib.pyplot as plt
parser= argparse.ArgumentParser()
parser.add_argument('--directory', type=str, help='Directory with associated trajectory files')

args = parser.parse_args()


directory = os.fsencode(args.directory)
    

pair_distance_1 = []
pair_distance_2 = []
for file in os.listdir(directory):
    print(file)
    pair_distance_1_sample = []
    pair_distance_2_sample = []
    filename = os.fsdecode(file)
    xyz = io.read(args.directory + '/' + filename, index=':')
    for frame in xyz:
        positions = frame.positions
        C_12 = positions[11,:]
        C_20 = positions[19,:]
        C_17 = positions[16,:]
        C_18 = positions[17,:]
        pair_distance_1_sample.append(np.linalg.norm(C_20 - C_12)) 
        pair_distance_2_sample.append(np.linalg.norm(C_18 - C_17)) 
        
    pair_distance_1.append(pair_distance_1_sample)
    pair_distance_2.append(pair_distance_2_sample)
time = np.arange(0, 2000, 10)

print(pair_distance_1[0])


print(pair_distance_1[1])


exit()

print(pair_distance_1)




plt.figure()


plt.plot(time, pair_distance_1[0], label = 'pair_1', color = 'black')
plt.plot(time, pair_distance_2[0], label = 'pair_2', color = 'blue')
plt.plot(time, pair_distance_1[31], label = 'pair_1', color = 'black')
plt.plot(time, pair_distance_2[31], label = 'pair_2', color = 'blue')
    

plt.show()