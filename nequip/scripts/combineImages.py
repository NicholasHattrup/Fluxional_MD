import os
import sys

# Get directory of all the image directories 

filepath = sys.argv[1]

numImages = int(sys.argv[2])

for i in range(numImages):
		with open(filepath + '/' + str(i) + '/' + 'IMAGE_' + str(i) + '_forces_500K.xyz', 'r') as f:
			with open(filepath + '/' + 'combinedForces.xyz', 'a') as f_plus:
				for line in f:
					f_plus.write(line)

		with open(filepath + '/' + str(i) + '/' + 'IMAGE_' + str(i) + '_trajectory_500K.xyz', 'r') as f:
                        with open(filepath + '/' + 'combinedPositions.xyz', 'a') as f_plus:
                                for line in f:
                                        f_plus.write(line)		

