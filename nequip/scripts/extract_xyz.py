import sys

file_path = sys.argv[1]

with open(file_path, 'r') as f:
	lst = f.read().splitlines()

numAtoms = int(lst[0])
totalFrames = len(lst)/(numAtoms + 2)

frame_num = 0
start = 0
stop =  numAtoms + 2



with open('MD.xyz','a') as f:
	while frame_num < totalFrames:
		f.write(lst[0] + '\n')
		properties = lst[start + 1]
		f.write(properties + '\n')
		for i in range(start + 2, stop):
			line = lst[i].split()
			f.write('  ' + line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\n')
			
			
		start = stop
		stop += numAtoms + 2
		frame_num += 1


