import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--base', help='Basename for files to append', type=str)
parser.add_argument('--num', help='Number of files to append', type = int)
args = parser.parse_args()

if args.base is None:
	raise Exception('Need a base name to append with!')
if args.num is None:
	raise Exception('I know this is annoying ... but for more user control please specify the number of files to specify')

base = args.base
N = args.num

# I assume files are 0 indexed (i.e. 9 files to append corresponds to base_0 .... base_8)

with open(base + '_combined.xyz', 'a') as file_large:
	for n in range(N):
		with open(base + '_' + str(n) + '.xyz', 'r') as file_small:
			lst = file_small.read().splitlines()
		for line in lst:
			file_large.write(line + '\n')
		

