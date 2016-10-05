# z_norm_v2.py by CMUYEHARA
# Z Normalizes wig files by chromosome arm

import os
import re
import sys
import math
import numpy as np

def mean_stdev(input_file):
	norm_dict = {}
	with open(input_file, 'r') as i:
		i.readline()
		line = i.readline()
		chrom = line.strip().split('\t')[0]
		#print (chrom)
		values = []
		while line:
			while '#' in line or line.strip().split('\t')[0] == chrom:
				if "#" in line:
					line = i.readline()
					continue
				line = line.strip().split('\t')
				dist = int(line[2]) - int(line[1])
				v = float(line[3])
				values[len(values):] = [v] * dist
				line = i.readline()
			values = np.asarray(values)
			mean = values.mean()
			stdev = values.std()
			norm_dict[chrom] = [mean, stdev]
			print ("{}\t{}\t{}\t{}".format(input_file, chrom, mean, stdev))
			#print (chrom, mean, stdev, line.strip().split('\t')[0])
			chrom = line.strip().split('\t')[0]
			values = []
			'\t'.join(line)
		return norm_dict
			
def z_score(norm_dict, input_file, output_file):
	with open(input_file, 'r') as i:
		with open(output_file, 'w') as o:
			line = i.readline()
			o.write(line)
			print (line)
			while True:
				line = i.readline()
				if not line:
					break
				elif '#' in line:
					o.write(line)
					continue
				line = line.strip().split('\t')
				chrom = line[0]
				line[3] = (float(line[3]) - norm_dict[chrom][0])/norm_dict[chrom][1]
				o.write('\t'.join([str(n) for n in line])+'\n')

def main():
	if len(sys.argv) != 3 or not re.match('.*(wig|wg)$',sys.argv[1]) or not re.match('.*(wig|wg)$', sys.argv[2]):
		print ("usage <input.wig> <output.wig>")
		print (sys.argv)
		quit()
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	norm_dict = mean_stdev(input_file)
	z_score(norm_dict, input_file, output_file)

if __name__ == "__main__":
	main()
