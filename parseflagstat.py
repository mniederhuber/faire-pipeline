# parseflagstat.py by cuyehara (cuyehara@email.unc.edu)
# Version 2.0

# NOTE - 10/4/16 - need to fix the directory structure to work w/ faire-pipeline
# should be able to take any directory as an argument and use this instead of the directory it's located in
# consider taking raw stat files to make it easier to work with

import os
import re
import sys


namingDict = {'qual_score':'.*q([0-9]+):',
			   'cleaned':'.*no[A-Za-z]+:',
			   'duplicates':'.*dupsMarked:'}

def getFileList(baseDir):
	pre = re.compile('(.+)_flagstats.txt')
	files = os.listdir(baseDir)
	statfiles = []
	for f in files:
		if pre.match(f):
			prefix = pre.match(f).group(1)
			filePath = baseDir + f
			statfiles.append([filePath, prefix])
	return statfiles

def makeMatrix(statfiles):
	matrix = [['Name', 'Total_Reads', 'Mapped', 'Percent_Mapped', 'MapQ_Cutoff', 'Reads_MapQ>Cutoff',
		'Percent>Cutoff', 'Duplicates', 'Percent_Duplicates','After_dup_removal', 'Mapped_to_YUHet', 'Final_Read_Count']]
	for f in statfiles:
		entry = getEntry(f)
		matrix.append([entry[col] for col in matrix[0]])
	return matrix

def getEntry(f):
	with open(f[0], 'r') as i:
		prefix = f[1]
		print (prefix)
		entry = {}
		while True:
			line = i.readline()
			if not line:
				break
			line = line.strip()
			#print (prefix)
			if re.match(prefix + ':', line):
				entry['Name'] = line[0:-1]
				entry['Total_Reads'] = i.readline().strip().split()[0]
				entry['Mapped'] = i.readline().strip().split()[0]
				print (entry['Mapped'])
			elif re.match(namingDict['qual_score'], line):
				mapq_cutoff = re.match(namingDict['qual_score'], line).group(1)
				entry['MapQ_Cutoff'] = mapq_cutoff
				entry['Reads_MapQ>Cutoff'] = i.readline().strip().split()[0]
			elif re.match(namingDict['duplicates'], line):
				line = i.readline()
				entry['Duplicates'] = i.readline().strip().split()[0]
			elif re.match(namingDict['cleaned'], line):
				entry['Final_Read_Count'] = i.readline().strip().split()[0]
			else:
				continue
	entry['Percent_Mapped'] = float(entry['Mapped'])/float(entry['Total_Reads']) * 100
	entry['Percent>Cutoff'] = float(entry['Reads_MapQ>Cutoff'])/float(entry['Mapped']) * 100
	entry['After_dup_removal'] = int(entry['Reads_MapQ>Cutoff']) - int(entry['Duplicates'])
	entry['Percent_Duplicates'] = float(entry['Duplicates'])/float(entry['Mapped']) * 100
	entry['Mapped_to_YUHet'] = int(entry['After_dup_removal']) - int(entry['Final_Read_Count'])
	return entry
		
def main():
	if len(sys.argv) == 2:
		baseDir = sys.argv[1]
	else:
		baseDir = os.getcwd()
	if baseDir[-1] != '/':
		baseDir += '/'
	statfiles = getFileList(baseDir)
	matrix = makeMatrix(statfiles)
	with open(baseDir + "collected_flag_statfiles.csv", 'w') as o:
		for m in matrix:
			line = ','.join([str(x) for x in m])
			o.write(line + '\n')

if __name__ == '__main__':
	main()
