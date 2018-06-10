#!/usr/bin/env python
from sys import argv
import sys

# definite command ption
#opts, args = getopt.getopt(sys.argv[1:], 'hm:s:a:o:h:')
print(argv)

file = argv[1]

output_file = argv[2]
#input: -m modification file(ID_FT.txt) -s sequence file(uniprot_human_reviewed.fasta) -a C -o ActiveCys
#for op, value in opts:
#	if op == '-f':
#		file = value
#	elif op == '-o':
#		output_file = value
#	elif op == '-h':
#		sys.exit()

# define a class to store the ms1 spectra
class MS1:
	def __init__(self, name, RT, ions_mz_list,ions_int_list):
		self.name = name
		self.RT = RT
		self.ions_mz_list = ions_mz_list
		self.ions_int_list = ions_int_list

# read the ms1 file
def read_ms1(file):
	name = ''
	items = []
	# set the initial index as 0, so that it will not append the Fasta class items when meeting the first ">"
	index = 0
	for line in file:
		# if the line starts with "S", and index is greater than 1, set the corresponding name,RT and spectra as a grouped "MS1" class
		if line.startswith("S"):
			# if the index is 0, we don't have name and seq pair, so just read in the name but not add anything the items
			if index >= 1:
				aninstance = MS1(name, RT, ions_mz_list,ions_int_list)
				items.append(aninstance)
			# every time we add a name and seq pair to the items, reset the seq as ''
			RT = ''
			ions_mz_list = ''
			ions_int_list = ''
			index += 1
			# name is the whole line
			name = line.strip().split('\t')[1]
			if float(name) % 2000==0:
				print(name)
		elif line[2:9] == 'RetTime':
			RT = float(line.strip().split('\t')[2])
		# the lines not starting with ">" are ion lines
		elif line[1].isdigit():
			line_split = line.strip().split(' ')
			ion_mz = ''.join([line_split[0],';'])
			ion_int = ''.join([line_split[1],';'])
			ions_mz_list += ion_mz
			ions_int_list += ion_int
	#add the final sequence
	aninstance = MS1(name, RT, ions_mz_list,ions_int_list)
	items.append(aninstance)
	return items



# use the ID_FT.txt as the input, which is created by extracting the ID and FT line from the uniprot full imformation db



ms1_file = open(file, 'r').readlines()
ms1_read = read_ms1(ms1_file)

output = open (output_file, 'w')
output.write("Scan_No")
output.write(',')
output.write("RT")
output.write(',')
output.write("ions_mz_list")
output.write(',')
output.write("ions_int_list")
output.write('\n')
for i in ms1_read:
	output.write(i.name)
	output.write(',')
	output.write(str(i.RT))
	output.write(',')
	output.write(str(i.ions_mz_list))
	output.write(',')
	output.write(str(i.ions_int_list))
	output.write('\n')
output.close()





































