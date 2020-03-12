import sys

def read_substitution_matrix(filename):
	'''
	This function will be able to read in an amino acid substitution matrix from a file
	that we need to apply to our comparisons of sequences.
	
	This function outputs a symmetrical sub. mat. dictionary containing the value of the substitution score.
	It is keyed by a frozenset of the substitutions, where keys for matches will be of length 1. 
	'''
	matrix = {} #make a dictionary
	readHeader = False 
	header = [] #an array
	rowi = 0 #start at zero
	for line in open("hw3align/" + filename, "r"):
		#We dont care about comments
		if line[0] == "#":
			continue
		#Read in the line and make sure to separate on any whitespace
		fields = line.rstrip().split()
		#Read in the header
		if readHeader == False:
			header = fields
			colLabels = fields #column labels
			readHeader = True #now its been read
			continue
		#Confirm the correct number of fields
		assert len(fields) == len(header)
		row = header[rowi]
		#Move through the fields and add each matrix entry
		for col, score in zip(header, fields):
			key = frozenset((row, col))
			if key not in matrix:
				matrix[key] = int(score) #make sure we report integers
			else:
				assert matrix[key] == int(score)
		#Move on to next row
		rowi += 1
	return matrix

blosum50 = read_substitution_matrix("BLOSUM50")
blosum62 = read_substitution_matrix("BLOSUM62")
matio = read_substitution_matrix("MATIO")
pam100 = read_substitution_matrix("PAM100")
pam250 = read_substitution_matrix("PAM250")
#optblosum50 = read_substitution_matrix("optimalblosum.txt")

matrix_names = ["BLOSUM50", "BLOSUM62", "MATIO", "PAM100", "PAM250"]
matrices = map(read_substitution_matrix, matrix_names)
