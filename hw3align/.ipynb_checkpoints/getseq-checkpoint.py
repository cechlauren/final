import sys

def readFasta(filename):
	'''
	This function will read in a sequence from a fasta file
	and return a string.
	'''
	seq = ""
	for line in open(filename, "r"):
		if line[0] == ">":
			continue
		seq = seq + line.rstrip().upper()
	return seq

#Still need to get the sequence names so I have a separate file containing the location of the fasta files
#as well as their names:
seqnames = []
for line in open("hw3align/seqnamesprep.txt", "r"):
	seqnames.append(line.rstrip())

#This output will take the sequence name that I grabbed above and attach the appropriate sequence string.
sequences = {}
for filename in seqnames:
	sequences[filename] = readFasta("hw3align/"+filename)
#This identifies the positive pairs (given in Pospairs.txt) and will make an array containing the positive pairs tuples
positivesfiles = [line.rstrip().split() for line in open("hw3align/Pospairs.txt")]
positives = [(sequences[fileA], sequences[fileB]) for fileA, fileB in positivesfiles]

#As above but with negative pairs
negativesfiles = [line.rstrip().split() for line in open("hw3align/Negpairs.txt")]
negatives = [(sequences[fileA], sequences[fileB]) for fileA, fileB in negativesfiles]
