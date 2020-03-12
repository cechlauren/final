import sys

import numpy as np




def sw(a, b, substitutionmatrix, startcost, extendcost):

	'''

	This function uses Smith-Waterman to find a maximal local alignment between two 
	sequences.
	We'll have a big table that has indices keeping track of the aligned 
	seq's beginning and ends. An ideal solution will only need to compute
	each sub-alignment once, while providing coverage of the whole
	search space.


	The input for this function is two strings (a,b) to align,  a

	substitution matrix, a gap start cost, and finally a gap extension cost.

	Most commonly in evolution, gaps will occur in bunches, so it would make 
	sense that once a break in alignment is made that it would be easier to 
	make multiple insertions or deletions. In contrast, to introduce the first 
	gap (the "start"/"open" gap) would be more costly since a break must occur
	in the DNA. 
	Therefore, a gap cost is: (d + (n-1)*e), where "n" is the length of

	the gap, "d" is the gap opening, and "e" is the gap extend.
	By this logic, the gap extension cost will only be only bepaid for 2+ 
	length gaps.



	This function will return:
	1) a tuple of the indices for beginning of local alignment in a and b 
	2) a tuple of the indices of the end of the alignment 
	3) the score of the aligning section 
	4) the aligned sections of a and b (gaps = "-").

	'''

	# initialize matrices and the trace matrices
	#Local alignment initializes with top left =0

	match = np.zeros(shape=(len(a) + 1, len(b) + 1)) #size of best score matrix according to size of seq.

	traceMatch = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)] #initialize trace

	Y = np.zeros(shape=(len(a) + 1, len(b) + 1)) #initialize matrix to keep track of set "Y" matches

	traceY = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)] #trace

	X = np.zeros(shape=(len(a) + 1, len(b) + 1)) #initialize matrix to keep track of "X" matches

	traceX = [["" for j in range(len(b)+ 1)] for i in range(len(a)+ 1)] #trace

	# we fill the first rows/columns of matrices to avoid alignments causing unusual edges

	for i in range(len(a) + 1):

		Y[i, 0] = float("-inf") #basically we dont want to match things that have no length

		if i>0: 

			X[i, 0] = startcost + extendcost*(i-1) 
			#if the length is greater than zero then we can start doing matches or gaps

	for j in range(len(b) + 1):

		X[0, j] = float("-inf") #basically we dont want to match things that have no length

		if j>0:

			Y[0, j] = startcost + extendcost*(j-1)
			#if the length is greater than zero then we can start doing matches or gaps

	#Fill in matrices by comparing each element in our sequences

	maxSeen = (float("-inf"), (-1, -1)) # this is (score, (i,j))

	for i in range(1, len(a) + 1): #go thru seq a

		for j in range(1, len(b) + 1):  #go thru seq b

			extendX = X[i-1, j] + extendcost  # extend gap in b in matrix X

			openX = match[i-1, j] + startcost # open gap in b in matrix X

			if openX > extendX:

				X[i,j] = openX #if start cost more than extend cost then score is start cost

				traceX[i][j] = "openX" #trace starting from start cost

			else:

				X[i,j] = extendX

				traceX[i][j] = "extendX"

			extendY = Y[i, j-1] + extendcost  # extend gap in a in matrix X

			openY = match[i, j-1] + startcost # open gap in a in matrix X

			if openY > extendY: 

				Y[i,j] = openY #same idea as above

				traceY[i][j] = "openY"

			else:

				Y[i,j] = extendY

				traceY[i][j] = "extendY"

			match[i,j] = match[i-1, j-1] + substitutionmatrix[frozenset((a[i-1], b[j-1]))]
			#sub.mat. contains an immutable set initialized with the given iterables
			#the score for the comparison of amino acid in seq a or seq b depends on those sub. mat. scores and
			#on the fact if they are indeed a match or not.

			traceMatch[i][j] = "match" #a match for tracing matrix

			if Y[i,j] > match[i,j]:

				match[i,j] = Y[i,j]

				traceMatch[i][j] = "closeY" #there are no more common subseq

			if X[i,j] > match[i,j]:

				match[i,j] = X[i,j]

				traceMatch[i][j] = "closeX" #there are no more common subseq

			if match[i,j] > maxSeen[0]:

				maxSeen = (match[i,j], (i, j)) #record the highest matching scores!

	#Start the traceback based on the best score between seq comparisons and return that subsequence

	i,j = maxSeen[1] #start at the highest score

	end = (i,j) 

	matrix, traceMatrix = match, traceMatch

	score = matrix[i,j] #this matrix told us the best scores available

	bestScore = score #we only care about the best score

	matchedA = "" #the string of matched seq a 

	matchedB = ""

	while score > 0:

		if traceMatrix[i][j] == "match":

			matchedA = a[i-1] + matchedA #start building the string based on matches

			matchedB = b[j-1] + matchedB

			i,j = i-1, j-1

		elif traceMatrix[i][j] == "closeX":

			matrix, traceMatrix = X, traceX #the end of a subseq recorded in matrix X

		elif traceMatrix[i][j] == "closeY":

			matrix, traceMatrix = Y, traceY

		elif traceMatrix[i][j] == "openY":

			matchedA = "-" + matchedA #introduction of a gap

			matchedB = b[j-1] + matchedB #we did not have a match

			i,j = i, j-1

			matrix, traceMatrix = match, traceMatch

		elif traceMatrix[i][j] == "extendY":

			matchedA = "-" + matchedA #introduction of the extension gaps

			matchedB = b[j-1] + matchedB #we still have no match

			i,j = i, j-1

		elif traceMatrix[i][j] == "openX":

			matchedA = a[i-1] + matchedA #we had no match

			matchedB = "-" + matchedB #introduction of start gap

			i,j = i-1, j

			matrix, traceMatrix = match, traceMatch

		elif traceMatrix[i][j] == "extendX":

			matchedA = a[i-1] + matchedA #we still have no match

			matchedB = "-" + matchedB #introduction of extension gaps

			i,j = i-1, j

		else: print("hmm something is wrong here")

		score = matrix[i,j] #ze final score!

	start = (i, j)

	#print(start, end, bestScore)

	#print(matchedA)

	#print(matchedB)

	return start, end, bestScore, matchedA, matchedB




