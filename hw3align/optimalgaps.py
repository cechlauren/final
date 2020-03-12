import sys
import numpy as np
from .sequences import *
from .getmatrix import blosum50
from .smithwaterman import sw


#This assignment asks for a TRUE POSITIVE RATE of 0.7, or 70% so:
TPR = 70

#Now, what's the best false positive rate that I can achieve?
#We'll vary gap opening from 1 to 20 and extension penalties from 1 to 5 with the BLOSUM50 matrix.

def getFpr(scorematrix, gap_start, gap_extend):
	'''
	
	This fucntion will find the false positive rate given a true positive rate of 0.7
	
	This function input is a scoring matrix and gap start/extension costs.
	This function returns the FALSE POSITIVE RATE.
	
	'''
	#Determine all true positive alignment scores to find the cutoff that sets true positive rate at 70%
	posScores = [sw(a,b, scorematrix,  gap_start, gap_extend)[2] for a,b in positives] #get scores of pos. pairs
	cutoff = np.percentile(posScores, 100-TPR) #get cutoff of those pos.scores based off TPR of 70%
	
	#Determine the false positive rate.
	negScores = [sw(a,b, scorematrix,  gap_start, gap_extend)[2] for a,b in negatives] #get scores of neg. pairs
	
	#FPR is the count of neg. scores above the cutoff we designated, divided by the number of true negatives.
	return sum(map(lambda x: x > cutoff, negScores))/len(negScores)


def optimizeGaps():
	'''
	
	This function will get the false positive rate when we have a true positive rate of 0.7 
	for the different values of gap_start and gap_extend costs with blosum50. 
	The output of the function are the results in a table that one can use to plot in R.
	
	'''
	with open("optimizeGapPenaltiesTest.txt", "w") as f:
		for  gap_start in range(-1, -21, -1): #start at -1, end at -21, step at -1 (in more neg direction , exclusive)
			for gap_extend in range(-1, -6, -1): #start at -1, end at -6, step at -1 (in more neg direction, exclusive)
				print("%s, %s:" % (gap_start, gap_extend))
				fpr = getFpr(blosum50,  gap_start, gap_extend) #using blosum50 bc asked for
				print(" ", fpr)
				f.write("\t".join(map(str, [gap_start, gap_extend, fpr])) + "\n") #write as a table
