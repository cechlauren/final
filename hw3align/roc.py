import sys
from .getseq import *
from .sequences import *
from .smithwaterman import sw
import matplotlib.pyplot as plt
import sklearn.metrics as skm

#Here we make our ROC plots based off the true positive and true negative data. I don't know how to make
#these without incorporating another program so I've implemented sklearn here. 

#ROC plots 

def makeRocPlotLC(scorematrix, gap_start, gap_extend, name):
	'''
    
	This function input is a scoring matrix and gap start/extend params.
    This function will find the ROC using known TP and TN data.
     
    
	'''
	#Score the true positives and true negatives...we've seen this before
	posScores = [sw(a,b, scorematrix, gap_start, gap_extend)[2] for a,b in positives]
	negScores = [sw(a,b, scorematrix, gap_start, gap_extend)[2] for a,b in negatives]
	#Sklearn's ROC stuff expects a single vector of truepositive/truenegative identities
	# and a single vector of scores. so need to make those vectors here:
	allScores = posScores + negScores #single vector of scores
	allClasses = [1]*len(posScores) + [0]*len(negScores) #single vector of tp and tn identities
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores) #now put them into sklearn program
	roc_auc = skm.auc(fpr, tpr) #the roc function takes tpr and false positive rate
    
	#plot. Adapted from https://stackoverflow.com/questions/56203889/how-to-get-the-optimal-threshold-from-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Curve for' + name + 'with gap start %s and extend %s' % \
	 (gap_start, gap_extend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.show()
	

    
#May be prudent to consider how sequence length influences the scores and therefore the tpr and fpr.
#Try to normalize for the string length, based off of question 3 in part 1
def makeRocPlot_normScoresLC(scorematrix, gap_start, gap_extend, name):
	'''
    
	This function input is a scoring matrix and gap start/extend params.
    This function will find the ROC using TP and TN data, scaling by smaller string length.
    This function output is a PDF plot. 
    
	'''
	#As above, score the TP and TN
	posScores = [sw(a,b, scorematrix, gap_start, gap_extend)[2]/min(len(a), len(b)) for a,b in positives]
	negScores = [sw(a,b, scorematrix, gap_start, gap_extend)[2]/min(len(a), len(b)) for a,b in negatives]
	#Sklearn's ROC stuff expects a single vector of truepositive/truenegative identities
	# and a single vector of scores. so need to make those vectors here:
	allScores = posScores + negScores #single vector of scores
	allClasses = [1]*len(posScores) + [0]*len(negScores) #single vector of tp and tn identities
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores) #plug into sklearn program
	roc_auc = skm.auc(fpr, tpr) #get roc output
    
	#plot it
    #Adapted from https://stackoverflow.com/questions/56203889/how-to-get-the-optimal-threshold-from-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Curve for' + name+ 'with gap start %s and extend %s\nwith normalized scores' % \
	 (gap_start, gap_extend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.show()
	
    

#in the instance that we provide a given set of scores for the true positives and true negatives 
def makeRocPlot_givenScores(posScores, negScores, gap_start, gap_extend, name):
	'''
    
	This function input is a list of TP and TP scores.
    This function output is a ROC and PDF plots.
    
	'''
	#Sklearn's ROC stuff expects a single vector of truepositive/truenegative identities
	# and a single vector of scores. so need to make those vectors here:
	allScores = posScores + negScores #single vector of scores
	allClasses = [1]*len(posScores) + [0]*len(negScores)
	fpr, tpr, threshold = skm.roc_curve(allClasses, allScores)
	roc_auc = skm.auc(fpr, tpr)
    
	#plot it
    #Adapted from https://stackoverflow.com/questions/56203889/how-to-get-the-optimal-threshold-from-roc-curve-in-python
	fig, axes = plt.figure(), plt.axes()
	plt.title('Receiver Operating Curve for %s\nwith gap start %s and extend %s' % \
	 (name, gap_start, gap_extend))
	plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	axes.set_aspect('equal', 'box')
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	#plt.show()
	fig.savefig("ROC_normScores_%s.pdf" % name)
