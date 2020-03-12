import sys
import optimalgaps
from optimalmatrix import getAlignments, optimizeMatrix_geneticAlg
from getmatrix import matrices, matrixNames, blosum50, matio
import roc


# will get opt gap penalties to minimize FPR when TPR is 0.7
if sys.argv[1][0:2] == '-A':
	optimizeGaps.optimizeGaps()

# compares scoring matrices with ROC curves
if sys.argv[1][0:2] == '-B':
	gapStart = -8
	gapExtend = -3
	for name, matrix in zip(matrixNames, matrices):
		roc.makeRocPlotLC(matrix, gapStart, gapExtend, name)

#will make min seq normalized plot
if sys.argv[1][0:2] == '-C':
	roc.makeRocPlot_normScoresLC(blosum50, gapStart, gapExtend, "blosum50")

# optimizes the blosum50, first gets alignment 
if sys.argv[1][0:2] == '-D':
	gap_start = -8
	gap_extend = -3
	scorematrix = blosum50
	true_pos_align, true_neg_align = getAlignments(scorematrix=blosum50,
		gap_start=-8, gap_extend=-3)
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		scorematrix=blosum50, mutationChance=0.5, mutationAmount=0.1,
		selectivePressure=1, N=100, totalStepsToStop=1000,
		stepsWithNoImprovement=100, librarySize=10, gap_start=-7,
		gap_extend=-4, true_pos_align, true_neg_align)
	bestMat_blosum50 = library[max(library.keys())]
	#ROC curve with alignments
	negScores = [scoreAlignment(a,b, bestMat_blosum50, gap_start, gap_extend) for a,b in true_neg_align]
	posScores = [scoreAlignment(a,b, bestMat_blosum50, gap_start, gap_extend) for a,b in true_pos_align]
	roc.makeRocPlot_givenScores(posScores, negScores, gap_start, gap_extend, "Optblosum50")
	#do one with rescoring
	roc.makeRocPlotLC(bestMat_blosum50, gap_start, gap_extend, "Optblosum50_rescored")

# optimize matio
	# get alignments
if sys.argv[1][0:2] == '-E':
	gap_start = -8
	gap_extend = -3
	scorematrix = blosum50 #since these are best alignments 
	true_pos_align, true_neg_align = getAlignments(scorematrix=blosum50,
		gap_start=-8, gap_extend=-3)
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		scorematrix=matio, mutationChance=0.5, mutationAmount=0.1,
		selectivePressure=1, N=100, totalStepsToStop=1000,
		stepsWithNoImprovement=100, librarySize=10, gap_start=-7,
		gap_extend=-4, true_pos_align, true_neg_align)
	optmatio = library[max(library.keys())]
	#ROC curve with alignments
	negScores = [scoreAlignment(a,b, optmatio, gap_start, gap_extend) for a,b in true_neg_align]
	posScores = [scoreAlignment(a,b, optmatio, gap_start, gap_extend) for a,b in true_pos_align]
	roc.makeRocPlot_givenScores(posScores, negScores, gap_start, gap_extend, "Optmatio")
	#do one with rescoring
	roc.makeRocPlotLC(optmatio, gap_start, gap_extend, "optmatio_rescored")
