from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import copy

mismatch = -4 #MISMATCH SCORE
match = 5 #MATCH SCORE
gapOpen = -10 #GAP OPENING SCORE
gapExtend = -1 #GAP EXTEND SCORE
sequences = ['AACAGTG','CGATCATCAGT','CAGTGAC','CATGACTCAG','CAGTCAGTCGA','CATGCATGAC'] #SEQUENCE VECTOR
size = len(sequences) #GETS THE NUMBER OF SEQUENCES TO DETERMINE THE MATRIX'S SIZE
matrix = [[0 for x in xrange(size)] for y in xrange(size)] # CREATES THE MATRIX

for x,item in enumerate(sequences):
	for z,j in enumerate(sequences[x+1:]):
		matrix[x][z+1+x] = pairwise2.align.globalms(item,j,match,mismatch,gapOpen,gapExtend)[0][2]
		matrix[z+1+x][x] = matrix[x][z+1+x]

mt = copy.deepcopy(matrix) #CREATES A COPY TO CALCULATE THE DISTANCE MATRIX

for x,item in enumerate(mt):
	for y,score in enumerate(item):
		if score != 0:
			mt[x][y] = round(1/score,5) #CALCULATES THE DISTANCE MATRIX

for item in mt:
	print item

#SO FAZER O KRUSKAL AGORA #