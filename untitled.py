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
			matrix[z+1+x][x] = 'EQ'#matrix[x][z+1+x]
for x,item in enumerate(sequences):
	for z,j in enumerate(item):
		if x == z:
			matrix[x][z] = 'EQ'
mt = copy.deepcopy(matrix) #CREATES A COPY TO CALCULATE THE DISTANCE MATRIX

for x,item in enumerate(mt):
	for y,score in enumerate(item):
		if type(score) is int:
			if score != 0:
				mt[x][y] = round(1/score,5) #CALCULATES THE DISTANCE MATRIX

#mt =[['EQ',3,3,4,5],['EQ','EQ',2,2,1],['EQ','EQ','EQ',-2,-1],['EQ','EQ','EQ','EQ',-3],['EQ','EQ','EQ','EQ','EQ']]

def getMin(mt,t):
	for x,item in enumerate(mt):
		if x == 0:
			minimo = min(item)
			indexy = item.index(minimo)
			indexx = x
		else:
			if minimo > min(item[x:]):
				if (x,item.index(min(item[x:]))) not in t:
					minimo = min(item[x:])
					indexy = item.index(minimo)
					indexx = x
	t.append((indexx,indexy))
	return indexx,indexy,minimo,tree

n=0
tree = []
while n != size:
	r = getMin(mt,tree)
	mt[r[0]][r[1]] = 'EQ'
	tree = r[3]
	print tree
	n+=1

def kruskal(matx):
	pass