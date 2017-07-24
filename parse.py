from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

pars = open("res.fasta","r")
tag = '16S'
nomes = []
info =[] #LIST WITH LISTS OF INFORMATIONS, EACH ITEM OF THIS LIST IS [ID,NAME,SEQ,SEQ'S LEN]
for sec in SeqIO.parse(pars, "fasta"):
    item = []
    item.append(sec.id)
    n = sec.description.split() 
    for x,i in enumerate(n):
        if i == tag:
            index = x    
    name = " ".join(n[:index])
    nomes.append(name)
    item.append(name)
    item.append(sec.seq)
    item.append(len(sec))
    info.append(item)
for item in info:
	print item

#for a in pairwise2.align.globalms("ACCGT", "ACG",5,-4,-10,-1):
#    print(format_alignment(*a))