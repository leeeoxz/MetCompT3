from Bio import Entrez, SeqIO
import random 

Entrez.email = "leonardoas95@gmail.com"

handle = Entrez.esearch(db="nucleotide", term = '"16S ribosomal RNA gene"[Tilt] NOT "partial sequence" [Tilt] NOT "uncultured" [Tilt] NOT "clone" [Tilt]', retmax = 20000)
record = Entrez.read(handle)
n = [int(round(random.uniform(0,9000))) for x in xrange(10)]
fast = []
for i in n:
    handle = Entrez.efetch(db="nucleotide",id=record["IdList"][i], rettype = "fasta")
    item = handle.read()
    fast.append(item)

f = open("res1.fasta","w")

for item in fast:
    f.write(item)
    f.write("\n")

"""
for DB in record[u'DbList']:
    handle = Entrez.einfo(db=DB)
    r1 = Entrez.read(handle)
    print "%s - %s - %s"%(DB,r1[u'DbInfo'][u'Description'],r1[u'DbInfo'][u'Count'])
    print "___________"
"""


