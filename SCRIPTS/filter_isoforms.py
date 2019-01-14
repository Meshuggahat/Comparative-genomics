from Bio import SeqIO
import re
import sys


# INPUT FILE
filin = sys.argv[1]
data = SeqIO.parse(filin,"fasta")
records = list(SeqIO.parse(filin,"fasta"))


# INPUT FASTA INFORMATIONS
sizes = [len(rec) for rec in records]
print("INPUT")
print("Nombre de séquences protéiques dans le fichier FASTA :", len(sizes)) 
print("Taille moyenne des protéines : ", sum(sizes)/len(sizes))


# RETRIEVING UNIQUE ISOFORM FOR EACH GENE AND STORE IN DICTIONNARY
res = {}

for record in data :
	# Retrieve gene name and isoform name
	m = re.search(r'gene:(\S*)\s+transcript:(\S*)\s+', record.description)
	gene = m.group(1)
	isoform = m.group(2)

	#print(record.description)
	#print("gene : " + gene + ", isoform : " + isoform)

	l = len(record.seq)
	# If gene is not in dictionnary ==> store it
	if gene not in res:
		res[gene] = [record, l]
	# If gene already in dictionnay ==> keep the record with the longest sequence
	else:
		if l > res[gene][1]:
			res[gene] = [record, l]


# OUTPUT FASTA INFORMATIONS
sizes = [rec[1] for rec in res.values()]
print("\nOUTPUT")
print("Nombre de séquences protéiques restantes : ", len(res.keys()))
print("Taille moyenne des protéines : ", sum(sizes)/len(sizes))


# WRITE FILTERED SEQUENCES IN FASTA FILE
filout = [rec[0] for rec in res.values()]
SeqIO.write(filout, "filtered.fasta", "fasta")