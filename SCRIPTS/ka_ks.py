from Bio import SeqIO
import os
import re
import sys
import subprocess


# ================================================================================================
# Parse protein sequences from fasta file (PEP)
# ================================================================================================
filinPEP = sys.argv[1]
dataPEP = SeqIO.parse(filinPEP,"fasta")
PEP = {}

for record in dataPEP :
	# Retrieve gene name
	m = re.search(r'gene:(\S*)\s+transcript:(\S*)\s+', record.description)
	gene = m.group(1)
	isoform = m.group(2)

	#print("gene : " + gene + ", isoform : " + isoform)

	l = len(record.seq)
	# If gene is not in dictionnary ==> store it
	if gene not in PEP:
		PEP[isoform] = [record, gene, l]

	# If gene already in dictionnay ==> keep the record with the longest sequence
	else:
		if l > PEP[isoform][2]:
			PEP[isoform] = [record, gene, l]


# ================================================================================================
# Parse DNA sequences from fasta file (CDS)
# ================================================================================================
filinCDS = sys.argv[2]
dataCDS = SeqIO.parse(filinCDS,"fasta")
CDS = {}

for record in dataCDS :
	# Retrieve gene name
	m = re.search(r'gene:(\S*)\s+', record.description)
	gene = m.group(1)

	#print("gene : " + gene)

	l = len(record.seq)
	# If gene is not in dictionnary ==> store it
	if gene not in CDS:
		CDS[gene] = [record, l]

	# If gene already in dictionnay ==> keep the record with the longest sequence
	else:
		if l > CDS[gene][1]:
			CDS[gene] = [record, l]


# ================================================================================================
# Files preparation
# ================================================================================================
yn00 = os.getcwd()+"/tmp_cds.ali.phy"
subprocess.run("awk -v file="+yn00+" '{gsub(\"XXXXX\",file); print $0}' SCRIPTS/tmp_yn00.ctl > tmp_yn00.ctl", shell = True)
filout = open("RESULTS/paml_result.txt", "w")
filout.write("family\tpep1\tpep2\tomega\n")

# ================================================================================================
# Create couple from same family and calculate Ka/ks
# ================================================================================================
filinFam = open(sys.argv[3], 'r')
filinFam.readline()
i = 0
li = []

for gene in filinFam :
	name = gene.split()[0]
	family = gene.split()[1]

	# If different family ==> create new list
	if family != i :
		li = []
		i = family

		###print("================================================================================")

	# If same family ==> calcultate ka/ks with other genes already in list
	else :
		for gen in li :
			
			# Create temporary files (PEP and CDS)
			SeqIO.write([PEP[name][0], PEP[gen][0]], "tmp_pep.fasta", "fasta")
			SeqIO.write([CDS[PEP[name][1]][0], CDS[PEP[gen][1]][0]], "tmp_cds.fasta", "fasta")

			# Alignment of PEP (clustal)
			subprocess.call("clustalw -quiet -align -infile=tmp_pep.fasta -outfile=tmp_pep.ali.aln", shell=True, stdout=open(os.devnull, 'wb'))
			
			# ALignment of CDS (pal2nal)
			subprocess.run("./SCRIPTS/pal2nal.pl tmp_pep.ali.aln tmp_cds.fasta -output paml > tmp_cds.ali.phy", shell=True)
			
			# Calculation of Ka/Ks (PAML)
			subprocess.call("yn00 tmp_yn00.ctl", shell=True, stdout=open(os.devnull, 'wb'))

			# Parse yn file to retrieve Ka/KS value
			values = str(subprocess.check_output("awk '/seq. seq./{flag=1;next}/LWLm methods/{flag=0}flag' yn | awk 'NF > 0'", shell=True))
			values = re.findall(r'-?\d+\.?\d*', values)

			# Write result in file
			if len(values) < 9 :
				filout.write(family+"\t"+name+"\t"+gen+"\t"+"NA"+"\t"+"NA"+"\t"+"NA"+"\n")
			else :
				filout.write(family+"\t"+name+"\t"+gen+"\t"+values[6]+"\t"+values[7]+"\t"+values[9]+"\n")

	# Add to list
	li.append(name)


# Remove temporary files
subprocess.run("rm 2YN.dN 2YN.t rst1 tmp tmp_cds.fasta tmp_pep.dnd tmp_yn00.ctl yn00.ctl 2YN.dS rst rub tmp_cds.ali.phy tmp_pep.ali.aln tmp_pep.fasta yn", shell = True)