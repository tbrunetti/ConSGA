import re

def sixMersFoldChange():

	#reads in all sequences and stored in seqences array
	sequences=[]
	with open("trimmedSequences.txt") as seq:
		for line in seq:
			singleSeq=line.split()
			sequences.append(singleSeq)

	#cleaned usable sequences
	seqs=[sequences[i][0] for i in range(0, len(sequences))]

	
	#finds all 6mers and then takes 10 nt upstream and downstream from site
	#and saves them in the store array
	sixMer=re.compile('GCTGCC')
	store=[]
	for i in range(0, len(seqs)):
		seqSelction=seqs[i]
		findIt=sixMer.findall(seqSelction)
		iterableF=sixMer.finditer(seqSelction)
		for match in iterableF:
			location=match.span()
			lowerBound=location[0]-10
			upperBound=location[1]+10
			temp=seqSelction[lowerBound:upperBound]
			#only get sequences length 26 to make sure there is a position at every string
			#to determine most conserved nucleotides at each position
			if len(temp)==26:
				store.append(temp)

	#prints out a final string in numbers of the most conserved nucleotide at each position
	oligoString=''
	for x in range(0, len(store[0])):
		adenine=0
		thymine=0
		guanine=0
		cytosine=0
		positionCounts=[]
		for j in range(0, len(store)):
			if store[j][x]=='A':
				adenine=adenine+1;
			elif store[j][x]=='T':
				thymine=thymine+1;
			elif store[j][x]=='G':
				guanine=guanine+1;
			elif store[j][x]=='C':
				cytosine=cytosine+1
		positionCounts.append(adenine)
		positionCounts.append(thymine)
		positionCounts.append(guanine)
		positionCounts.append(cytosine)

		mostConserved=positionCounts.index(max(positionCounts))		
		oligoString=oligoString+str(mostConserved)


	#converts oligoString into nucleotides
	finalOligo=''
	for i in range(0, len(oligoString)):
		if oligoString[i]=='0':
			finalOligo=finalOligo+'A'
		elif oligoString[i]=='1':
			finalOligo=finalOligo+'T'
		elif oligoString[i]=='2':
			finalOligo=finalOligo+'G'
		elif oligoString[i]=='3':
			finalOligo=finalOligo+'C'

	print len(finalOligo)
	print finalOligo

if __name__=='__main__':
	sixMersFoldChange();