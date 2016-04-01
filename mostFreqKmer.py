#finds the most frequent kmer
#takes input file allCombinationMotifs.txt output from possibleCombos.py


def main():

	combinations=[]
	sequences410Trimmed=[]

	#imports data of all possible 6-mer combinations of sequences
	with open('allCombinationMotifs.txt') as inputData:
		for line in inputData:
			data=line.split()
			combinations.append(data)

	#imports data from Furlong's data 410bp sequence
	with open('trimmedSequences.txt') as sequences:
		for line in sequences:
			data2=line.split()
			sequences410Trimmed.append(data2)

	#clean up procedure so data is usuable
	cleanedCombos=[]
	for i in range(0, len(combinations)):
		cleanedCombos.append(combinations[i][0])

	#clean up procedure so data is usable
	cleanedSeqs=[]
	for j in range(0, len(sequences410Trimmed)):
		cleanedSeqs.append(sequences410Trimmed[j][0])

	#probability of seeing each k-mer assuming nucleotide independence
	probabilityOfKmer=[]
	#probabilities taken from...........
	probAT=float(0.30)
	probCG=float(0.20)
	for x in range(0, len(cleanedCombos)):
		tempProb=float(1.0)
		for z in range(0, len(cleanedCombos[0])):
			if cleanedCombos[x][z]=='A' or cleanedCombos[x][z]=='T':
				tempProb=tempProb*probAT
			elif cleanedCombos[x][z]=='G' or cleanedCombos[x][z]=='C':
				tempProb=tempProb*probCG
		probabilityOfKmer.append(tempProb)		
	print 'The highest probability is '+str(max(probabilityOfKmer))+' and lowest probability is '+str(min(probabilityOfKmer))

	totalDataA=0;
	totalDataT=0;
	totalDataG=0;
	totalDataC=0;
	for i in range(0, len(cleanedSeqs)):
		totalDataA=totalDataA+cleanedSeqs[i].count('A')
		totalDataT=totalDataT+cleanedSeqs[i].count('T')
		totalDataC=totalDataC+cleanedSeqs[i].count('C')
		totalDataG=totalDataG+cleanedSeqs[i].count('G')
	totalDataNucs=totalDataA+totalDataT+totalDataG+totalDataC

	print 'The percentage of A in the data is '+str(totalDataA/float(totalDataNucs))
	print 'The percentage of T in the data is '+str(totalDataT/float(totalDataNucs))
	print 'The percentage of G in the data is '+str(totalDataG/float(totalDataNucs))
	print 'The percentage of C in the data is '+str(totalDataC/float(totalDataNucs))


	#total number of 6-mers from Furlong Dataset
	#found using detect_novel_tfs.py
	total6MersData=float(148225)


	#counts and stores the total number of times a combination occurs in the sequences
	occurencesOfMotif=[]
	newFile=open('totalOccurences.txt', 'w')
	file150=open('occurence150.txt', 'w')
	file200=open('occurence200.txt', 'w')
	file250=open('occurence250.txt', 'w')
	file300=open('occurence300.txt' ,'w')
	foldChange2=open('foldGreaterThan2.txt', 'w')
	foldChange3=open('foldGreaterThan3.txt', 'w')
	foldChange4=open('foldGreaterThan4.txt', 'w')

	for x in range(0, len(cleanedCombos)):
		totals=0
		for z in range(0, len(cleanedSeqs)):
			totals=totals+cleanedSeqs[z].count(cleanedCombos[x])
		newFile.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n') 
		occurencesOfMotif.append(totals)
		if totals>=150:
			file150.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')
		if totals>=200:
			file200.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')
		if totals>=250:
			file250.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')
		if totals>=300:
			file300.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')

		if (totals/(total6MersData*probabilityOfKmer[x]))>=2.0:
			foldChange2.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')

		if (totals/(total6MersData*probabilityOfKmer[x]))>=3.0:
			foldChange3.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')

		if (totals/(total6MersData*probabilityOfKmer[x]))>=4.0:
			foldChange4.write(str(cleanedCombos[x])+' , '+str(totals)+' , '+str(probabilityOfKmer[x])+' , '+ str(total6MersData*probabilityOfKmer[x])+' , '+str(totals/(total6MersData*probabilityOfKmer[x]))+'\n')

if __name__=="__main__":
	main();