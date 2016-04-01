import csv
import numpy as np

def main():
	
	#specifies length of k-mer to make from data
	kmerLength=15

	#when breaking down sequences into k-mers specifies how much overlap wanted for one k-mer to the next
	kmerOverlap=kmerLength
	
	#arrays store data read from csv
	writtenData=[]

	#opens and appends data into appropirate array
	with open('Furlong_Mef2_Targets_final.csv') as inputData:
		for line in inputData:
			data=line.split(',')
			writtenData.append(data)


	#column 5 assuming there is a column zero in the dataset
	#store 410 bp region with MEF2 embedded in region 201-211
	sequencesOnly=[writtenData[x][5] for x in range(len(writtenData))]
	print len(sequencesOnly)
	
	#column 7 assuming there is a column zero in the dataset
	#stores the secific mef2 sequence in 410 bp region testing
	mef2SitesSeq=[writtenData[x][7] for x in range(len(writtenData))]
	print len(mef2SitesSeq)
	
	#lists that will hold sequences that are of a certain length
	#greater than or equal to 410 bp in length
	sequencesOnlyTrimmed=[]
	mef2SitesSeqTrimmed=[]

	for w in range(0, len(sequencesOnly)):
		if len(sequencesOnly[w])>=410:
			sequencesOnlyTrimmed.append(sequencesOnly[w])
			mef2SitesSeqTrimmed.append(mef2SitesSeq[w])

	print "The length of sequencesOnlyTrimmed is "+str(len(sequencesOnlyTrimmed))
	print "The length of mef2SitesSeq is "+str(len(mef2SitesSeqTrimmed))

	#outputs a new file with the trimmed MEF2 sequences from Furlong and the associated MEF2 sites
	createNewFile=open('trimmedSequences.txt', 'w')
	mef2Seqs=open('mef2Seqs.txt', 'w')
	for i in range(0, len(sequencesOnlyTrimmed)):
		createNewFile.write(str(sequencesOnlyTrimmed[i])+'\n')
	
	for j in range(0, len(mef2SitesSeqTrimmed)):	
		mef2Seqs.write(str(mef2SitesSeqTrimmed[j])+'\n')


	#check original sequence lengths before trimmingor 
	#for j in range(0, len(writtenData)):
	#	print len(writtenData[j][5])
	
	#counts the total number of each nucleotide in each sequence
	totNucPerString={}

	for x in range(0, len(sequencesOnlyTrimmed)):
		#temporary array the keeps track of total nucs per sequence
		nucCounts=[]
		adenine=0
		thymine=0
		cytosine=0
		guanine=0
		for i in range(0, 410):
			if sequencesOnlyTrimmed[x][i]=='A':
				adenine=adenine+1
			elif sequencesOnlyTrimmed[x][i]=='T':
				thymine=thymine+1
			elif sequencesOnlyTrimmed[x][i]=='C':
				cytosine=cytosine+1
			elif sequencesOnlyTrimmed[x][i]=='G':
				guanine=guanine+1
		nucCounts.append(adenine)
		nucCounts.append(thymine)
		nucCounts.append(cytosine)
		nucCounts.append(guanine)
		totNucPerString[x]=nucCounts

	#total number of each nucleotide at each position of total sequences
	totalNucPerPosition={}

	#probability of each nucleotide appearing at each position of a sequence
	totalNucProbability={}
	
	for n in range(0, 410):
		#temporary counts of nucs per position in sequence
		positionCounts=[]
		adenine=0
		thymine=0
		cytosine=0
		guanine=0
		for l in range(0, len(sequencesOnlyTrimmed)):
			if sequencesOnlyTrimmed[l][n]=='A':
				adenine=adenine+1
			elif sequencesOnlyTrimmed[l][n]=='T':
				thymine=thymine+1
			elif sequencesOnlyTrimmed[l][n]=='C':
				cytosine=cytosine+1
			elif sequencesOnlyTrimmed[l][n]=='G':
				guanine=guanine+1
		positionCounts.append(adenine)
		positionCounts.append(thymine)
		positionCounts.append(cytosine)
		positionCounts.append(guanine)
		totalNucPerPosition[n]=positionCounts
		#clears out position counts so can be used for probability storage
		positionCounts=[]
		positionCounts.append(adenine/float(len(sequencesOnlyTrimmed)))
		positionCounts.append(thymine/float(len(sequencesOnlyTrimmed)))
		positionCounts.append(cytosine/float(len(sequencesOnlyTrimmed)))
		positionCounts.append(guanine/float(len(sequencesOnlyTrimmed)))
		totalNucProbability[n]=positionCounts
		
	#just a check to make sure it works
	#print "The probabilities associated with postion 5 is "+str(totalNucProbability[5])
	#print "The probability of finding a guanine at position 5 is "+str(totalNucProbability[5][3])
	
	#Entropy calculation of each position
	positionEntropy=[]

	for q in range(0, len(sequencesOnlyTrimmed)):
		tempEntropy=0
		for e in range(0, 4):
			if totalNucProbability[q][e]==0:
				tempEntropy=tempEntropy
			else:
				tempEntropy=tempEntropy+totalNucProbability[q][e]*np.log2(totalNucProbability[q][e])
		positionEntropy.append(tempEntropy*-1)

	print "The lowest entropy is "+str(min(positionEntropy))+" and the highest entropy is "+str(max(positionEntropy))

	#print out index and entropy of any value below 1.5
	#also print total nucleotide counts in those positions
	for g in range(0, len(positionEntropy)):
		if positionEntropy[g]<1.5:
			print positionEntropy[g]
			print g
			print totalNucPerPosition[g]

#-----------------------------------------------------------------------------------------------------------------------------#
#break down a seqeunce into k-mers
	
	kmerStore=[]
	#creates k-mers of length kmerLength nucleotides long
	f=open('initial_pop.txt', 'w')
	for j in range(0, len(sequencesOnlyTrimmed)):
	#required do does not override kmer length when iterating through sequences making kmers	
		counter=300+kmerOverlap
		#410 refers to the lenth of each sequence in the data set
		#range here 150-250 covers mef2 sequence
		for i in range(300, 405, kmerOverlap):
			kmerStore.append(sequencesOnlyTrimmed[j][i:counter])
			f.write(str(sequencesOnlyTrimmed[j][i:counter])+'\n')
			counter=counter+kmerOverlap
	
	print len(kmerStore) 

if __name__=='__main__':
	main();