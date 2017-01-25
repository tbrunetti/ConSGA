import numpy as np
from random import randrange
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def main():
	
#--------------------------------INITIALIZE POPULATION------------------------------------------
	print 'INITIALIZING POPULATION...'

	#input the number of generations the GA will run
	generations=3000

	global arrayLength
	arrayLength=100

	global kmerLength
	kmerLength=15
	
	#population of "choromosomes"; 1 chrom=a 10 sequence in this case
	#reads in file out put from detect_novel_tfs.py
	initPopChrom=[]
	with open('initial_pop.txt') as startingInput:
		for line in startingInput:
			data=line.split()
			initPopChrom.append(data)

	#selects first 500 10-mers from random initPopChrom
	totalPopulation=[initPopChrom[x][0] for x in range(len(initPopChrom))]
	print 'The size of the initial population is '+ str(len(totalPopulation))+'\n'

	#index maps to seq number (len should be=initial pop)
	#number stored maps to index of first base in k-mer
	global storeAlignPos
	storeAlignPos=[]
	
	#stores alignment of all 500 k-mers
	global alignmentArray
	alignmentArray={}

	#uses alignmentArray and storeAlignPos output from def successiveAlign()
	global bestFitAtEachGen
	bestFitAtEachGen=[]

	#stores the best summation of entropy at each generation of the GA
	#index is generation number, value is lowest summation entropy
	global bestEntropySumOfEachGen
	bestEntropySumOfEachGen=[]


#------------------------SETTING UP INITIAL ALIGNMENT OF ALL 500 k-mers---------------------------------
	def initialAlign(arrayLength, kmerLength):
		print 'RANDOMIZING INITIAL ALIGNMENT...'
		#This should only be done after the initialization of the population
		#After that, alignment matrix should not be changed except by mutation and/or crossover

		
		for j in range(0, len(totalPopulation)):
			#picks a number to align the first base of a give k-mer
			#number picked maps to column position in alignArray
			#of first base pair in k-mer
			randAlign=randrange(0, (arrayLength-kmerLength))
			storeAlignPos.append(randAlign)
			for i in range(0, kmerLength):
				if randAlign in alignmentArray:
					retrieve=alignmentArray.get(randAlign)
					retrieve.append(totalPopulation[j][i])
					randAlign=randAlign+1
				else:
					alignmentArray[randAlign]=[totalPopulation[j][i]]
					randAlign=randAlign+1

#**Test print statement**
		#print 'The length of the array storeAlignPos after random alignment is '+ str(len(storeAlignPos))+'\n'

		return alignmentArray, storeAlignPos
	
		

#---------------------------SETTING UP ALIGNMENT BASED ON FITNESS-----------------------------------	
	
	def successiveAlign(bestFitness):
		
		#clears out list and dictionary from previous generation
		storeAlignPos=[]
		alignmentArray={}
		maxArrayPlacement=arrayLength-kmerLength
		print "ALIGNING NEXT GENERATION..."	
		for i in range(0, len(totalPopulation)):
			#this chunck of for loop checks for best fit population from previous generation and
			#aligns it in the alignment array as previous
			#if not in the best fit generation, it is randomly assigned a new alignment number
			for j in range(0, len(bestFitness)):
				if i==bestFitness[j][0]:
					randAlign=bestFitness[j][1]
					#print "Sequence found in best fit group...reassessing alignment"
				else:
					randAlign=randrange(0, maxArrayPlacement)		
			storeAlignPos.append(randAlign)

			# same alignment method as in def initialAlign()
			for x in range(0, kmerLength):
				if randAlign in alignmentArray:
					retrieve=alignmentArray.get(randAlign)
					retrieve.append(totalPopulation[i][x])
					randAlign=randAlign+1
				else:
					alignmentArray[randAlign]=[totalPopulation[i][x]]
					randAlign=randAlign+1			
#**Test print statement**
		return alignmentArray, storeAlignPos

#-------------------------CALCULATING FITNESS OF ALIGNMENT-------------------------------------------


	def fitnessCalc(alignmentArray, storeAlignPos, bestFitAtEachGen, bestEntropySumOfEachGen):
		
		#Calculate the entropy of each column in alignmentArray
		print 'CALCULATING FITNESS...'

		#dictionary of total nucleotide counts at each index
		#in alignmentArray
		countsDict={}
		for i in range(0, len(alignmentArray)):
			#initialization and clearing of list
			tempStore=[]
			#returns the list from the index of alignmentArray dictionary
			listFromAlignArray=alignmentArray.get(i)
			if listFromAlignArray==None:
				tempStore.append([0, 0, 0, 0])
			else:
				#counts the number of times a particular nucleotide appears at index i
				#and stores in a list tempStore
				#appends total A's as 0
				tempStore.append(0)
				#appends total T's as 0
				tempStore.append(0)
				tempStore.append(listFromAlignArray.count('G'))
				tempStore.append(listFromAlignArray.count('C'))
				#stores the counts from tempStore from index i of alignmentArray
				countsDict[i]=tempStore

		#probabilities of nucleotide at each index
		probabilityNuc={}
		for i in range(0, len(countsDict)):
			tempStore=[]
			if countsDict.get(i)==None:
				probabilityNuc[i]=[0, 0, 0, 0]
			else:
				totalListSum=float(sum(countsDict.get(i)))
				if totalListSum==0:
					probabilityNuc[i]=[0, 0, 0, 0]
				else:
					for l in range(0, 4):
						listRetrieve=countsDict.get(i)
						tempStore.append(listRetrieve[l]/totalListSum)
					probabilityNuc[i]=tempStore

		#stores the entropy calculation at each index
		entropyCalc=[]
		for j in range(0, len(probabilityNuc)):
			tempEntropy=[]
			tempCalc=0
			getList=probabilityNuc.get(j)
			for i in range(0, 4):
				#condition exists since the log2 of 0 will be undefined, so prevent errors
				if getList[i]!=0:
					tempCalc=tempCalc+(getList[i]*(np.log2(probabilityNuc[j][i])))
				#maximum pentality for having nothing in the position of the sequence
				elif getList[0]==0 and getList[1]==0 and getList[2]==0 and getList[3]==0:
					tempCalc=-2.0
		 	entropyCalc.append(tempCalc*-1)

		#pick the k consecutive positions that lower the sum of the entropies

		#stores the sum of entropies in 10-mer window
		#index of list corresponds to first base in the 10-mer window
		entropySums=[]
		for p in range(0, len(entropyCalc)-10):
			totalSum=0
			for w in range(0, 10):
				totalSum=totalSum+entropyCalc[p+w]
			entropySums.append(totalSum)

		#find the minimum sum of entropy window and returns the entropy value and index of window
		minValue=min(entropySums)
		location=entropySums.index(minValue)
		bestEntropySumOfEachGen.append(minValue)


		print 'The minimum entropy total is '+str(minValue)+' and occurs at index '+str(location)+' of entropySums'
	

		#stores the sequence number and alignment position of sequences
		#within the best 10-mer window of lowest entropy as a tuple
		bestFit=[]
		for i in range(0, len(storeAlignPos)):
			if storeAlignPos[i]<location+10 and storeAlignPos[i]>location-10:
				bestFit.append((i, storeAlignPos[i]))
		#stores the best fit list for each generation with the bestFit tuples
		bestFitAtEachGen.append(bestFit)

		return bestFitAtEachGen, bestEntropySumOfEachGen, location

#------------------------------INDUCE RANDOM POINT MUTATIONS-----------------------------------

	def mutation(arrayLength, kmerLength, bestFitness, storeAlignPos):
		#clears out list and dictionary from previous generation
		
		origPositions=storeAlignPos
		storeAlignPos=[]
		alignmentArray={}
		maxArrayPlacement=(arrayLength-kmerLength)
		#can shift to the right or left up to 3 positions in either direction
		mutationalShift=[-3, -2, -1, 1, 2, 3]
		print "APPLYING MUTATIONS...AND REALIGNING..."	
		for i in range(0, len(totalPopulation)):
			#this chunck of for loop checks for best fit population from previous generation and
			#aligns it in the alignment array as previous
			#if not in the best fit generation, it is randomly assigned a new alignment number
			for j in range(0, len(bestFitness)):
				if i==bestFitness[j][0]:
					randAlign=bestFitness[j][1]
					#print "Sequence found in best fit group...reassessing alignment"
				else:
					#determines if selected for mutation
					willMutate=randrange(0, 2)
					#if equals 1 it will mutate
					if willMutate==1:
						indexSelect=randrange(0, 6)
						shift=mutationalShift[indexSelect]
						randAlign=origPositions[i]+shift
						#check to make sure it will not exceed length of alignment array
						#or does not create a negative alignment number
						if randAlign>maxArrayPlacement:
							randAlign=maxArrayPlacement-(randAlign-maxArrayPlacement)
						elif randAlign<0:
							randAlign=randAlign*(-1)
					#if equals 0 will not mutate 
					else:
						randAlign=origPositions[i]
			storeAlignPos.append(randAlign)

			# same alignment method as in def initialAlign()
			for x in range(0, kmerLength):
				if randAlign in alignmentArray:
					retrieve=alignmentArray.get(randAlign)
					retrieve.append(totalPopulation[i][x])
					randAlign=randAlign+1
				else:
					alignmentArray[randAlign]=[totalPopulation[i][x]]
					randAlign=randAlign+1			
#**Test print statement**
		return alignmentArray, storeAlignPos
			
#-------------------------------------RUNNING THE GA-----------------------------------------
	alignmentArray, storeAlignPos=initialAlign(arrayLength, kmerLength);
	bestFitAtEachGen, bestEntropySumOfEachGen, location=fitnessCalc(alignmentArray, storeAlignPos, bestFitAtEachGen, bestEntropySumOfEachGen);
	bestEnt=min(bestEntropySumOfEachGen)
	bestPos=bestEntropySumOfEachGen.index(bestEnt)
	bestFitness=bestFitAtEachGen[bestPos]
	locOfBestAlignment=location

	for i in range(0, generations):
#**Test print statement**
		#print min(bestEntropySumOfEachGen)	
		alignmentArray, storeAlignPos=mutation(arrayLength, kmerLength, bestFitness, storeAlignPos)
		bestFitAtEachGen, bestEntropySumOfEachGen, location=fitnessCalc(alignmentArray, storeAlignPos, bestFitAtEachGen, bestEntropySumOfEachGen);
		if min(bestEntropySumOfEachGen)<bestEnt:
			bestEnt=min(bestEntropySumOfEachGen)
			bestPos=bestEntropySumOfEachGen.index(bestEnt)
			bestFitness=bestFitAtEachGen[bestPos]
			locOfBestAlignment=location

	print "The lowest entropy is "+str(bestEnt)+"with location of " +str(locOfBestAlignment)+" and the fitness alignment is "+str(bestFitness)
	
	#print "The best Entropy Of Sum At Each generation is"+str(bestEntropySumOfEachGen)
#------------------------------------ANALYSIS AND PLOTS-------------------------------------
	#plots for best sum of entropy at each generation with a new
	#random alignment at each generation, no selection
	#run under def main()

	#prints the sequences that are in the final best fit alignment with the lowest sum of entropy
	#for x in range(0, len(bestFitness)):
	#	sequenceID=bestFitness[x][0]
	#	print totalPopulation[sequenceID]
	
	#stores the final best alignment
	finalAlignment={}


	for x in range(0, len(bestFitness)):
		align=bestFitness[x][1]
		seqID=bestFitness[x][0]
		for z in range(0, kmerLength):
			if align in finalAlignment:
				retrieval=finalAlignment.get(align)
				retrieval.append(totalPopulation[seqID][z])
				align=align+1
			else:
				finalAlignment[align]=[totalPopulation[seqID][z]]
				align=align+1	
	
	#stores total nucleotide counts for best alignment
	countsBestFinal={}
	for i in range(min(finalAlignment), (max(finalAlignment)+1)):
		#initialization and clearing of list
		tempStore=[]
		if finalAlignment.get(i)==None:
			tempStore=[0, 0, 0, 0]
		else:
			#returns the list from the index of alignmentArray dictionary
			listFromFinalArray=finalAlignment.get(i)
			#counts the number of times a particular nucleotide appears at index i
			#and stores in a list tempStore
			#make A's zero to mask
			tempStore.append(0)
			#makes T's zero to mask
			tempStore.append(0)
			tempStore.append(listFromFinalArray.count('G'))
			tempStore.append(listFromFinalArray.count('C'))
		#stores the counts from tempStore from index i of alignmentArray
		countsBestFinal[i]=tempStore
		
	#probabilities of nucleotide at each index
	probabilityBestNuc={}
	finalString=""
	for i in range(min(finalAlignment), (max(finalAlignment)+1)):
		tempStore=[]
		totalListSum=float(sum(countsBestFinal.get(i)))
		if totalListSum==0:
			tempStore.append(0)
		else:
			for l in range(0, 4):
				listRetrieve=countsBestFinal.get(i)
				tempStore.append(listRetrieve[l]/totalListSum)
		probabilityBestNuc[i]=tempStore
		mostLikelyNuc=probabilityBestNuc[i].index(max(probabilityBestNuc[i]))
		#if there is no G or C a dash will be used as a placeholder
		if mostLikelyNuc==0:
			finalString=finalString+'-' 
		elif mostLikelyNuc==1:
			finalString=finalString+'-'
		elif mostLikelyNuc==2:
			finalString=finalString+'G'
		elif mostLikelyNuc==3:
			finalString=finalString+'C'
	print probabilityBestNuc
	print finalString

	
if __name__ == '__main__':
	main()
