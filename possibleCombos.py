
#lists all possible combinations of nucleotide of length specified by motifLength
import itertools
def possibleCombos():
	motifLength=6

	cartProduct=list(itertools.product('ACTG', repeat=motifLength))
	

	combos=[]
	for i in range(0, len(cartProduct)):
		cleanedString=""
		for j in range(0, motifLength):
			cleanedString=cleanedString+str(cartProduct[i][j])
		combos.append(cleanedString)
	print 'The total different combination of motifs of length '+ str(motifLength) + " is " + str(len(combos))

	#writes combinations in an output file
	f=open('allCombinationMotifs.txt', 'w')
	for j in range(0, len(combos)):
		f.write(combos[j]+'\n')

if __name__=='__main__':
	possibleCombos();