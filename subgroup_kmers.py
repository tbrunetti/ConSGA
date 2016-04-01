#this programs subdivides the kmers generated from detect_novel_tfs.py using 
#initial_pop.txt
import math
def main():
	#pertinent to change the length and overlap to calculate subgroups correctly
	kmerLength=15
	kmerOverlap=0

	storeKmers=[]
	dataFile=open('initial_pop.txt', 'r')
	for line in dataFile:
		storeKmers.append(line)
	
	#remove EOL from string kmer
	storeKmers=[storeKmers[x].strip('\n') for x in range(len(storeKmers))]

	#math.floor will always round a floating point down
	numberOfkmersPerSeq=math.floor(float(410)/(kmerLength-kmerOverlap))


if __name__=='__main__':
	main();