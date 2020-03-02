# ConSeGA - <ins>Con</ins>sensus <ins>Se</ins>qeunce-based <ins>G</ins>enetic <ins>A</ins>lgorithm
---------------------------------------------------------------------------------------------------
A genetic algorithm used to identify regions of regularity in the genome of Drosophila melanogaster.  This algorithm has been tested on previously published ChIP data for the transcription factor, Myocyte-enhancer factor 2 (MEF2).

1.  [Introduction and Theory](#introduction)

2.  [How does it work?](#what-is-the-theory-behind-ConSeGA-and-how-does-it-work)

3.  [Software Dependencies](#software-dependencies)


## Introduction
Genetic algorithms (GA) are algorithms in which the basic principles behind Darwinian evolution are applied algorithmically to a computational probelm of interest.  This approach can generate results by simply applying the idea of "survival of the fittest". 

   1. **_Chromosomes_** store genetic information
   2.  Individuals in a population that are deemed **_“most genetically fit” prevail_**
   3.  Those that do not, either die off or undergo **_mutation/crossover_**
   4.  Those that are the most **_elite_** continue their lineage in the next generation
   5.  This process is repeated for several generations, until the population **_evolves and equilibrium is reached_**

These five concepts can be applied programmatically to many scenarios in order to solve complex problems.  Some examples of where genetic algorithms have been very successful are (but not limited to):
  *  Software bug detection and repair ([Forrest et al 2009](https://www.cs.cmu.edu/~clegoues/docs/legoues-gecco09.pdf), [Le Goues et al 2013](https://ieeexplore.ieee.org/document/6035728))
  *  Medicine and Biology ([Ghaheri et al 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4678452/), [Hackenberger, 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6509630/), )
  *  Bioinformatics & data science ([Manning et al 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3813526/), [Yang and Honavar](https://link.springer.com/chapter/10.1007/978-1-4615-5725-8_8))
  *  Film, movies, and gaming ([Trescak et al 2012](https://dl.acm.org/doi/abs/10.1145/2407336.2407338?casa_token=HRSGhuKhdp8AAAAA:kJ8C1wpn3-MltHUemdlQHzojyoqLPzeqa_W1wfd6OluF-CG8L_5OOZR5hpq7VqCnzz3Qw-JWfmV8QQ), [Sanjuan et al 2007](https://link.springer.com/chapter/10.1007/978-3-540-71805-5_52))
  
## What is the theory behind ConSeGA and how does it work?

ConSeGA is a genetic algorithm approach for detecting consensus sequences in biological data sets.  

<p align="center"><img src="https://github.com/tbrunetti/ConSeGA/blob/master/pictures/overviewGA.png" width="500"></p>



   *  The algorithm defines the set of set of "individuals" or "chromosomes" as: a set of strings of 15 characters in length comprised of only A, T, G, and C.  In the case of the algorithm this was derived from MEF2-ChIP data ().
   *  The method to determine how "fit" an individual is, or fitness fucntion, is determed by the window of string with a sum of the 10 lowest consecutive entopy values.
   *  Any string that is not in the most fit windows is subject to random mutation.  A mutation in this case, is a random position shift in the alignment array of +/-3 base pairs.
   *  Elitism is use to ensure the most fit window of string is automatically incorporated into the next generation of string alignments.
   *  This process is repeated or "evolved" for 3000 generations.


## Software Dependencies
ConSeGA was written in Python version 2.7.  The following Python dependecies are required:
* NumPy
* Matplotlib
