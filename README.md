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

These five concepts can be applied programmatically to many scenarios in order to solve complex problems.  Some examples of where genetic algorithms have been very successful are:
  *  Software bug detection and repair ([Forrest et al 2009](https://www.cs.cmu.edu/~clegoues/docs/legoues-gecco09.pdf), [Le Goues et al 2013](https://ieeexplore.ieee.org/document/6035728))
  *  
  
## What is the theory behind ConSeGA and how does it work?
Chromosomes:  A set of strings, ACTG, 15 characters in length derived from the MEF2-ChIP data 

Fitness Function:  The window with the sum of the 10 lowest consecutive entropy values

Mutation:  String shift amount in alignment

Evolve for Generations:  Elitism for 3000 iterations of the algorithm


## Software Dependencies
ConSeGA was written in Python version 2.7.  The following Python dependecies are required:
* NumPy
* Matplotlib
