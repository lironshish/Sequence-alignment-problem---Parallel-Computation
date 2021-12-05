# Sequence alignment problem - Parallel Computation

### The project deals with a comparison between series of letters.
### The goal is to find similarities between biological molecules, each letter represents a chemical entity and a series of letters represents the structure of a molecule.
### The project solves the problem of series comparison and is called in biology "sequence alignment".
### To calculate the alignment score we compare two series:
* ### If the letters are the same we will mark the pair with $
* ### If the letters belong to a first type group we will mark the pair with %
* ### If the letters belong to a second type group we will mark the pair with the sign #
* ### Else we will mark ' '

## The groups:
* ### First type group:
  ### NDEQ NEQK STA
  ### MILV QHRK NHQK
  ### FYW HY MILF

* ### Second type group:
  ### SAG ATV CSA
  ### SGND STPA STNK
  ### NEQHRK NDEQHK SNDEQK
  ### HFY FVLIM

## The formula for calculating the score:
### alignment score = W1 * numberOfDollars - W2 * NumberOfPercents - W3 * NumberOfHashes - W4 * NumberOfSpaces 
* W1, W2, W3, W4 are non-negative integers that will appear in the input 

## Offset and Mutant Sequence
### Write Seq2 below Seq1 in a possible offset when offset >= 0. Seq2 letters must not appear beyond the end of Seq1.


## Problem Description:
### For given strings Seq1, Seq2 where Seq2 is shorter, we will look for the offset and the MS(k) mutation of Seq2 for which it will be accepted the maximum alignment score when comparing the Seq2 mutation to Seq1.
### We will define the Mutant Sequence to be marked in MS(k) as the series of letters obtained by adding a hyphen after the k sign in seq when k = 1, 2,… (strlen (seq)).

## Utils
### The project was written using MPI, OpenMP and Cuda in Linux environment

