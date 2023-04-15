# Needleman-Wunsch

## Description
This stems from a project in my introductory bioinformatics course and serves an introduction to dynamic programming. 

This python script should obtain the optimal global alignment for a pair of sequences. 

## DP Steps
(1)`Initialization of the matrix`
(2) `Filling matrix via the scoring system`
(3) `Backtracking to determine the optimal path`

### Backtracking Rules:

1) Always take the diagonal when the diagonal is either (1) the highest score or (2) tied for the highest score. 

2) If the diagonal is not the highest score, take the up if it is either (1) the highest score or (2) tied for the highest score. 

3) Take the left if the top two conditions do not hold true. 

## Example Inputs
Sample Input: ./aharris334_nwAlign.py <seq1_nw.fa> <seq2_nw.fa>