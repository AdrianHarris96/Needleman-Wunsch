# Needleman-Wunsch

The python script should obtain the optimal global alignment for a pair of sequences. 

This stems from a project in my introductory bioinformatics course and serves an introduction to dynamic programming. 

The larger objective is divided into subproblems: initialization of the matrix, matrix filling via the scoring system, and backtracking to determine the optimal path. 

The backtracking rules for this program:

1) Always take the diagonal when the diagonal is either (1) the highest score or (2) tied for the highest score. 

2) If the diagonal is not the highest score, take the up if it is either (1) the highest score or (2) tied for the highest score. 

3) Take the left if the top two conditions do not hold true. 

Sample Input: ./aharris334_nwAlign.py <seq1_nw.fa> <seq2_nw.fa>