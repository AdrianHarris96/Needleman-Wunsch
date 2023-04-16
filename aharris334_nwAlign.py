#!/usr/bin/env python3

import sys

#Example Input: ./<gtusername>_nwAlign.py seq1_nw.fa seq2_nw.fa
#Code contains comments for print debugging. 

filename1 = sys.argv[1]
filename2 = sys.argv[2]

#Simply open each fasta file and add seq to empty string.
file1 = open(filename1, "r")
seq1 = ""
seq2 = ""
for line in file1:
  if line.startswith(">"):
    pass
  else:
    line = line.strip()
    seq1 += line

file2 = open(filename2, "r")
for line in file2:
  if line.startswith(">"):
    pass
  else:
    line = line.strip()
    seq2 += line

#Code does not work if the shorter sequence is used in vertical orientation. This ensures that does not happen. 
if len(seq2) < len(seq1):
  tempseq = seq1
  seq1 = seq2
  seq2 = tempseq
else:
  pass


#Scoring system
match = 1
mismatch = -1
gap = -1

#Creation of empty matrix.
#Need this if I plan to move via loops row-by-row, populating the matrix later. 
lengthSeq1 = len(seq1)+1
lengthSeq2 = len(seq2)+1
matrix = []
for x in range(lengthSeq2):
  newList = []
  count = 0
  for y in range(lengthSeq1):
    if x == 0:
      item = y * gap
      newList.append(item)
    else:
      if count == 0:
        item = x * gap
        newList.append(item)
        count += 1
      else:
        item = 0
        newList.append(item)
  matrix.append(newList)

'''
Populating the empty matrix from earlier.
We first check horizontal sequence with the index [y-1] is equal to the vertical sequence with index[x-1]. Confusing, I know!
If so, the diag value is the diagonal value + match value. If not, we say it is diagonal value + mismatch. 
Up and left are always just gaps, so each time it iterates, it sets them to the value above or to the left plus the gap value, respectively. 
Then, we fill up our old matrix with values and our newMatrix with "d", "u" or "l" with that priority using relational operators. 
'''
newMatrix = []
for x in range(1, lengthSeq2):
  newList2 = []
  for y in range(1, lengthSeq1):
    if seq1[y-1] == seq2[x-1]:
      diag = matrix[x-1][y-1] + match
    else:
      diag = matrix[x-1][y-1] + mismatch
    up = matrix[x-1][y] + gap
    left = matrix[x][y-1] + gap
    if diag >= up and diag >= left:
      newList2.append("d")
      matrix[x][y] = diag
    elif up >= left:
      newList2.append("u")
      matrix[x][y] = up
    else:
      newList2.append("l")
      matrix[x][y] = left
  newMatrix.append(newList2)

#for row in matrix:
  #print(row)

#for row in newMatrix:
  #print(row)

count = 0
for row in newMatrix:
  for item in row:
    count +=1

#rowLength (i.e. number of columns) and length of newMatrix (i.e. number of rows).
#-1 is done because we will use the range function in reverse, and we have to pay attention to indexing.
rowLength = int(count / len(newMatrix))
position1 = len(newMatrix) - 1
position2 = rowLength - 1
newString = str(newMatrix[position1][position2])

'''
Breaks are used so we do not continually iterate through the range functions. May not be necessary, considering I am continually setting the position1 and position2.
newString is set from the start with whatever letter is in that corner of the matrix ("d", "u" or "l").
Depending on what it is we will move diagonal, up or left by adjusting position1 and position2.
newString is continually populated with each move. 
'''
for x in range((position1), -1, -1):
  for y in range((position2), -1, -1):
    if newMatrix[position1][position2] == "d":
      position1 -= 1
      position2 -= 1
      newString += newMatrix[position1][position2] 
      break
    elif newMatrix[position1][position2] == "u":
      position1 -= 1
      newString += newMatrix[position1][position2] 
      break
    else:
      position2 -= 1
      newString += newMatrix[position1][position2] 
      break

#print(newString)

#Not sure why, but my code falls apart without this. It could be the way I initially set up the matrices. 
#It prevents it from ending on a "u" or "l".
if (len(seq1)+len(seq2)) % 2 == 0:
  newString += "d"

#Converting string to list
seq1List = []
seq2List = []
for x in range(len(seq1)):
  seq1List.append(seq1[x])
for x in range(len(seq2)):
  seq2List.append(seq2[x])

# Once both even and old combinations of sequences are made equal in terms of their combined length, this length variable above ensures no index error with the next loop. 
length = len(newString)-1

#Simply iterate through string, pop for both if diagonal or one element depending on if it is "u" or "l".
#Using the insert() was messy, so I thought to populate two new list concurrently. 
count = 0
finalSeq1List = []
finalSeq2List = []
for x in range(length):
  if newString[x] == "d":
    finalSeq1List.append(seq1List.pop())
    finalSeq2List.append(seq2List.pop())
    pass
  elif newString[x] == "u":
    finalSeq1List.append("-")
    finalSeq2List.append(seq2List.pop())
  else:
    finalSeq1List.append(seq1List.pop())
    finalSeq2List.append("-")

#print(finalSeq1List)
#print(finalSeq2List)

#The list should be the same length after the "pop"-loop, so I iterate using the length of one, and create another list, telling me if there is a match, mismatch or gap. 
#We also create the score here.
matchList = []
score = 0

for x in range(len(finalSeq1List)):
  if finalSeq1List[x] == finalSeq2List[x]:
    matchList.append("|")
    score += match
  elif finalSeq1List[x] == "-" or finalSeq2List[x] == "-":
    matchList.append(" ")
    score += gap
  else:
    matchList.append("*")
    score += mismatch

#Conversion back to strings, but remember we created newString earlier via backtracking. 
#I have intentionally left these strings backwards.
seq1String = ("".join(str(char) for char in finalSeq1List))
matchString = ("".join(str(char) for char in matchList))
seq2String = ("".join(str(char) for char in finalSeq2List))

print(seq1String[::-1])
print(matchString[::-1])
print(seq2String[::-1])

print("Alignment score:", score)

#Reverse everything to the correct orientation and print
