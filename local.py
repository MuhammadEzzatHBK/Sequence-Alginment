# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:47:04 2020

@author: Muhammad Ayman Ezzat
"""
import numpy as np
def local_alignment(seq1,seq2,m,mm,g):
    """ Implements the Smith-Waterman local alginment algorithm.
    Args:
        seq1(str):The first sequence [Case insensitive]
        seq2(str):The second sequence [Case insensitive]
        m(float):The score given on a sequence match occurrence
        mm(float):The score given on a sequence mismatch occurrence
        g(float):The score given on a gap placement occurrence
        
    Output:
       1- The maximum local alignment score based on the provided socring system.
       2- The two best aligned local sequences.
       3- The Calculated Smith-Waterman algorithm matrix.
       
    Returns:
        maximum : Maximum Alginment score based on the provided scoring system.
    """
    #Matrix initialization
    match_score = float(m)
    mismatch_score = float(mm)
    gap_score = float(g)
    seq1 = "_" + seq1
    seq2 = "_" + seq2
    matrix = np.zeros((len(seq1),len(seq2)))
    directions=""
    seq1=seq1.upper()
    seq2=seq2.upper()
    
    #Matrix Calculation
    for i in range(1,len(seq1)):
        for j in range(1,len(seq2)):
            if(seq1[i] == seq2[j]) :
                D = match_score + matrix[i-1][j-1]
            else :
                D = mismatch_score + matrix[i-1][j-1]
            R = gap_score + matrix[i][j-1]
            C = gap_score + matrix[i-1][j]
            key = np.argmax([D,R,C])
            if(key == 0):
                directions +="D"
            elif(key ==1):
                directions +="R"
            else:
                directions +="C"
            matrix[i][j] = np.max([D,R,C])
            if(matrix[i][j] < 0):
                matrix[i][j] = 0
                
                
    #Traceback
    directions_matrix = np.reshape(list(directions),(len(seq1)-1,len(seq2)-1))
    directions_matrix = np.vstack([["_"]*directions_matrix.shape[1],directions_matrix])
    directions_matrix = np.column_stack([["_"]*directions_matrix.shape[0],directions_matrix])
    alignment1=[]
    alignment2=[]
    
    maximum = np.max(matrix)
    for x in range(1,len(seq1)):
        for y in range(1,len(seq2)):
            if(matrix[x][y] == maximum):
                i = x
                j = y
 
    while(1):
        if(directions_matrix[i][j]== "D"):
            alignment1.append(seq1[i])
            alignment2.append(seq2[j])
            i-=1
            j-=1
        elif(directions_matrix[i][j] == "C"):
            alignment1.append(seq1[i])
            alignment2.append("*")
            i-=1
        else:
            alignment2.append(seq2[j])
            alignment1.append("*")
            j-=1
            
        if(i < 0):
            break
        if(j < 0):
            break
        if(matrix[i][j] == 0):
           
            break
    
    alignment1.reverse()
    alignment2.reverse()
    
    print("Maximum alignment score : ",maximum,"\n")
    print("Algined Sequences : ")
    print(alignment1)
    print(alignment2)
    print("Note that * indicates a gap\n")
    print("Matrix : ")
    print(matrix)
    
    return maximum

print("Local Sequence Alginment !")
print("Please enter the two sequences to align [Case insensitive] !")
S1 = input("First Sequence : ")
S2 = input("Second Sequence : ")
print("Please provide your scoring weights")
M = input("Matching Score : ")
MM = input("Mismatching Score : ")
G = input("Gap Score : ")
print("\nLocal Alginment : ")
local_alignment(S1,S2,M,MM,G)