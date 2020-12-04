# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:47:04 2020

@author: Muhammad Ayman Ezzat
"""

import numpy as np

def global_allignment(seq1,seq2,m,mm,g):
    """ Implements the Needleman-Wunsch global allignment algorithm.
    Args:
        seq1(str):The first sequence [Case insensitive]
        seq2(str):The second sequence [Case insensitive]
        m(float):The score given on a sequence match occurrence
        mm(float):The score given on a sequence mismatch occurrence
        g(float):The score given on a gap placement occurrence
        
    Output:
       1- Final Allginment score based on the provided scoring system.
       2- The resulting two sequences after allgining the initial two sequences.
       3- The Calculated Needleman-Wunsch algorithm matrix.
       
    Returns:
        matrix[len(seq1)-1][len(seq2)-1] : Final Allginment score based on the provided scoring system.
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
    
    #Matrix Set-UP
    for i in range(len(seq1)) :
        matrix[i][0] = i*gap_score
    for i in range(len(seq2)) :
        matrix[0][i] = i*gap_score

    #Matrix Calculation
    for i in range(1,len(seq1)) :
        for j in range(1,len(seq2)) :
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
            
    #Traceback
    directions_matrix = np.reshape(list(directions),(len(seq1)-1,len(seq2)-1))
    directions_matrix = np.vstack([["_"]*directions_matrix.shape[1],directions_matrix])
    directions_matrix = np.column_stack([["_"]*directions_matrix.shape[0],directions_matrix])
    allignment1=[]
    allignment2=[]
    i = len(seq1)-1
    j = len(seq2)-1
    
    
    while(1):
        if(directions_matrix[i][j] == "D"):
            allignment1.append(seq1[i])
            allignment2.append(seq2[j])
            i-=1
            j-=1
        elif(directions_matrix[i][j] == "C"):
            allignment1.append(seq1[i])
            allignment2.append("*")
            i-=1
        else:
            allignment2.append(seq2[j])
            allignment1.append("*")
            j-=1
            
        if(i < 0):
            break
        if(j < 0):
            break



    #Printing Results

    allignment1.reverse()
    allignment2.reverse()
    
    print("Final allignment score : ",matrix[len(seq1)-1][len(seq2)-1],"\n")
    print("Allgined Sequences : ")
    print(allignment1[1:])
    print(allignment2[1:])
    print("Note that * indicates a gap\n")
    print("Matrix : ")
    print(matrix)
    
    return matrix[len(seq1)-1][len(seq2)-1]
    

#Program

print("Global Sequence Alginment !")
print("Please enter the two sequences to allign [Case insensitive] !")
S1 = input("First Sequence : ")
S2 = input("Second Sequence : ")
print("Please provide your scoring weights")
M = input("Matching Score : ")
MM = input("Mismatching Score : ")
G = input("Gap Score : ")
print("\nGlobal Allginment : ")
global_allignment(S1,S2,M,MM,G)



