# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:33:26 2020

@author: Muhammad Ezzat
"""

flag = input("Welcome to the sequence alignment program!\nFor global alignment enter 0\nFor local alignment enter 1\n")
if(flag == 0):
    runfile("global.py")
else:
    runfile("local.py")